# Assumption stress tests. Source this file for helpers; sourcing NEVER simulates.
# From the project root (or remove Lyu_2018/ when inside that folder):
# Rscript --vanilla Lyu_2018/R/assumption_stress.R 10000 32 2018
# This calls the existing simulate_table3(), run_basyc() and select_mts().

stress_directory <- function() {
  if (file.exists("Lyu_2018/R/mtpi2.R")) "Lyu_2018" else "."
}

# Rows are doses; columns are cycles. b contains the three cycle-1 risks.
# r is the first relative increase; r*c is the second relative increase.
stress_probability <- function(b, r = 1.1, c = 1) {
  stopifnot(length(b) == 3, all(is.finite(b)), all(b > 0 & b <= 1),
            all(diff(b) >= 0), length(r) == 1, is.finite(r), r > 0,
            length(c) == 1, is.finite(c), c > 0)
  t <- 0:2
  pmin(outer(b, r^t * c^(t * (t - 1) / 2)), 1)
}

stress_scenarios <- function() {
  data.frame(scenario = 1:10,
    label = c(paste("Paper S", 1:6, sep = ""),
              "Early dose jump", "Late dose jump", "Reversed S1 cycles", "Peak then recovery"),
    b1 = c(.15,.15,.08,.48,.08,.15,.20,.25,.1815,.20),
    b2 = c(.30,.27,.16,.66,.16,.46,.275,.255,.363,.27),
    b3 = c(.45,.45,.24,.84,.58,.53,.285,.33,.5445,.40),
    r = c(rep(1.1, 7),1.02,10/11,1.3), c = c(rep(1,9),.5),
    purpose = c(rep("Original Table 3", 6),
      "Intensity violation; cycle trend preserved",
      "Opposite dose curvature; cycle trend preserved",
      "Cycle trend violated; intensity ordering preserved",
      "Nonmonotone cycles; global intensity violation"))
}

stress_truth <- function(specifications = stress_scenarios()) {
  lapply(seq_len(nrow(specifications)), function(i) {
    s <- specifications[i, ]
    stress_probability(c(s$b1,s$b2,s$b3),s$r,s$c)
  })
}

# Enumerate true risks independently of the trial's fitted PAVA values.
stress_sequences <- function(p, target = .30, safety = .35) {
  d <- as.matrix(expand.grid(d1=1:3,d2=1:3,d3=1:3))
  q <- apply(d,1,function(x) mean(p[cbind(x,1:3)]))
  safe <- apply(d,1,function(x) all(p[cbind(x,1:3)] < safety))
  down <- apply(d,1,function(x) all(diff(x) <= 0))
  optimum <- function(allowed) {
    if (!any(allowed)) return(rep(FALSE,nrow(d)))
    allowed & abs(abs(q-target)-min(abs(q[allowed]-target))) < 1e-12
  }
  data.frame(sequence=apply(d,1,paste,collapse=","),intensity=rowSums(d),
             q=q,safe=safe,nonincreasing=down,
             optimal_down=optimum(safe & down),optimal_all=optimum(safe))
}

stress_audit <- function(p) {
  z <- stress_sequences(p)
  pairs <- which(outer(z$intensity,z$intensity,"<"),arr.ind=TRUE)
  gap <- z$q[pairs[,1]]-z$q[pairs[,2]]
  safe_pair <- z$safe[pairs[,1]] & z$safe[pairs[,2]]
  down_pair <- z$nonincreasing[pairs[,1]] & z$nonincreasing[pairs[,2]]
  local_pair <- safe_pair & down_pair &
    abs(z$q[pairs[,1]]-.30)<=.05 & abs(z$q[pairs[,2]]-.30)<=.05
  get_sequences <- function(ok) paste(z$sequence[ok],collapse="; ")
  distance <- function(ok) if(any(ok)) min(abs(z$q[ok]-.30)) else NA_real_
  data.frame(dose_monotone=all(apply(p,2,diff)>=-1e-12),
    cycle_nondecreasing=all(apply(p,1,diff)>=-1e-12),
    cycle_nonincreasing=all(apply(p,1,diff)<=1e-12),
    inversions=sum(gap>1e-12),pairs=nrow(pairs),
    safe_inversions=sum(gap[safe_pair]>1e-12),safe_pairs=sum(safe_pair),
    down_inversions=sum(gap[down_pair]>1e-12),
    local_safe_down_inversions=sum(gap[local_pair]>1e-12),
    max_inversion=max(0,gap),
    optimal_down=get_sequences(z$optimal_down),optimal_all=get_sequences(z$optimal_all),
    distance_down=distance(z$optimal_down),distance_all=distance(z$optimal_all),
    structural_gap=distance(z$optimal_down)-distance(z$optimal_all))
}

run_assumption_stress <- function(nsim=10000L,cores=32L,seed=2018L,
                                  reference=NULL) {
  started <- proc.time()[["elapsed"]]
  folder <- stress_directory()
  if(is.null(reference)) reference <- readRDS(file.path(folder,"data","table3.rds"))
  # Load the SAME simulation engine in this function's environment.
  # No altered allocation, safety, PAVA, tie or fallback algorithm is supplied.
  source(file.path(folder,"simulate.R"),local=environment())
  value <- function(key) {
    out <- reference$settings$value[reference$settings$setting==key]
    stopifnot(length(out)==1)
    out
  }
  specs <- stress_scenarios()
  probability <- stress_truth(specs)
  # Prove the common formula exactly reproduces the six benchmark surfaces.
  for(s in 1:6) {
    old <- subset(reference$truth,scenario==s)
    old <- old[order(old$cycle,old$dose),]
    stopifnot(max(abs(matrix(old$p,3,3)-probability[[s]]))<1e-12)
  }
  result <- simulate_table3(nsim=nsim,cores=cores,seed=seed,scenarios=7:10,
    truth=probability,weights=value("pava_weights"),nearest_tie=value("nearest_tie"),
    safety_min_n=as.integer(value("safety_min_n")),
    mts_screen_min_n=as.integer(value("mts_screen_min_n")),
    no_mts_fallback=value("no_mts_fallback"),global_monitor=value("global_monitor"))
  result$engine_runtime_seconds <- result$runtime_seconds
  result$scenario_parameters <- specs
  result$truth_audit <- cbind(specs,do.call(rbind,lapply(probability,stress_audit)))
  result$sequence_truth <- do.call(rbind,lapply(7:10,function(s) {
    cbind(scenario=s,stress_sequences(probability[[s]]))
  }))
  # Additional scores do not feed back into allocation or recommendation.
  z <- result$sequence_truth
  idx <- match(paste(result$trials$scenario,result$trials$selected_mts),
               paste(z$scenario,z$sequence))
  for(name in c("safe","optimal_down","optimal_all")) {
    result$trials[[paste0("selected_",name)]] <- !is.na(idx) & z[[name]][idx]
  }
  result$trials$selected_q <- z$q[idx]
  result$trials$unsafe_recommendation <- result$trials$selected_mts!="None" &
    !result$trials$selected_safe
  result$trials$optimal_down_seen <- result$trials$optimal_all_seen <- FALSE
  for(s in 7:10) for(type in c("down","all")) {
    optimal <- z$sequence[z$scenario==s & z[[paste0("optimal_",type)]]]
    seen <- unique(result$paths$trial[result$paths$scenario==s &
                                      result$paths$sequence %in% optimal])
    rows <- result$trials$scenario==s
    result$trials[rows,paste0("optimal_",type,"_seen")] <-
      result$trials$trial[rows] %in% seen
  }
  result$cells$true_p <- mapply(function(s,d,j) probability[[s]][d,j],
                              result$cells$scenario,result$cells$dose,result$cells$cycle)
  result$cells$unsafe <- result$cells$true_p>=.35
  result$job <- Sys.getenv(c("SLURM_JOB_ID","SLURM_JOB_PARTITION","SLURM_JOB_NODELIST",
                            "SLURM_CPUS_PER_TASK"))
  names_to_check <- c("maximum_patients","cohort_size","cycles","target","EI","xi",
    "safety_min_n","mts_screen_min_n","pava_weights","nearest_tie","clinical_ties",
    "no_mts_fallback","global_monitor","global_exclusion_start","global_evaluation",
    "outcome_generator","probability_cap")
  result$rule_check <- data.frame(setting=names_to_check,
    benchmark=vapply(names_to_check,value,character(1)),
    stress=result$settings$value[match(names_to_check,result$settings$setting)])
  stopifnot(all(result$rule_check$benchmark==result$rule_check$stress))
  files <- c("simulate.R","R/basyc.R","R/mtpi2.R","R/mts.R","R/assumption_stress.R")
  result$source_md5 <- setNames(unname(tools::md5sum(file.path(folder,files))),files)
  result$runtime_seconds <- unname(proc.time()[["elapsed"]]-started)
  result
}

summarize_assumption_stress <- function(result) {
  do.call(rbind,lapply(7:10,function(s) {
    t <- subset(result$trials,scenario==s)
    cells <- subset(result$cells,scenario==s)
    seen <- t$optimal_down_seen
    data.frame(scenario=s,n=nrow(t),pcs=mean(t$correct_selection),
      optimal_down=mean(t$selected_optimal_down),optimal_all=mean(t$selected_optimal_all),
      unsafe=mean(t$unsafe_recommendation),none=mean(t$selected_mts=="None"),
      early_stop=mean(t$stopped_early),fallback=mean(t$fallback_used),
      coverage=mean(seen),conditional_pcs=if(any(seen)) mean(t$selected_optimal_down[seen]) else NA_real_,
      mean_enrolled=mean(t$n_enrolled),mean_cycles=mean(t$patient_cycles),
      mean_dlts=mean(t$dlts),unsafe_exposure=sum(cells$n[cells$unsafe])/sum(cells$n))
  }))
}

if(sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  nsim <- if(length(args)>=1) as.integer(args[1]) else 10000L
  cores <- if(length(args)>=2) as.integer(args[2]) else 32L
  seed <- if(length(args)>=3) as.integer(args[3]) else 2018L
  result <- run_assumption_stress(nsim,cores,seed)
  output <- file.path(stress_directory(),"data","assumption_stress.rds")
  temporary <- paste0(output,".tmp")
  saveRDS(result,temporary)
  if(!file.rename(temporary,output)) stop("Could not save ",output)
  print(summarize_assumption_stress(result),row.names=FALSE)
  cat("Saved",nrow(result$trials),"trials. Total runtime:",result$runtime_seconds,"seconds.\n")
}
