#!/bin/bash
# Run locally: bash unity.sh [trials per scenario] [cores] [seed] [table3|stress]
# Example: bash unity.sh 10000 auto 2018  (all six scenarios)
# Upload -> Slurm on Unity -> one RDS file -> render the existing Rmd locally.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
PROJECT_DIR=$PWD
NSIM=${1:-10000}
CORES=${2:-auto}
SEED=${3:-2018}
STUDY=${4:-table3}
case "$STUDY" in
  table3) RESULT=table3; REPORT=reproduce_table2; SCENARIOS='original six' ;;
  stress) RESULT=assumption_stress; REPORT=simulation_exploration; SCENARIOS='four stress' ;;
  *) echo 'Study must be table3 or stress.' >&2; exit 1 ;;
esac
JOB_TIME=${UNITY_TIME:-01:00:00}
PARTITION=${UNITY_PARTITION:-stat}
MAX_CORES=${UNITY_MAX_CORES:-32}
for value in "$NSIM" "$MAX_CORES" "$SEED"; do
  [[ "$value" =~ ^[1-9][0-9]*$ ]] || { echo 'Use positive integers for trials, cores and seed.' >&2; exit 1; }
done
[[ "$CORES" == auto || "$CORES" =~ ^[1-9][0-9]*$ ]] || exit 1
[[ "$JOB_TIME" =~ ^[0-9]+:[0-5][0-9]:[0-5][0-9]$ && "$PARTITION" =~ ^[A-Za-z0-9_-]+$ ]] || exit 1

# Reuse one live connection so uploads/downloads do not repeat authentication.
# Password and Duo are entered privately; no credentials are saved.
OWN_SOCKET=0
if [[ -n "${UNITY_SOCKET:-}" ]]; then
  SOCKET=$UNITY_SOCKET
  ssh -S "$SOCKET" -O check asc-unity
else
  OWN_SOCKET=1
  SOCKET_DIR=$(mktemp -d /tmp/lyu-unity-XXXXXX)
  SOCKET="$SOCKET_DIR/ssh"
  echo 'Connect to Unity: enter your password and Duo responses when prompted.'
  ssh -M -S "$SOCKET" -o ServerAliveInterval=30 -fnN asc-unity
fi
cleanup() {
  if (( OWN_SOCKET )); then
    ssh -S "$SOCKET" -O exit asc-unity >/dev/null 2>&1 || true
    rm -f "$SOCKET"
    rmdir "$SOCKET_DIR" 2>/dev/null || true
  fi
}
trap cleanup EXIT

# Request up to 32 currently idle CPUs on ONE healthy node in the shared pool.
# Slurm still makes the actual allocation; do not use free CPUs without a job.
if [[ "$CORES" == auto ]]; then
  CORES=$(ssh -S "$SOCKET" -o BatchMode=yes asc-unity \
    "sinfo -N -h -p $PARTITION -o '%t %C'" | \
    awk -v cap="$MAX_CORES" '/^(idle|mix)[[:space:]]/ {split($2,a,"/"); if(a[2]>best) best=a[2]} END {n=(best<cap?best:cap); print(n>0?n:1)}')
fi
printf 'Requesting %s CPU cores in partition %s; %s scenarios, %s trials each.\n' "$CORES" "$PARTITION" "$SCENARIOS" "$NSIM"
ssh -S "$SOCKET" asc-unity 'mkdir -p trial-coding-exercise'

# Update changed source files. Keep generated outputs and Mac packages local.
# There is no --delete; unrelated remote files are preserved.
rsync -az -e "ssh -S $SOCKET" \
  --exclude='.git/' --exclude='.Rproj.user/' --exclude='.DS_Store' \
  --exclude='.Renviron*' --exclude='.Rhistory' --exclude='.RData' --exclude='.Rdata' \
  --exclude='renv/library/' --exclude='renv/staging/' --exclude='renv/cache/' \
  --exclude='renv/local/' --exclude='renv/python/' --exclude='.bootstrap-library/' \
  --exclude='*.html' --exclude='*/results/' --exclude='*_cache/' --exclude='*_files/' \
  --exclude='Lyu_2018/data/*.rds*' \
  --exclude='Drug_combo_Boin/data/*.RData' --exclude='Drug_combo_Boin/data/*.Rdata' \
  "$PROJECT_DIR/" asc-unity:trial-coding-exercise/

# The stress run inherits its clinical rules from this exact saved benchmark.
if [[ "$STUDY" == stress ]]; then
  rsync -az -e "ssh -S $SOCKET" Lyu_2018/data/table3.rds \
    asc-unity:trial-coding-exercise/Lyu_2018/data/table3.rds
fi

# --wait returns only after the Slurm job finishes, and propagates job failure.
# R runs on an allocated compute node. Do not edit/upload again during this job.
ssh -S "$SOCKET" asc-unity "bash -s -- $NSIM $CORES $SEED $JOB_TIME $PARTITION $STUDY" <<'REMOTE'
set -euo pipefail
cd "$HOME/trial-coding-exercise"
if [[ "$6" == stress ]]; then
  SCRIPT=Lyu_2018/R/assumption_stress.R
  LOG=Lyu_2018/assumption_stress.log
else
  SCRIPT=Lyu_2018/simulate.R
  LOG=Lyu_2018/table3.log
fi
sbatch --wait --parsable --job-name="lyu-$6" --partition="$5" \
  --nodes=1 --ntasks=1 --cpus-per-task="$2" --mem=8G --time="$4" \
  --output="$LOG" \
  --wrap="bash -lc 'set -e; module load r/4.5.1-5zezbzn; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1; Rscript --vanilla $SCRIPT $1 \${SLURM_CPUS_PER_TASK} $3'"
REMOTE

# Download only the successfully saved data; render on the Mac afterward.
rsync -az -e "ssh -S $SOCKET" \
  "asc-unity:trial-coding-exercise/Lyu_2018/data/$RESULT.rds" \
  "Lyu_2018/data/$RESULT.rds.download"
mv "Lyu_2018/data/$RESULT.rds.download" "Lyu_2018/data/$RESULT.rds"
Rscript --vanilla -e "rmarkdown::render('Lyu_2018/$REPORT.Rmd', quiet = TRUE)"
echo "Finished: Lyu_2018/$REPORT.html"
