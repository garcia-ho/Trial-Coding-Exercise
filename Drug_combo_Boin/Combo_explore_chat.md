# Combination-dose exploration: BOIN two-neighbor rules, wider neighborhoods, and parallel paths

**Context handoff for the Bayesian Trial Design / R workspace**  
**Subject:** The scientific and algorithmic distinctions between standard two-dimensional BOIN allocation, consideration of additional adjacent drug combinations (as in some model-based combination-dose frameworks), and parallel dose-finding paths.  
**Scope:** Background and open design questions from our discussion. This note deliberately does **not** prescribe an experiment, simulation workflow, software architecture, or R implementation.

## 1. Research motivation

In a two-agent phase I dose-finding trial, the candidate treatment space is a matrix of combinations of levels of drug A and drug B. We want to understand the consequences of **restricting each BOIN escalation/de-escalation decision to at most two axis-adjacent combinations** versus **allowing more nearby, potentially untested combinations**. The issue is particularly relevant in early cohorts, when most dose combinations have no DLT observations. The discussion also includes the **parallel-path extension** proposed by Lin and Yin for broader exploration and toxicity-contour identification.

These are three distinct ideas that should not be conflated:

1. **Size/geometry of the admissible neighborhood:** Which dose combinations may be considered as the *next* treatment location?
2. **Information borrowing:** Are toxicity estimates for a candidate informed by outcomes at other dose combinations?
3. **Number of simultaneous paths:** Does the design pursue a single cohort-allocation path or multiple concurrent paths?

The first, second, and third can vary independently in principle; standard BOIN and the Lin–Yin parallel-path extension instantiate different combinations of these design choices.

## 2. Basic notation and ordering

Let drug A have levels \(A_1,\ldots,A_J\) and drug B have levels \(B_1,\ldots,B_K\). At combination \((j,k)\):

- \(\pi_{jk}=P(\mathrm{DLT}\mid A_j,B_k)\): unknown true DLT probability.
- \(n_{jk}\): number of patients treated.
- \(y_{jk}\): number of observed DLTs.
- \(D_{jk}=(n_{jk},y_{jk})\): local observed data.
- \(\widehat\pi_{jk}=y_{jk}/n_{jk}\), defined for \(n_{jk}>0\): observed DLT proportion.
- \(\phi\): target DLT probability.

The common sampling model is

\[
Y_{jk}\mid n_{jk},\pi_{jk}\sim\mathrm{Binomial}(n_{jk},\pi_{jk}).
\]

A usual partial-order assumption is that toxicity cannot decrease if **one drug increases while the other is fixed**:

\[
\pi_{j+1,k}\geq\pi_{jk},\qquad \pi_{j,k+1}\geq\pi_{jk}.
\]

The two *exchange* combinations \((j+1,k-1)\) and \((j-1,k+1)\) are **not** ordered relative to \((j,k)\) by these assumptions, because one drug increases as the other decreases. In contrast, \((j+1,k+1)\) *is* ordered above \((j,k)\) if both monotonicities hold, but represents a simultaneous increase in both drugs. Consequently, **partial ordering alone does not uniquely force the exact two-neighbor BOIN rule**: excluding exchange moves and excluding simultaneous two-drug increases are distinct choices.

**Source:** Chapter 4, pp. 69–71, especially Figure 4.1, *Model-Assisted Bayesian Designs for Dose Finding and Optimization* (`references/Model_Assisted_Bayesian_Designs_Chapter4.pdf`).

## 3. What original Lin and Yin BOIN actually specifies

The original article (published online in 2015, journal issue 2017) calls the proposed method a **two-dimensional optimal interval design** (often abbreviated **INT** in its simulations); it is an extension of Bayesian optimal interval ideas to drug combinations.

At the current combination \((j,k)\), it **defines**, rather than derives in detail, the two admissible sets:

\[
\boxed{\mathcal A_E^{\mathrm{BOIN}}=\{(j+1,k),(j,k+1)\}}
\]

\[
\boxed{\mathcal A_D^{\mathrm{BOIN}}=\{(j-1,k),(j,k-1)\}}.
\]

These sets have **at most two candidates**, with only valid dose levels retained at boundaries. The current combination's cumulative DLT proportion is compared with the optimal interval \(I^*=(\lambda_e,\lambda_d)\). The paper's notation expresses this interval as \((\phi-\Delta_L,\phi+\Delta_U)\).

- Below the interval: select an admissible escalation candidate with the **largest posterior probability of having toxicity inside \(I^*\)**.
- Above the interval: apply the corresponding de-escalation rule.
- Inside the interval: stay at the same combination.
- If several admissible candidates attain the maximum probability, **randomly choose among those tied candidates with equal probability**.
- If there is no admissible candidate, retain the current dose, subject to other applicable safety/stopping rules.

The original paper uses **Jeffreys' independent prior** for each combination:

\[
\pi_{jk}\sim\mathrm{Beta}(0.5,0.5),\qquad
\pi_{jk}\mid D_{jk}\sim\mathrm{Beta}(y_{jk}+0.5,n_{jk}-y_{jk}+0.5).
\]

For target \(\phi=0.30\), its simulation uses the optimal interval approximately \((0.236,0.359)\). The textbook's later presentation uses a \(\mathrm{Beta}(1,1)\) local prior instead; **do not silently interchange those prior choices when interpreting the original article**.

**What the article does and does not say:** It explicitly defines the two-candidate neighborhoods and tie-breaking rules. It **does not give a detailed standalone argument demonstrating that exactly two adjacent candidates are optimal relative to four, eight, or some other neighborhood**. The explanation that axis-only movement respects unambiguous, one-drug-at-a-time escalation/de-escalation is a reasonable *interpretation*, not a quoted theorem proving this precise neighborhood is necessary.

**Source:** Lin & Yin, §2, journal pp. 2157–2158 (`references/Lin_Yin_BOIN_Combination_Original.pdf`), especially the paragraphs defining \(\mathcal A_E,\mathcal A_D\), allocation, tie breaking, boundaries, and the Jeffreys prior.

## 4. Untested neighbors and why random selection arises early

For an untested combination, \(n_{jk}=0\), so its local posterior remains its common prior. If both admissible candidates are untested, they have identical local posteriors and therefore identical probabilities of lying within \(I^*\):

\[
P(\pi_{j+1,k}\in I^*\mid D)=P(\pi_{j,k+1}\in I^*\mid D)
\]

when both have \((n,y)=(0,0)\) and the same prior. They tie and—if an escalation is required and both are admissible—are selected with probability \(1/2\) each.

**Concrete early-cohort illustration (cohort size 3):**

1. Start at \((1,1)\) as prescribed by the original single-path design.
2. Observe 0/3 DLTs, below the example lower boundary \(0.236\).
3. Escalation is indicated to \((2,1)\) or \((1,2)\). Both are untested; their target-interval posterior probabilities tie; choose randomly with equal probability.
4. If the newly chosen combination also has 0/3 DLTs, its two upward neighbors may likewise be untested and another random tie can arise.

**Inference:** Such randomized tie breaking **can be frequent at the beginning**, particularly along a path of repeated escalation into unvisited combinations. It is **not guaranteed after every early cohort**: the design might stay or de-escalate; one candidate may be out of bounds; one or both candidates may already have information; or candidates may have unequal posterior interval probabilities. Equal posterior probabilities can also arise after testing, not just for two untested combinations. The original paper gives the tie rule but **does not quantify the frequency of early random ties**.

A useful conceptual distinction: this is **random tie breaking caused by insufficient candidate-specific information**, **not** a designed preference for exploring the least-visited dose. With an independent local model, adding *more untested candidates* without further allocation rules can simply create **more equally scored candidates**.

**Source:** Lin & Yin, §2, pp. 2157–2158, particularly the equal-probability tie rule and local Jeffreys prior.

## 5. Why wider neighborhoods are a separate design issue

The uploaded textbook gives a **particular model-based copula-type framework**, not a universal rule for all CRM variants. In that example, when considering escalation from \((j,k)\), the candidate list includes up to **four** combinations:

\[
\mathcal A_E^{\mathrm{model\ example}}=
\{(j+1,k),\;(j,k+1),\;(j+1,k-1),\;(j-1,k+1)\}.
\]

The last two are **exchange moves**: one drug goes up while the other goes down. The framework uses estimated joint toxicity probabilities and a posterior rule to judge movement. A model can assign those exchange combinations estimated toxicity probabilities relative to the current location, because it borrows information across tested combinations through shared dose–toxicity parameters.

By contrast, standard BOIN's local beta-binomial updates provide **no joint model** telling us whether an exchange move increases or decreases toxicity. Broadening the neighborhood therefore is **not equivalent merely to adding two more candidates to the original BOIN allocation set**, unless an additional interpretation/rule is supplied for admissibility and direction of movement.

Other possible geographic neighbors, such as \((j+1,k+1)\), raise a **different** concern: moving both agents upward at once, despite their partial-order relationship, is not the same one-drug-at-a-time step used by standard BOIN. The original Lin–Yin paper acknowledges that the grid can have up to **eight adjacent combinations**; it does not compare all eight as candidate sets within the original single-path BOIN rule.

**Important wording:** It is accurate to say “some model-based combination frameworks consider four nearby candidate combinations, including exchange moves.” It is **not** accurate to say “all CRM or POCRM combination designs have exactly four neighboring candidates.” The actual admissible movement policy varies by model/design and implementation.

**Source:** Chapter 4, §4.2, pp. 71–73 for the copula-type candidate lists and joint model; Lin & Yin, Introduction, p. 2155 for up to eight adjacent combinations.

## 6. Information borrowing vs. neighborhood size

The methods are distinguished along two axes:

| Property | Original single-path BOIN / INT | Example joint model-based framework |
|---|---|---|
| Sequential DLT evidence | Local \((n_{jk},y_{jk})\) | All accumulated dose-combination DLT data |
| Candidate toxicity evaluation | Independent beta-binomial posterior per combination | Model-implied joint dose–toxicity estimates |
| When a candidate is untested | Posterior remains the common prior | May still have an informative model prediction |
| At most how many candidates in cited escalation set? | 2 axis neighbors | 4, including 2 exchange moves |
| Main methodological trade-off | Simple decisions that avoid joint-model misspecification | Cross-dose information can guide comparisons, but relies on model assumptions |

**Subtle point:** “No information borrowing in BOIN” refers to **sequential allocation**. At trial completion, Lin and Yin apply **bivariate isotonic regression** to order-restrict the toxicity estimates and select an MTD **among tested combinations**. Final order-restricted pooling is different from using a parametric model to borrow information during adaptive allocation.

The original paper's final MTD selection explicitly restricts eligible combinations to \(\mathcal N=\{(j,k):n_{jk}\neq0\}\). It does **not** guarantee reliable estimation at every untested grid point. Thus **finding one MTD** and **characterizing the full MTD contour / two-dimensional surface** are distinct scientific objectives.

**Sources:** Lin & Yin, §2, pp. 2157–2158; Chapter 4, §§4.2–4.3.1, pp. 71–75.

## 7. Parallel exploration is already discussed in the original BOIN paper

Lin and Yin propose **parallel dose-finding paths** in **§4.2, pp. 2163–2165**, labeled **INTP** in their simulation tables. This is a particularly relevant comparator because it is a broader-exploration extension *within the same interval-design family*, not just a comparison to a joint logistic/copula/POCRM model.

Their stated motivation is better coverage of the two-dimensional dose space and better characterization of **toxicity contours / multiple MTDs**.

**Key mechanism stated in the article:** If the current combination is below the target interval and the admissible escalation set contains untried combinations, the method can **assign a cohort to each of the two upward directions**, \((j+1,k)\) and \((j,k+1)\), **in parallel**. Paths recommending the same combination are **merged**. Subsequent paths may split or merge as data accumulate. This **does not** require adding the two *exchange* neighbors to the original \(\mathcal A_E\) set.

In the article's trial illustration, the parallel paths begin simultaneously at \((1,2)\) and \((2,1)\); the number of paths reaches **four** at the third enrollment sequence. The reported parallel allocation framework covers more combinations than the single path, with different patient concentration and safety consequences.

**Authors' reported trade-off:** Their comparison in §4.2 finds that the **single-path** design tends to concentrate more patients at the MTD(s), while the **parallel-path** design spreads patients more broadly across the matrix and can have higher probabilities of identifying multiple MTDs and better characterization of the dose–toxicity surface. Broader exploration can also mean more patients allocated to combinations above the MTD. These are outcomes under the simulated settings in the paper, not universal guarantees.

**Important distinction:**

- **Broader single-path candidate set:** one next cohort selects among more potential nearby sites (requires a policy for exchange/multidrug moves).
- **Parallel BOIN/INTP:** multiple cohorts can proceed along **both** original axis-upward directions, with path splitting/merging (greater coverage via multiple simultaneous paths).

These are **different exploration mechanisms** even though both can visit more previously untested dose combinations.

**Source:** Lin & Yin, §4.2, pp. 2163–2165, particularly the paragraph introducing the simultaneous two-direction assignment, Table 5 (parallel illustration), and the INT vs. INTP comparisons in Table 6.

## 8. Scientific distinctions to preserve in future discussion

The central question is not just whether “more untested neighbors” increases exploration. At least four conceptual outcomes must be kept distinct:

- **Early tie frequency:** Under independent priors, how much of early candidate selection is determined by equal posterior probabilities rather than observed toxicity evidence?
- **Spatial coverage vs. MTD concentration:** A wider admissible neighborhood or multiple paths may reach more grid points, but may direct fewer patients to a particular near-target combination.
- **Safety and partial ordering:** Exchange directions have unknown toxicity ordering; jointly escalating both drugs is ordered but makes a larger move. Broader coverage need not imply comparable overdose exposure.
- **Objective alignment:** A design appropriate for **finding one MTD** may not be appropriate for **estimating a multi-combination MTD contour**.

The original sources offer methods and simulated evidence, but **do not establish that an expanded four- or eight-neighbor model-assisted rule is superior to the two-neighbor BOIN rule**, nor do they provide the definitive rationale/optimality proof for restricting BOIN specifically to two neighbors. Such an expanded rule is conceptually **a new allocation-policy variant**, not a quoted feature of the original BOIN article.

## 9. Local reference files and source map

The PDFs are included in the accompanying package in the **`references/`** directory. These are copies of the user's uploaded originals, not rewritten or revised texts.

1. **Lin, R. & Yin, G.** “Bayesian optimal interval design for dose finding in drug-combination trials.” *Statistical Methods in Medical Research*, **26**(5), 2155–2167 (2017 journal issue; first published 2015). DOI: **10.1177/0962280215594494**. File: `references/Lin_Yin_BOIN_Combination_Original.pdf`. Focus on Introduction p. 2155; §2 pp. 2157–2158; §4.2 pp. 2163–2165.
2. **Model-Assisted Bayesian Designs for Dose Finding and Optimization**, Chapter 4, “Drug-Combination Trials,” pp. **69–83**. File: `references/Model_Assisted_Bayesian_Designs_Chapter4.pdf`. Focus on Figure 4.1 and partial ordering, §4.2 model-based methods (esp. model-based four-candidate moves), §4.3.1 BOIN combination, and §4.3.3 waterfall design for identifying the MTD contour.

*Textbook/context differences worth preserving:* The chapter's BOIN combination subsection uses a \(\mathrm{Beta}(1,1)\) prior and a particular two-neighbor rule. The original Lin–Yin journal article uses \(\mathrm{Beta}(0.5,0.5)\), labels the single-path design **INT**, and presents parallel-path **INTP**. The four-neighbor set in the textbook is an example from its **model-based copula-type method**, not a specification of the POCRM algorithm generally.
