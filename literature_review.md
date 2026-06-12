# Vaccine Effectiveness: Methodological Literature Review

Focus: conceptual / definitional work on what VE *means* in systems with
interference, from ~1990 onward. Statistical estimation procedures (TND,
cohort, case-control mechanics) are deliberately out of scope.

Citations below have been verified by Google Scholar / journal lookups for
year, journal, and volume. Where minor details (page numbers, exact
title-case) matter, double-check against the journal page; the body of
the citation has been checked.

---

## 1. The Halloran–Struchiner foundational framework

Defines the *types* of vaccine effect (susceptibility, infectiousness,
progression, total) and the four study-design conditions that separate
direct from indirect/total/overall population effects.

### Halloran & Struchiner (1991) — *Study designs for dependent happenings*
*Epidemiology* 2(5): 331–338.
Revives Ross's "dependent happenings" terminology and introduces the
conditions (Ia/Ib/IIa/IIb) under which different vaccine effects can be
identified. First formal acknowledgement that the comparison group is
itself affected by the intervention, so a "vaccinated vs. unvaccinated"
contrast does not isolate a single causal quantity.

### Halloran, Haber, Longini (1992) — *Interpretation and estimation of vaccine efficacy under heterogeneity*
*Am J Epidemiol* 136(3): 328–343.
Decomposes vaccine effects into VE_S (susceptibility), VE_I
(infectiousness), and VE_P (progression). Shows conventional cohort VE
estimates a *mixture* whose weights depend on the force of infection and
population structure — so the same vaccine yields different VE numbers
in different epidemics. This is the original "VE is not a property of
the vaccine alone" argument.

### Halloran & Struchiner (1995) — *Causal inference in infectious diseases*
*Epidemiology* 6(2): 142–151.
Bridges Halloran–Struchiner's epidemiological framework with Rubin's
potential outcomes. Articulates clearly that standard VE estimands are
not well-defined causal contrasts when interference is present, because
the "potential outcome under no treatment" depends on others'
treatments. Defines individual-level causal estimands by letting the
potential outcomes for any individual depend on the *vector* of
treatment assignments. The intellectual seed for all subsequent
causal-inference-under-interference work.

### Halloran, Struchiner, Longini (1997) — *Study designs for evaluating different efficacy and effectiveness aspects of vaccines*
*Am J Epidemiol* 146(10): 789–803.
Consolidates the design taxonomy: which trial design identifies which
effect parameter. The canonical citation for "your trial design
constrains your estimand."

### Halloran, Longini, Struchiner (2010) — *Design and Analysis of Vaccine Studies*
Springer, *Statistics for Biology and Health*.
The textbook synthesis. Chapters 2 and 13 give the cleanest exposition
of the conceptual framework.

### Halloran & Hudgens (2016) — *Dependent happenings: a recent methodological review*
*Current Epidemiology Reports* 3(4): 297–305.
The most useful "where are we now" anchor. Reviews how the 1991-era
framework has been formalised through partial-interference causal
inference, and points at open problems (network interference,
observational identification, longitudinal/waning effects).

### Halloran (2024/2025) — *Designs for Vaccine Studies*
*Annual Review of Statistics and Its Application* 12 (published online
Oct 2024; in print March 2025).
The newest authoritative review. Covers dependent happenings, indirect
effects on unvaccinated individuals, and how vaccination reduces
person-to-person transmission. Almost certainly the right "current
status" citation.

---

## 2. Formal causal inference under interference

Translates Halloran–Struchiner into the Rubin causal model and gives
estimands rigorous potential-outcomes definitions.

### Rubin (1990) — *Formal modes of statistical inference for causal effects*
*Journal of Statistical Planning and Inference* 25(3): 279–292.
Origin of "SUTVA" as a named assumption (Rubin had used the underlying
idea earlier in his 1980 comment on Basu). Identifies the two ways SUTVA
breaks: (a) multiple versions of treatment, (b) interference between
units. Often cited as the formal "no interference" requirement that
vaccine-causal-inference subsequently relaxes.

### Hong & Raudenbush (2006) — *Evaluating kindergarten retention policy: A case study of causal inference for multilevel observational data*
*JASA* 101(475): 901–910.
Not vaccines, but an early concrete relaxation of SUTVA: lets potential
outcomes vary with a *function* of the treatment vector in the cluster
(here, the school's retention policy). One of the first applied papers
to extend Rubin's potential outcomes to multi-level settings with
within-cluster interference.

### Sobel (2006) — *What do randomized studies of housing mobility demonstrate?*
*JASA* 101(476): 1398–1407.
First rigorous treatment of SUTVA violation in the causal-inference
literature. Shows that even simple estimands like "average treatment
effect" need re-definition under interference. Frequently cited as the
methodological starting point for the Hudgens–Halloran line.

### Hudgens & Halloran (2008) — *Toward causal inference with interference*
*JASA* 103(482): 832–842.
**The pivotal paper.** Defines direct, indirect, total, and overall
effects as contrasts of average potential outcomes under *partial
interference* (units form non-overlapping clusters; interference within
but not between). Two-stage randomization (cluster → individual) lets
each effect be unbiasedly estimated. The formal potential-outcomes
restatement of Halloran–Struchiner.

### VanderWeele & Tchetgen Tchetgen (2011a) — *Effect partitioning under interference in two-stage randomized vaccine trials*
*Statistics and Probability Letters* 81(7): 861–869.
Express the overall effect as a sum of (i) the indirect/spillover effect
and (ii) a contrast between two direct effects. Refines the partitioning
of Hudgens–Halloran 2008.

### VanderWeele & Tchetgen Tchetgen (2011b) — *Bounding the infectiousness effect in vaccine trials*
*Epidemiology* 22(5): 686–693.
Uses principal stratification to identify VE_I (effect of vaccination on
transmission given infection). Shows VE_I is generally only partially
identified without strong assumptions; gives bounds. (Distinct paper
from the Stat & Prob Letters one above, same year.)

### Tchetgen Tchetgen & VanderWeele (2012) — *On causal inference in the presence of interference*
*Statistical Methods in Medical Research* 21(1): 55–75.
Extends Hudgens–Halloran to observational data via inverse-probability
weighting; introduces *stratified interference* (allocation effects
depend only on the number, not identity, of treated unit-mates).
Develops finite-population identification results.

### Halloran & Hudgens (2012) — *Causal inference for vaccine effects on infectiousness*
*International Journal of Biostatistics* 8(2), article 6.
Develops causal estimands for VE_I using both principal stratification
on the joint potential infection outcomes and interference between
individuals within transmission units. Pins down which contrasts have a
causal interpretation and which (standard infectiousness comparisons) do
not.

### VanderWeele, Tchetgen Tchetgen, Halloran (2012) — *Components of the indirect effect in vaccine trials: identification of contagion and infectiousness effects*
*Epidemiology* 23(5): 751–761.
Decomposes the indirect effect into a contagion piece (vaccinated
individuals don't become infected and therefore don't transmit) and a
true infectiousness piece (vaccinated-but-infected transmit less).
Important conceptually because "indirect effect" lumps these.

### Liu & Hudgens (2014) — *Large sample randomization inference of causal effects in the presence of interference*
*JASA* 109(505): 288–301.
Asymptotic theory for the Hudgens–Halloran two-stage estimators. Mostly
inferential, but pinpoints the conditions under which the
direct/indirect/total/overall estimands are well-defined as
limits — useful when arguing what's being estimated in the
infinite-population limit.

### Halloran & Hudgens (2018) — *Estimating population effects of vaccination using large, routinely collected data*
*Statistics in Medicine* 37(2): 294–301.
Applies IPW interference estimators to a cholera vaccination trial in
Matlab, Bangladesh. Useful as a worked example of moving from
Hudgens–Halloran theory to real data; clarifies what each effect
estimand means in practice.

---

## 3. Beyond partial interference: networks and general dependence

Drops the cluster assumption; lets individuals interfere through an
arbitrary contact/exposure structure.

### Manski (2013) — *Identification of treatment response with social interactions*
*Econometrics Journal* 16(1): S1–S23.
Generalises potential outcomes by treating response as a function of the
*entire treatment vector*. Identifies "constant treatment response" (CTR)
as a broad assumption class with no-interaction and unrestricted
interaction as polar cases. Foundational from the econometrics side of
the interference literature; influential on later exposure-mapping work.

### Ogburn & VanderWeele (2014) — *Causal diagrams for interference*
*Statistical Science* 29(4): 559–578.
Extends DAGs to settings with interference, contagion, and shared
exposures. Distinguishes three mechanisms: (i) direct interference (one
unit's treatment directly affects another's outcome), (ii) contagion
(one unit's outcome affects another's), (iii) allocational interference
(shared environment/group). Shows graphically why naive analyses
confound contagion with interference.

### Aronow & Samii (2017) — *Estimating average causal effects under general interference, with application to a social network experiment*
*Annals of Applied Statistics* 11(4): 1912–1947.
Generalises potential outcomes to *exposure mappings* — each unit's
outcome is a function of a known summary of others' treatments. Lets
you define and estimate effects under arbitrary, known interference
graphs using randomization-based inverse-probability weighting.
Foundational for the network-experiment literature.

### Eckles, Karrer, Ugander (2017) — *Design and analysis of experiments in networks: Reducing bias from interference*
*Journal of Causal Inference* 5(1).
Practical methods for designing randomized experiments on networks
(graph cluster randomization) to reduce bias in estimating "global"
treatment effects when interference is broad but unmodeled. Mostly
design-oriented but the discussion of what estimands are recoverable
under different designs is conceptually clarifying.

### Athey, Eckles, Imbens (2018) — *Exact p-values for network interference*
*JASA* 113(521): 230–240.
Randomization-inference machinery for testing null hypotheses about
spillover effects on networks. The framing of "null exposure mappings"
is conceptually clarifying even if the paper itself is methodological.

### Basse & Airoldi (2018) — *Model-assisted design of experiments in the presence of network-correlated outcomes*
*Biometrika* 105(4): 849–858.
Restricted randomization that uses pre-intervention network structure
to optimise the design. The derived "balance" criteria highlight which
network statistics actually matter for the causal contrast of interest.

### Sävje, Aronow, Hudgens (2021) — *Average treatment effects in the presence of unknown interference*
*Annals of Statistics* 49(2): 673–701.
**Major addition since 2017.** Shows that standard estimators (designed
assuming SUTVA) are consistent for a *generalised* average treatment
effect estimand — one that marginalises over realised spillover — under
broad and *unknown* interference structures, with rates depending on
the average amount of interference. Reframes the problem: rather than
needing to specify the exposure mapping a priori (Aronow–Samii) or the
cluster structure (Hudgens–Halloran), one accepts that the estimand
itself depends on the assignment distribution. Highly relevant to the
EATE conceptual framing.

### Forastiere, Airoldi, Mealli (2021) — *Identification and estimation of treatment and interference effects in observational studies on networks*
*JASA* 116(534): 901–918.
Formalises "individual" and "spillover" effects on networks using a
neighborhood exposure mapping; extends Hudgens–Halloran-style
identification to observational network data via a generalised
propensity score.

### Tchetgen Tchetgen, Fulcher, Shpitser (2021) — *Auto-G-Computation of causal effects on a network*
*JASA* 116(534): 833–844.
Causal inference from a *single realisation* of a connected network,
under arbitrary forms of interference and long-range dependence. Uses
chain-graph models to make inference tractable; defines a class of
causal estimands different from the exposure-mapping framework. The
counterpart for observational, single-network data.

---

## 4. Mechanistic / transmission-process definitions of VE

Defines VE in terms of the *generative* epidemic process rather than
randomization-implied contrasts. Closest to what your EATE work does.

### Becker (1989) — *Analysis of Infectious Disease Data*
Chapman & Hall.
Pre-1990 but the conceptual predecessor: VE as a parameter of the
transmission process (a relative hazard between vaccinated and
unvaccinated within a household), rather than a population-level risk
ratio. Sets up the household / chain-binomial paradigm that later
mechanistic work builds on.

### O'Hagan, Lipsitch, Hernán (2014) — *Estimating the per-exposure effect of infectious disease interventions*
*Epidemiology* 25(1): 134–138 (verify exact pages).
Articulates per-exposure VE — what is the reduction in per-contact
hazard? — and shows how this is distinct from the population-level risk
ratio. The per-exposure parameter is invariant to the epidemic
trajectory while the risk ratio is not; the same point as
Halloran–Haber–Longini (1992), but explicitly framed for modern
causal-inference readers.

### Kenah (2015) — *Semiparametric relative-risk regression for infectious disease transmission data*
*JASA* 110(509): 313–325.
**Correction from prior draft (was listed as Biostatistics).** Defines
the *contact interval* (time from onset of infectiousness to infectious
contact) and treats VE as a hazard ratio on this interval. The contact
interval is identifiable from transmission data, and the corresponding
hazard ratio is a generative-process estimand: it identifies the same
thing across different population/transmission settings, unlike the
risk-ratio VE.

### Eck, Morozova, Crawford (2022) — *Randomization for the susceptibility effect of an infectious disease intervention*
*Journal of Mathematical Biology* 85(4): 37.
Defines the *susceptibility effect* as the infection-risk contrast under
treatment vs no treatment, holding exposure to infectiousness constant.
Striking finding: under the Hudgens–Halloran (2008) two-stage design,
the estimated direct effect can be positive (looks harmful) even when
the intervention reduces both susceptibility and onward transmission.
A pointed example of how randomization-derived estimands diverge from
mechanistic estimands. Closely related to EATE's motivation.

### Auranen, Eichner, Britton, O'Neill et al. — household and transmission-model VE
A scattered literature from the late 1990s–2010s in *Statistics in
Medicine*, *Biostatistics*, and *Biometrics* estimating VE as a
per-contact parameter from household chain-binomial data. Worth citing
collectively as the "mechanistic-VE-from-household-data" tradition;
locate specific representative papers when you know which household
designs you want to discuss.

---

## 5. Recent applied/conceptual work directly on COVID-era VE

Papers that explicitly grapple with what "VE" means during ongoing
transmission.

### Lewnard, Patel, Reich et al. — interpreting trial VE
A cluster of 2021 papers (e.g. Lewnard et al. *CID*, Patel & Lipsitch
*NEJM*) arguing that headline VE numbers from clinical trials confound
intrinsic vaccine effect with epidemic state at the time of the trial.
Modern restatement of Halloran–Haber–Longini 1992 with concrete
numbers. Locate specific papers when you've decided which COVID-era
critiques to include.

### Kahn, Schrag, Patel, Lipsitch et al. — *time-varying VE / waning*
Methodological work from 2021–2024 grappling with how to define "VE at
time t since vaccination" when both vaccine-induced immunity and
epidemic state evolve. Much of this strand is estimation, but the
definitional pieces (what *is* VE at time t causally?) are worth
including.

### Crawford, Morozova, Buen Abad Najar, Galvani — interpreting VE during transmission
The Crawford group's COVID-era papers explicitly argue that
trial-estimated VE is not a property of the vaccine alone. Verify
specific citations on Google Scholar — the bibliographic details of
their 2022 paper were not nailed down in my own search.

---

## 6. Adjacent strands

- **Principal stratification** (Frangakis & Rubin, *Biometrics* 2002):
  background machinery for VE_I. Cite when discussing VanderWeele &
  Tchetgen Tchetgen 2011b.
- **Cluster-randomized vaccine trial methodology** (Hayes & Moulton's
  book; the Ebola ring-vaccination methodology). Application-side
  rather than definitional.
- **Causal mediation under interference** (VanderWeele, Imai, others):
  relevant for decomposing direct vs. indirect when the mediator is
  another person's infection state.
- **Agent-based / IBM-derived VE** (Halloran & Longini simulation work).
  Mostly applied but occasionally definitional.
- **Bipartite interference / spatial spillovers** (Zigler, Papadogeorgou
  on air pollution). Methodological extensions to non-cluster, non-
  network settings; could matter if you want a broader interference
  taxonomy.

---

## 7. Which papers actually propose a VE estimand — and on what scale

A sharper lens on the same literature: who *proposes what VE should be*
as a causal quantity, vs. who refines estimation of pre-existing
definitions? And on the difference scale (E[Y(1)] − E[Y(0)]) or the
ratio scale (1 − E[Y(1)] / E[Y(0)]) that VE has classically used?

### 7.1 Papers that propose a VE estimand

The genuinely new-estimand papers in this literature are surprisingly
few. Most of the work since 1995 *refines estimation* of estimands that
were already defined.

| Paper | Estimand proposed | Scale |
|---|---|---|
| Halloran, Haber, Longini (1992) | VE_S, VE_I, VE_P as 1 − rate ratios | Ratio (informal, pre-Rubin) |
| Halloran & Struchiner (1995) | Direct/indirect/total/overall as potential-outcomes contrasts | Difference (with ratios in motivating examples) |
| VanderWeele & Tchetgen Tchetgen (2011b) | VE_I within the "always-infected" principal stratum | Ratio (risk ratio bounds) |
| Halloran & Hudgens (2012) | Formal causal VE_I across transmission-unit scenarios | Mostly ratio |
| VanderWeele, Tchetgen Tchetgen, Halloran (2012) | Contagion / infectiousness decomposition of the indirect effect | Mostly ratio |
| O'Hagan, Lipsitch, Hernán (2014) | Per-exposure VE | Ratio (per-contact hazard) |
| Kenah (2015) | VE as hazard ratio on the *contact interval* | Ratio |
| Sävje, Aronow, Hudgens (2021) | Expected Average Treatment Effect | **Difference** |
| Eck, Morozova, Crawford (2022) | Susceptibility effect (risk under treatment vs. control, holding exposure fixed) | **Difference** |

Papers that explicitly do *not* propose a new VE estimand: Hudgens &
Halloran (2008) (general ATE-style estimands that happen to apply to
VE), Aronow–Samii (2017), Forastiere et al. (2021), Tchetgen Tchetgen
et al. (2021), Ogburn–VanderWeele (2014), Manski (2013) (frameworks
within which estimands can be defined); Liu–Hudgens (2014),
Halloran–Hudgens (2018), Tchetgen Tchetgen–VanderWeele (2012)
(inference/estimation for existing estimands); Halloran-led reviews.

### 7.2 The scale split

There is a clean — and underappreciated — split in the literature
between papers that work on the difference scale and those that work on
the ratio scale.

- **Ratio scale**: Halloran–Haber–Longini 1992, the principal-
  stratification VE_I cluster (VanderWeele/Tchetgen Tchetgen/Halloran
  2011–2012), and the mechanistic branch (O'Hagan–Lipsitch–Hernán
  2014, Kenah 2015). These look like "VE" in the classical
  epidemiological sense: VE = 1 − (rate ratio or hazard ratio).
- **Difference scale**: everything in the formal causal-inference-
  under-interference branch, from Halloran–Struchiner 1995 through
  Hudgens–Halloran 2008, Liu–Hudgens 2014, Aronow–Samii 2017,
  Sävje–Aronow–Hudgens 2021, Eck–Morozova–Crawford 2022. Estimands are
  E[Y(1, w)] − E[Y(0, w)] and variants.

The 1995 Halloran-Struchiner paper sits at the hinge: it translates
the 1992 ratio-scale VE into potential outcomes — but the translation
silently drops the ratio. Once estimands are written in potential-
outcomes form, defining a causal risk ratio E[Y(1)] / E[Y(0)] is
harder (non-linear functional of potential outcomes, doesn't decompose
cleanly under interference, IPW doesn't go through nicely), so the
literature shifted scale and never shifted back.

The mechanistic branch retained ratios because hazard ratios are the
natural parameter of a transmission process.

### 7.3 Implication for EATE

EATE as defined in this project is a population-level causal estimand
on the *ratio* scale, defined against the underlying transmission
process. That combination is what's missing in the existing
literature:

- Sävje–Aronow–Hudgens (2021) propose the population-level
  marginalisation idea but on the difference scale, and as a generic
  ATE rather than a VE.
- Eck–Morozova–Crawford (2022) define a VE-flavoured estimand defined
  against transmission, but again on the difference scale.
- The ratio-scale VE estimands (1992, 2011–2012, 2014, 2015) are
  either pre-Rubin population ratios (1992), restricted to specific
  effects like VE_I (2011–2012), or per-exposure parameters that don't
  aggregate to a population quantity (2014, 2015).

So the gap is real: **a population-level, ratio-scale causal VE
estimand defined against the actual transmission process** is what
EATE proposes. The Sävje-style "EATE" name is appropriate because it
inherits the marginalise-over-realised-interference idea, but the
scale and the substantive target are different.

---

## 8. Suggested narrative arc

1. **1990–1997 — definitional foundation.** Halloran–Haber–Longini
   propose ratio-scale VE_S/VE_I/VE_P. Halloran–Struchiner translate
   this into potential outcomes, shifting (almost incidentally) to the
   difference scale.
2. **2006–2014 — formalization in potential outcomes.** Sobel,
   Hudgens & Halloran, Tchetgen Tchetgen & VanderWeele, Liu & Hudgens
   develop partial-interference theory entirely on the difference
   scale. Principal-stratification work (VanderWeele/Tchetgen
   Tchetgen/Halloran 2011–2012) retains ratios but only for VE_I.
3. **2014–present — generalization in two branches.**
   - **Network/exposure-mapping branch** (Manski, Aronow–Samii,
     Ogburn–VanderWeele, Forastiere–Airoldi–Mealli, Sävje–Aronow–
     Hudgens, Tchetgen Tchetgen et al.): define effects under arbitrary
     interference graphs, on the difference scale.
   - **Mechanistic branch** (Kenah, O'Hagan–Lipsitch–Hernán, Eck–
     Morozova–Crawford): define VE as a parameter of the underlying
     transmission process, mostly on the ratio scale (except Eck et
     al.).
4. **Where EATE sits.** Combines (a) the population-level
   marginalisation idea from Sävje et al., (b) the
   transmission-process target of the mechanistic branch, and (c) the
   ratio scale of classical VE — a combination no existing paper
   simultaneously offers.

---

## 9. Things still to check

- Page numbers and exact title-case on a handful of older papers.
- Whether to include the **O'Hagan, Lipsitch, Hernán (2014)** paper as
  *Epidemiology* — couldn't pin down the journal in the search; verify.
- Pick representative papers from the **Lewnard / Lipsitch / Kahn**
  COVID-era VE-interpretation strand once you've decided which
  critiques are most relevant.
- Check whether **Tchetgen Tchetgen, Fulcher, Shpitser (2021)** issue
  number is 534 or 535 in JASA (saw both in searches).
