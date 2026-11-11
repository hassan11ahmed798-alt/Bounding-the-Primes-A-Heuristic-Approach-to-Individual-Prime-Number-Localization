# Preliminary Technical Research Study of the PN → PI Model: Mathematical Structure, A/B/C/F Calibration, and Scalability Limits

**Status:** Source-code-grounded audit draft, **not** a final peer-reviewed paper or a proof of universal containment.  
**Reference implementation:** `PrimeCloud_v18_2_PN_to_PI_Original_Mode_Limits.py`.  
**Non-modification policy:** No changes to model coefficients, calibration formulas, pruning, anomaly handling, or interval-extraction logic.  
**Prepared:** November 11, 2026.

## Abstract

This study examines a heuristic procedure for bounding the prime of index \(n\). A PN model generates candidate locations, a PI model tests those candidates against an approximation to the prime-counting function, and a final contiguous integer interval is constructed from selected admissible components. The study distinguishes published theoretical prime bounds from the model's non-certified empirical interval. It documents the A/B/C/F calibration routes and identifies implementation-level obstacles to numerical scalability. In particular, the current PI-component extractor scans integers individually and explicitly rejects search regions whose upper endpoint exceeds \(2^{53}\). Consequently, successful calibration at \(N=2^{50}\) does not by itself establish that the implementation can finish bounding \(p_{2^{51}}\). This draft does not establish a mathematically validated replacement for that extractor.

## 1. Scope and notation

- \(n\): target prime index, a positive integer.
- \(p_n\): the \(n\)-th prime, with \(p_1=2\).
- \(N\), or `TRAIN_END`: greatest prime index permitted in calibration.
- \(\pi(x)\): the number of primes not exceeding \(x\).
- \(R(x)\): the reference function implemented as `R_vec`; it must not be identified unconditionally with \(\pi(x)\).
- \(R^{-1}(n)\): a numerically approximated inverse of the reference function, used to form the initial PN candidate.
- \(F_{PN}(n)\): a prime-location approximation formed from the inverse reference, a Ridge correction, and a median adjustment.
- \(F_{PI}(x)\): an approximation to the prime-counting function incorporating \(R\), a PI bias, a Ridge correction, and a median adjustment.
- \(C_{\mathrm{EMP},PN}\), \(C_{\mathrm{EMP},PI}\): empirical margin coefficients; they are **not** established universal constants.
- \(L_n,U_n\): theoretical endpoints returned by `certified_nth_prime_bounds`; the function name alone does not verify that every bound's hypotheses and numerical implementation are valid.
- `RAW`: the cloud before application of the empirical `C_EMP` margin.
- `post-RAW`: the application of the empirical PI margin after constructing the raw cloud.
- `checkpoint`: a discrete branch-state update location.
- `branch`: a candidate state with displacement, slope, phase/wave, and state corrections; pruning removes certain states under the existing rules.

## 2. Computational pipeline: from target index to interval

1. **Input validation.** `_require_exact_int` rejects floating-point inputs whose exact integer identity cannot be assured beyond \(2^{53}\), while accepting Python integers.
2. **Theoretical baseline.** `certified_nth_prime_bounds(n)` returns \([L_n,U_n]\) before calibration.
3. **Mode selection.** `--calibration-mode auto` selects A, B, C, or F based on \(N\), not directly on \(n\).
4. **Calibration.** A/B call `compute_calibration` with stored prime data; C/F call `calibration_from_pattern` using `primecount` queries.
5. **PN center.** `jump_predict` computes `inv_R_scalar(float(n_target))` and adds `cal_pn['c_vec'](...)`.
6. **PN branch walk.** `jump_walk_pn` generates raw candidates; the base geometry uses individual candidate values rather than an independently expanded PN-only envelope.
7. **PI branch walk.** `jump_walk_pi` generates snapshots around a center, using `CP_PI` and logarithmic scaling.
8. **Raw acceptance.** For a PN candidate \(x_i\), test whether \(\pi_{\rm raw,lo}(x_i)\le n\le\pi_{\rm raw,hi}(x_i)\).
9. **Empirical acceptance.** Expand PI counting bounds, in counting units rather than \(x\)-units, by \(M_{PI}(x)=\mu C_{\rm EMP,PI}(\log x)^2\).
10. **Integer components.** `_pi_raw_and_empirical_components` checks every integer in the search range and returns all disconnected acceptance components before selection.
11. **Component selection.** `_pi_components_seeded_by_pn` retains PI components that intersect accepted PN seeds.
12. **Final interval.** If components exist, return the smallest contiguous integer interval spanning all selected components, including gaps; otherwise invoke the theoretical fallback and report the diagnostic status.

## 3. Equations explicitly represented in the source code

### 3.1 PN reference

Within `jump_predict`,

\[
x_0=R^{-1}(n),\qquad x_{\rm target}=x_0+c_{PN}(n).
\]

**Important qualification:** For C/F, `calibration_from_pattern` defines `c_vec` as a zero function. Thus the PN correction available under A/B must not be assumed to carry over unchanged to C/F. Separate equivalence testing is required.

### 3.2 PI reference

The code defines `Fb_pi(x)` using `R_vec(x)` and a bias term of the form

\[
B_{PI}(x)=A_{PI}\sqrt{\frac{x}{\log(\max(x,3))}}+B_{PI,0}.
\]

A Ridge correction based on features associated with 60 selected zeros of the Riemann zeta function and a residual median adjustment are then added. **Zeta-zero-inspired features do not by themselves constitute a rigorously derived probabilistic model or an error bound from the explicit formula.** Their empirical selection must be distinguished from analytic number-theoretic claims.

### 3.3 Acceptance tests

\[
\operatorname{RAW}(x)\iff \ell_{PI}(x)\le n\le u_{PI}(x).
\]

\[
\operatorname{EMP}(x)\iff
\ell_{PI}(x)-\mu C_{\mathrm{EMP},PI}(\log x)^2\le n
\le u_{PI}(x)+\mu C_{\mathrm{EMP},PI}(\log x)^2.
\]

The base implementation uses \(\mu=1\) unless a different value is supplied. **`C_EMP_PN` does not expand PN candidate points in this particular geometry path**; the code explicitly distinguishes this behavior.

### 3.4 Contiguous output interval

If selected components are \(J_1,\ldots,J_k\), the returned set is

\[
I_{\rm final}=\left[\min_i\inf J_i,\;\max_i\sup J_i\right]\cap\mathbb Z.
\]

This operation fills gaps. It does **not** establish that \(p_n\) lies in any empirical component in the first place.

## 4. Original automatic mode thresholds

Source: `AUTO_A_MAX`, `AUTO_B_MAX`, `AUTO_C_MAX`, and the CLI dispatch in the reference implementation.

| Training condition | Mode | Implementation route |
|---|---|---|
| \(N\le2\times10^9\) | A | Store prime data; direct calibration |
| \(2\times10^9<N\le2\times10^{11}\) | B | Compressed storage; stream `C_EMP` evaluation |
| \(2\times10^{11}<N\le2\times10^{13}\) | C | `primecount` and sampled `C_EMP` measurement |
| \(N>2\times10^{13}\) | F | `primecount` and extrapolation of `C_EMP` from previous calibration patterns |

**Caveat:** The current CLI's B route does **not** eliminate storage of prime data. These thresholds are original design settings, not resource limits experimentally certified for 12 GB RAM.

## 5. What the benchmarks actually measured

The user supplied 24 successful calibration records for \(N=10^4\) through \(10^6\), followed by 24 more for \(N=2\times10^6\) through \(10^8\). These demonstrate **calibration completion**, not final-interval coverage. At \(N=10^8\):

| Mode | Calibration time (s) | Peak RSS (bytes) | `C_EMP_PN` | `C_EMP_PI` |
|---|---:|---:|---:|---:|
| A | 592.2505 | 974241792 | 37.6705213 | 1.2844031 |
| B | 533.2362 | 983617536 | 37.6705213 | 1.2844031 |
| C | 94.1016 | 299421696 | 39.5797404 | 1.3652741 |
| F | 11.3549 | 298979328 | 65.7980142 | 2.6028685 |

Mode F cannot be preferred merely because it is faster: the margin coefficients differ, and final-interval coverage has not been evaluated on an independent validation set. Nor may the timing trend observed at \(10^8\) be extrapolated to \(2^{50}\) without measurement.

## 6. Confirmed implementation-level scalability obstacles

### 6.1 Exact-integer threshold at \(2^{53}\)

`pi_integer_components` and `_pi_raw_and_empirical_components` explicitly reject a search region if its upper endpoint exceeds \(2^{53}\). The stated rationale is that consecutive integers above this limit cannot always be distinguished by `float64`. At \(n=2^{51}\), the anticipated magnitude of \(p_n\) is approximately \(8.55\times10^{16}\), exceeding \(2^{53}\). Thus successful calibration does not ensure that these extraction routines can return a final interval.

### 6.2 Complexity proportional to search width

In `_pi_raw_and_empirical_components`, the code contains a loop equivalent to:

```python
for x in range(low, high + 1):
    pi_lo, pi_hi = pi_cloud_at(float(x))
```

The scan therefore costs at least \(\Omega(U-L+1)\) cloud evaluations, potentially more after snapshot lookup and branch evaluation. A sufficiently wide actual search region is an algorithmic bottleneck that cannot be resolved merely by accelerating `primecount` calls.

### 6.3 Mode F cost at very large training indices

`calibration_from_pattern` times one `nth_prime` call at \(N\), then chooses a grid of 500–20,000 points using `GRID_TIME_BUDGET=900` seconds. This is an estimate based on a single call; mean call costs may vary across the grid. Ridge fitting, PI calibration, snapshot generation, and component extraction add further work. Therefore 900 seconds is **not** a global runtime upper bound.

### 6.4 Floating-point representation of candidate positions

`jump_predict` converts `n_target` to `float` in `inv_R_scalar(float(n_target))`, and candidate arrays use `np.float64`. At \(x\approx8.55\times10^{16}\), adjacent representable `float64` numbers are separated by approximately 16. Recovering integer endpoints more precisely requires separately justified numerical handling.

### 6.5 Storage versus computational complexity

`gen_primes_compact` allocates at least \(N-1\) bytes for its gap array, excluding auxiliary arrays. Hence A/B are impractical at \(N=2^{50}\); however, this is not the same bottleneck as Mode F's large-scale cost.

## 7. Performance-improvement path without changing the model

**Stage 0 — Measurement.** Run `profile_pn_pi_stable.py` in a separate process with identical inputs. Record per-function timings, peak memory, and prediction outputs. Do not adjust model coefficients while profiling.

**Stage 1 — Conservative optimizations.** Cache only truly identical operations; remove repeated evaluations only when inputs and outputs are demonstrably equivalent; investigate vectorization without changing sorting, pruning predicates, or threshold rules. Require output equivalence, not merely numerical proximity.

**Stage 2 — Large-range exact component extraction.** Develop an extractor that avoids testing every integer, but prove that it cannot skip any accepted/rejected transition. Interval arithmetic or subdivision with validated local bounds on PI cloud functions may be investigated, provided each subregion has sound enclosure bounds. Sparse sampling without guarantees is not equivalent to exhaustive scanning.

**Stage 3 — Precision beyond `float64`.** Preserve integer endpoints with Python arbitrary-precision integers or multiprecision arithmetic where necessary, distinguishing continuous function evaluation from discrete integer decisions. Merely removing the \(2^{53}\) guard is unacceptable.

**Stage 4 — Calibration scaling.** Investigate safe reuse of valid calibration caches and parallel execution of independent `primecount` queries, subject to total-memory limits and exact comparison of calibration outputs where possible. Cache validity must incorporate training size, mode, and fitting settings.

## 8. Acceptance criteria for any optimization

1. Preserve all original settings, coefficients, and mode thresholds.
2. Never use \(p_k\) for \(k>N\) in calibration or prediction.
3. Do not discard any admissible PI component due to approximation or interpolation.
4. Preserve pruning, branch, and anomaly outcomes, or explicitly document a mathematically necessary and separately authorized correction.
5. Compare original and optimized outputs on a suite of inputs for which the baseline completes.
6. Measure total runtime—not just calibration—and peak memory.
7. Do not present empirical containment as a theorem; report fallback cases separately.
8. Never convert a large exact integer endpoint to `float64` and assume the original integer can be recovered.

## 9. Established and unestablished findings

**Established from code inspection:** Original mode thresholds; A/B/C/F dispatch logic; the \(2^{53}\) guard; integer-by-integer PI scanning; Mode F's adaptive grid; raw versus empirical acceptance; and the theoretical fallback.

**Established from user-provided observations:** All four calibration modes completed through \(N=10^8\); a run targeting \(n=2^{51}\) with \(N=2^{50}\) reportedly required more than four hours. No full profiler trace for that large run was supplied.

**Not established:** Universal containment of \(p_n\); superiority over every published prime bound; optimal mode transition thresholds; the exact source of the four-hour runtime; successful completion of final interval extraction at \(n=2^{51}\); or a quantified speedup without changing the algorithm.

## 10. Work required for the final research paper

This draft is **not** a symbol-by-symbol or line-by-line derivation of the entire implementation. Completion requires extracting, checking, and mathematically analyzing equations and implementation assumptions in `R_vec`, `inv_R_vec`, `D_corr`, `build_pi_features`, `build_pn_features`, `fast_ridge_correction`, `compute_calibration`, `calibration_from_pattern`, `thresh_pi_at`, `jump_walk_pi`, `jump_walk_pn`, `_pi_raw_and_empirical_components`, `_joint_post_raw_geometry`, `certified_nth_prime_bounds`, and `jump_predict`. It also requires verifying the hypotheses of each theoretical bound against its primary published source and auditing floating-point error.

## Preliminary source record

1. **Local source code:** `PrimeCloud_v18_2_PN_to_PI_Original_Mode_Limits.py`, particularly the functions and constants identified above. This is an implementation reference, **not** an independent mathematical proof.
2. **User-provided benchmark records:** Training sizes \(10^4\)–\(10^8\). These are environment-specific empirical observations.
3. **Outstanding primary-source verification:** The `axler_bounds` and `certified_nth_prime_bounds` formulas must be checked against the original papers by Christian Axler on bounds for the \(n\)-th prime, including theorem numbers and precise applicability conditions. That verification has **not** yet been completed in this draft.


---

# Extended mathematical and reproducibility analysis

## 11. Reference function: exact definitions and approximation hierarchy

For real $x>2$, the implementation defines

$$\operatorname{li}_2(x)=\operatorname{Ei}(\log x)-\operatorname{Ei}(\log 2)=\int_2^x\frac{dt}{\log t}.$$

This is the *offset* logarithmic integral, not the principal-value logarithmic integral from zero. The distinction is an additive constant and matters when reproducing calibration intercepts. The derivative is $\operatorname{li}_2'(x)=1/\log x$. The Newton inverse implemented in `inv_li_vec` therefore has the update

$$x_{j+1}=x_j-\big(\operatorname{li}_2(x_j)-n\big)\log x_j.$$

The implemented reference is a **truncated, interpolated** Möbius correction:

$$R_{\rm impl}(x)=\operatorname{li}_2(x)+\widetilde D(x),\qquad D_{12}(x)=\sum_{k\in K}\frac{\mu(k)}{k}\operatorname{li}_2(x^{1/k}),$$

where $K=\{2,3,5,6,7,10,11,13,14,15,17,19\}$ and the signs $\mu(k)$ are explicitly encoded in `_MU_T`. Here $\widetilde D$ is piecewise linear interpolation of $D_{12}$ **in $t=\log x$** over 200,000 nodes between $\log(2.5)$ and $\log(10^{18})$. The omitted Möbius terms and interpolation error have **not** been enclosed rigorously; the implementation is not identical to the infinite Riemann $R$ function.

For $t_i\le t\le t_{i+1}$, write $d_i=D_{12}(e^{t_i})$ and $\theta=(t-t_i)/(t_{i+1}-t_i)$. Then

$$\widetilde D(e^t)=(1-\theta)d_i+\theta d_{i+1},\quad \frac{d}{dx}\widetilde D(x)=\frac{d_{i+1}-d_i}{(t_{i+1}-t_i)x}.$$

This derivative is valid only *inside* a cell. At grid nodes, the left and right derivatives need not agree. `inv_R_vec` uses this cellwise derivative, damped Newton steps and a Brent fallback. Its termination tolerance is tied to `float64` spacing, not an independently proven absolute error bound on the exact real root. The table is clamped by `np.interp` outside its domain, while selected entry points check the upper range. **Research implication:** any theorem involving $R_{\rm impl}$ must explicitly include both truncation and interpolation errors.

## 12. Ridge features, fitting objective, and what zeta zeros do not prove

For a zero ordinate $\gamma_k$ and $x>2$, the PI features are

$$\phi_{k,c}(x)=\frac{\sqrt{x}}{\gamma_k\log x}\cos(\gamma_k\log x),\qquad
\phi_{k,s}(x)=\frac{\sqrt{x}}{\gamma_k\log x}\sin(\gamma_k\log x).$$

There are 60 ordinates and therefore 120 columns. PN features have the same functional form but are evaluated at the inverse-reference location. `build_pi_features` and `build_pn_features` allocate their feature matrices in `float32`; `fast_ridge_correction` subsequently evaluates corrections in `float64`. These two precisions are not interchangeable, especially for large phase arguments.

For design matrix $X\in\mathbb R^{m\times120}$ and residual vector $y$, Ridge with `fit_intercept=False` minimizes

$$\widehat\beta_\lambda=\arg\min_\beta\left(\|X\beta-y\|_2^2+\lambda\|\beta\|_2^2\right).$$

Where $X^TX+\lambda I$ is invertible, its mathematical solution is

$$\widehat\beta_\lambda=(X^TX+\lambda I)^{-1}X^Ty.$$

`calibration_from_pattern` selects $\lambda$ from 30 logarithmically spaced candidates between $10^{-1}$ and $10^3$ using three-fold cross-validation on at most 15,000 sampled rows (`RandomState(42)`). This documents **how** the coefficients are fitted; it does not prove that zeta-zero features predict future prime residuals or that the fitted coefficients are unique to a theoretical explicit formula. The finite trigonometric basis is an empirical regression architecture.

The PI baseline contains

$$F_{b,PI}(x)=R_{\rm impl}(x)+a_{PI}\sqrt{\frac{x}{\log(\max(x,3))}}+b_{PI}.$$

`Fz_pi` adds the Ridge dot product. The complete branch model is more elaborate than this baseline; the baseline alone is **not** the final PI cloud. The derivative of the bias, if needed analytically for $x>3$, is

$$\frac{d}{dx}\left(a_{PI}\sqrt{\frac{x}{\log x}}\right)=\frac{a_{PI}}{2\sqrt{x\log x}}\left(1-\frac1{\log x}\right).$$

The source's `dbias_pi_dx` instead uses a symmetric finite difference with $h=\max(10^{-5}x,1)$; therefore derivative-based statements must account for its discretization error.

## 13. Distinct notions of correctness

**Arithmetic correctness** means every computed integer endpoint encloses the exact real quantity the numerical procedure claims to bound. **Algorithmic equivalence** means an optimized routine produces the same raw predicates, connected components, selection decisions, and final interval as the reference implementation. **Empirical coverage** means the output contains the known $p_n$ on a specified independent test set. **Universal containment** is the statement

$$\forall n\in\mathcal D:\quad L_{\rm model}(n)\le p_n\le U_{\rm model}(n),$$

which requires a proof valid for every $n$ in the stated domain. The first three do not imply the fourth. Likewise, claiming a universally tighter interval than a published theorem requires proving both containment and

$$U_{\rm model}(n)-L_{\rm model}(n)<U_{\rm published}(n)-L_{\rm published}(n)$$

for the entire claimed domain, with clear treatment of small $n$ and fallback cases.

A subtle distinction in the current implementation: the `RAW` PI predicate is tested at floating PN seeds; the `EMP` predicate adds $C_{\rm EMP,PI}\log^2x$ **in prime-count units**; full integer PI components containing compatible seeds are selected; and the output fills gaps between selected components. A wider final contiguous envelope is not itself a proof that its true prime index is covered. The selected components and filled gaps must be reported separately in any experiment.

## 14. Complexity model and a defensible acceleration target

Let $G(N)$ be the number of `nth_prime` queries in the C/F grid, $T_{\rm nth}(k)$ their runtime, $S$ the number of PI snapshots, $B$ the effective number of branches evaluated per snapshot, and $W=U-L+1$ the **actual integer scan width**. A descriptive cost decomposition is

$$T_{\rm total}=T_{\rm setup}+\sum_{j=1}^{G(N)}T_{\rm nth}(k_j)+T_{\rm fit}+T_{\rm walks}+T_{\rm extraction}.$$

The present exhaustive extractor has $T_{\rm extraction}=\Omega(W)$ and, if each PI-cloud query requires at most $O(B+\log S)$ operations under its data structures, a corresponding conditional upper estimate $O(W(B+\log S))$. The upper estimate is **not** a certified bound without profiling the actual lookup and branch operations. A one-call timing estimate for `GRID_TIME_BUDGET=900` seconds cannot guarantee $\sum_jT_{\rm nth}(k_j)\le900$ seconds.

A possible *research* replacement for scanning is a validated branch-and-bound interval extractor. On each integer block $I=[a,b]\cap\mathbb Z$, compute proven enclosures $\underline\ell_I\le\ell(x)\le\overline\ell_I$ and $\underline u_I\le u(x)\le\overline u_I$ and proven bounds on $M(x)$ for every $x\in I$. A block can be declared wholly accepted only if the **worst-case** inequalities imply acceptance throughout; wholly rejected only if the bounds imply rejection throughout; otherwise subdivide until a singleton is reached. Such an algorithm can preserve every component **only if the enclosures are mathematically sound**, including branch/snapshot discontinuities and floating-point rounding. No such enclosure implementation has been validated in the current study.

For $x>2^{53}$, a sound extractor must retain integer coordinates as arbitrary-precision integers and evaluate continuous functions with error bounds sufficient to decide integer membership. Arbitrary-precision integer endpoints alone do not repair a cloud function that first casts $x$ to `float64`.

## 15. Mathematical edge cases that require explicit validation

1. **Index-to-float conversion:** `inv_R_scalar(float(n))` loses unit precision for sufficiently large integer indices; a correct Python integer supplied at the CLI does not prevent this downstream loss.
2. **Endpoint rounding:** outward rounding must be applied to mathematically enclosed real endpoints, not to unvalidated floating approximations.
3. **Domain clipping:** intersecting a heuristic search with a theorem interval is safe only if the theorem interval is genuinely proven and correctly computed; clipping alone cannot certify the heuristic selection.
4. **No-seed case:** the fallback is a distinct output path; mixing it into an empirical success count would inflate reported performance.
5. **Disjointness:** reporting one filled envelope hides internal gaps, so record both original components and final integer width.
6. **Training leakage:** `nth_prime` queries at indices greater than $N$ must be forbidden in fitting, coefficient selection, anomaly calibration, and margin tuning; they may be used only by a clearly separated evaluation routine.
7. **Calibration sampling:** a deterministic grid and deterministic Ridge split improve reproducibility, but do not turn extrapolation into a bound.
8. **Scale transitions:** original A/B/C/F cutoffs are execution policy, not mathematical discontinuities in prime behavior; no experiment to $10^8$ certifies a cutoff at $2\times10^9$ or beyond.
9. **Two different Axler constants in source:** `axler_bounds` uses `11.508` in its lower formula, whereas `certified_nth_prime_bounds` uses `11.321`. These functions must not be described as numerically identical. The primary theorem and applicability range for **each** must be checked before publication.
10. **Unproven complexity at $2^{50}$:** the reported four-hour runtime is an observation, not a measured per-function profile. Assigning all of it to the grid, PI scan, or primecount would be speculation.

## 16. Experimental protocol suitable for a research submission

Use a frozen commit/hash of the reference implementation and publish the Python, NumPy, SciPy, scikit-learn, mpmath, numba and primecountpy versions, CPU, memory, and operating system. For each calibration size $N$, run A/B/C/F separately under the same constraints and record wall time, peak RSS, exact success/failure, all calibration coefficients, number of grid points, `nth_prime` call counts, branch counts, snapshot counts, raw PI component counts, empirical component counts, fallback reason, and final interval width. For holdout indices $n>N$, compute $p_n$ **only in the independent evaluation phase**. Report the numerator and denominator of empirical coverage and the width ratio relative to the rigorously evaluated theorem baseline. Run multiple repetitions and show variability, not only a best-case timing. Keep failures, empty sets, and timeouts in the denominator where applicable.

For a speed optimization, use a differential test suite: compare reference and optimized implementations at every checkpoint and on all final predicates. At small enough $W$, compare the full sets of integer members, not just endpoints. Include cases near $2^{53}$, table boundary $10^{18}$, calibration-mode thresholds, zero selected components, and multiple disconnected components. Only after correctness equivalence is demonstrated should speedup $T_{\rm old}/T_{\rm new}$ be reported.

## 17. Source provenance and publication limitations

**Primary implementation evidence:** the locally inspected `PrimeCloud_v18_2_PN_to_PI_Original_Mode_Limits.py`, notably `li_val`, `inv_li_vec`, `_D_exact`, `D_corr`, `R_vec`, `inv_R_vec`, `build_pi_features`, `build_pn_features`, `calibration_from_pattern`, `_pi_raw_and_empirical_components`, `_joint_post_raw_geometry`, `certified_nth_prime_bounds`, and the CLI mode dispatch. **Observed evidence:** the user's benchmark records through $N=10^8$ and the reported high-index runtime. **External literature requiring primary-text verification:** Christian Axler, *New Estimates for the nth Prime Number*, arXiv:1706.03651; original Dusart and Rosser–Schoenfeld inequalities; the Riemann explicit formula and Möbius inversion; numerical interval-arithmetic guarantees. A code comment or citation title is not sufficient to validate a theorem's constant, range, or implementation.

**Publication status:** This expanded study is a reproducible mathematical *audit and research manuscript*, not a claim of universal containment or a complete line-by-line formal verification. No performance improvement to the model itself is asserted or applied here.



# PART II — EXPANDED FORMAL ANALYSIS AND REPRODUCIBILITY

## 18. Research question and exact logical status

The scientific question is whether a finite, calibrated PN-to-PI candidate generator can output a reliably narrow interval containing the unknown nth prime, without consulting any prime of index greater than the training cutoff N. This is a **computational hypothesis**, not a theorem. Three logically independent claims must be tested: (i) the implementation executes to completion, (ii) its empirical interval contains the prime in an independent test, and (iii) the interval is universally guaranteed to contain the prime. Neither (i) nor (ii), nor their conjunction on finitely many examples, proves (iii). A fallback to a theorem interval must be labeled separately and must never be counted as empirical model containment.

For each test index n, record an outcome tuple (N,n,mode,calibration_hash,raw_PN_count,raw_PI_components,expanded_PI_components,selected_components,L,U,theorem_L,theorem_U,fallback,elapsed,peak_RSS). If p_n is available to the **evaluation process**, define H(n)=1 if L<=p_n<=U and 0 otherwise. On a fixed test set T, empirical coverage is sum(H(n))/|T|, with failures and missing predictions reported explicitly. Width comparison uses W_model=U-L and W_theorem=theorem_U-theorem_L; an interval of zero width needs careful interpretation of inclusive integer cardinality U-L+1.

## 19. Fundamental objects and domains

The primes p_1=2,p_2=3,... form a strictly increasing integer sequence. pi(x)=#{p prime:p<=x} for real x. Consequently pi(p_n)=n and pi(x)=n for p_n<=x<p_(n+1). A model approximation F_PI(x) need not be integer-valued or monotone. R(x) in the implementation is a computational reference assembled from logarithmic-integral/Möbius ingredients and an interpolated correction; equality R(x)=pi(x) is **not assumed**. The inverse R^{-1} is numerical and may not be a globally unique mathematical inverse unless monotonicity and domain restrictions are proved.

A prime-index training cutoff N means that prime values p_k with k>N are forbidden as input to fitting, coefficient selection, threshold tuning, or state initialization. The target index n>N is allowed as a known integer. An analytic published inequality in n is not leakage. Calling nth_prime(n) inside prediction would be leakage; doing so in a strictly separate evaluator is permissible. A result may still be contaminated by a built-in calibration table whose provenance overlaps a holdout set; independence requires tracing its generation.

## 20. Arithmetic identities versus computational approximations

For x>1 the logarithmic integral li(x) is conventionally the principal-value integral of 1/log t, with its singularity at t=1. The classical Riemann prime-counting approximation is R_classical(x)=sum_{k>=1} mu(k)/k li(x^(1/k)), where mu is the Möbius function. Truncating the sum to a finite selection K changes the function. The source uses a selected finite Möbius-term list, numerical li evaluation, and an interpolated D correction. Thus write R_impl(x), not an unqualified equality with the infinite Riemann R-function. The source has more than one declaration of certain correction globals: in Python, the later executed assignment takes effect. An audit must distinguish declared names from runtime bindings.

For a differentiable reference R with R'(x)!=0, the Newton iteration for R(x)=n is x_(j+1)=x_j-(R(x_j)-n)/R'(x_j). This identity does **not** imply convergence: convergence requires local regularity, a sufficiently good starting point, and appropriate handling of derivative sign and numerical error. Likewise a symmetric difference [f(x+h)-f(x-h)]/(2h) has truncation error O(h^2) under suitable smoothness, but roundoff and nonsmooth interpolation knots can dominate. Any assertion of an exact derivative for piecewise-linear-in-log interpolation at a knot must choose one-sided derivatives or state non-differentiability.

## 21. Calibration and regression: mathematical interpretation

For a feature matrix X in R^(m x d), response y in R^m and lambda>=0, the standard Ridge objective is min_beta ||X beta-y||_2^2+lambda||beta||_2^2. When X^T X+lambda I is invertible, beta=(X^T X+lambda I)^(-1)X^T y. This formula explains regularization in general; the actual code's intercept treatment, scaling, fit procedure, and features must be read from the exact source before identifying its numerical coefficients with this expression. Spectral features based on a finite selection of zeta zeros are model features, not a proof of the explicit formula's remainder. Fitting a median offset changes empirical centering, not guaranteed error coverage.

If C_EMP is defined from a maximum normalized residual on training observations, its coverage statement is restricted to the sampled training set and to the exact denominator and filtering policy. In particular, the form |p_k-F_PN(k)| <= C_EMP,PN log^2(p_k) on a measured set is not a theorem for k>N. Changing from full scans to grids changes what 'maximum' means. Mode A and B returned identical C_EMP values in the supplied experiments through 10^8, but this does not prove equality at all N or equality of final intervals.

## 22. Exact set-theoretic structure of the cloud

Let P_raw={x_1,...,x_m} be the finite PN candidates, and for integer x in the chosen scan domain D define raw lower/upper cloud values ell(x),u(x). Define A_raw={x in D:ell(x)<=n<=u(x)}. Let delta(x)=multiplier*C_EMP,PI*(log x)^2 under the code's normalization. Define A_emp={x in D:ell(x)-delta(x)<=n<=u(x)+delta(x)}. Assuming finite ell,u, C_EMP>=0 and multiplier>=0, A_raw subseteq A_emp **pointwise on the same domain**. This inclusion is mathematical and does not require monotonicity of the cloud. The code's actual implementation additionally rejects nonfinite values; comparisons at floating seeds and integers must not be silently conflated.

An integer component is a maximal consecutive run [a,b] intersect Z within the acceptance set. If the selected components are [a_i,b_i], the final contiguous envelope is [min_i a_i,max_i b_i], which fills all internal gaps. Its membership can be a strict superset of the selected set. If no component is selected, the code's theorem fallback is a distinct branch. Intersecting a proven theorem interval with an empirical interval can only preserve true-prime containment if the empirical interval already contains the prime; a theorem interval alone does not make the intersection safe.

## 23. The 2^53 barrier and a quantitative scaling bound

Binary64 represents every integer exactly only through 2^53. For integers larger than 2^53, converting to float can merge adjacent integers; Python's arbitrary-precision integer input does not fix a downstream float cast. The extractor explicitly raises when high>2^53. Even before that limit, the scan evaluates pi_cloud_at(float(x)) for each x, so the number of evaluations is exactly high-low+1 if it finishes. A scan width W requires Omega(W) calls regardless of how efficiently each call is coded. If W=589231908788, a hypothetical one million evaluations per second would require 589231.9 seconds (~163.7 hours); this is a conditional illustration, not a measured throughput. The user's reported four-hour end-to-end runtime has not been decomposed by profiling.

An acceleration must distinguish calibration costs from prediction, branch walking, PI cloud construction, integer acceptance, component selection, and printing. A replacement for per-integer enumeration may use interval subdivision **only** after rigorous interval enclosures of every branch and switching discontinuity are established. If the cloud predicate is discontinuous or has narrow isolated components, checking only endpoints or a coarse mesh can miss valid regions. Numerical validation requires directed rounding or a justified error budget, not simply increasing float precision without proof.

## 24. What must be proved for a universal containment theorem

One sufficient route is to prove, for every target n in a stated domain, that the true p_n lies in a search domain D_n and satisfies the empirical PI acceptance predicate; then prove the PN seeding/selection procedure retains the component containing p_n; finally prove the computed integer endpoints outwardly enclose the mathematical component. Each step is independent. In particular, even a perfect theorem bound D_n does not prove that the learned PI cloud accepts p_n, nor that the PN seed hits its component. A purely empirical C_EMP cannot bridge this logical gap without a theorem on out-of-sample residuals. Any claimed strict improvement over the narrowest published theorem interval also requires a valid literature comparison with matching domains and rounding conventions.

## 25. Experimental data supplied in this conversation

The following are **reported** calibration outcomes, not independent reruns by this manuscript. At N=100000000: A (592.2505 s, 974241792 peak RSS bytes, C_PN=37.6705212785, C_PI=1.2844031390); B (533.2362 s, 983617536 bytes, identical C values); C (94.1016 s, 299421696 bytes, C_PN=39.5797403780, C_PI=1.3652740741); F (11.3549 s, 298979328 bytes, C_PN=65.7980142293, C_PI=2.6028684970). All four completed at N=2M,5M,10M,20M,50M,100M in the user's records. These results certify neither 12 GB feasibility at the original 2B/200B/20T thresholds nor any optimal transition point. The reported ability to train at 2^50 and the reported >4-hour prediction of index 2^51 are user observations, not profiled or independently verified here.

## 26. Reproducibility, scientific integrity, and limitations

The reference source SHA-256 is **d6a80264bb4d5b0b90f84092f4f84936cc2d6ac7be50df0aa36e54f3cfb7930c**. Source listing below is included to remove ambiguity about the exact implementation. A manuscript including a program listing is not automatically a proof that every function is correct. To reproduce results, archive environment/package versions, seed states, CPU model, resource limits, invocation, output JSON, and a source hash. Keep all failed runs. For external mathematical results, cite primary theorem statements with precise assumptions and index thresholds; the source code's use of a named 'certified' routine does not independently establish those assumptions. The present manuscript is a comprehensive *implementation documentation and mathematical audit*, not a peer-reviewed proof of universal correctness.

## 27. Primary mathematical background references (not a verification of every coded constant)

1. Riemann, B. (1859), *Ueber die Anzahl der Primzahlen unter einer gegebenen Grösse*, original memoir (prime-counting explicit-formula context).
2. Möbius, A. F. (1832), *Über eine besondere Art von Umkehrung der Reihen* (Möbius inversion).
3. Rosser, J. B. and Schoenfeld, L. (1962), *Approximate formulas for some functions of prime numbers*, Illinois Journal of Mathematics 6, 64–94.
4. Dusart, P. (2010), *Estimates of some functions over primes without R.H.*, arXiv:1002.0442.
5. Axler, C. (2019), *New estimates for the nth prime number*, Journal of Integer Sequences 22, Article 19.4.2; preprint arXiv:1706.03651.
6. Higham, N. J. (2002), *Accuracy and Stability of Numerical Algorithms*, 2nd ed., SIAM (floating-point error analysis).
7. IEEE 754-2019, *IEEE Standard for Floating-Point Arithmetic* (binary64 integer exactness).

**Citation limitation:** The manuscript does not assert that the exact constants 11.321 and 11.508 in this code have been matched to the correct primary theorem, applicability range, and outward-rounding implementation. This remains a blocking publication verification task.


# PART III — SOURCE-LEVEL FUNCTION AND CLASS INDEX

The index below enumerates every top-level function and class, with exact source line ranges. The attached listing is authoritative for each parameter, local variable, return expression, and control-flow branch. The index does not imply that every behavior is mathematically certified.

- `gen_primes` — function, lines 24–49; implementation in Part IV.
- `CompactPrimes` — class, lines 73–142; implementation in Part IV.
- `gen_primes_compact` — function, lines 144–182; implementation in Part IV.
- `_NORM` — function, lines 190–199; implementation in Part IV.
- `_NORM_S` — function, lines 201–203; implementation in Part IV.
- `li_val` — function, lines 205–207; implementation in Part IV.
- `li_vec` — function, lines 209–214; implementation in Part IV.
- `inv_li_vec` — function, lines 216–224; implementation in Part IV.
- `_D_exact` — function, lines 245–249; implementation in Part IV.
- `_ensure_D` — function, lines 251–255; implementation in Part IV.
- `D_corr` — function, lines 257–259; implementation in Part IV.
- `_check_D_range` — function, lines 261–270; implementation in Part IV.
- `R_vec` — function, lines 272–273; implementation in Part IV.
- `inv_R_vec` — function, lines 275–310; implementation in Part IV.
- `_li_asymptotic_val` — function, lines 312–320; implementation in Part IV.
- `inv_li_vec_fast_bulk` — function, lines 323–356; implementation in Part IV.
- `fast_ridge_correction` — function, lines 357–370; implementation in Part IV.
- `inv_R_scalar` — function, lines 393–394; implementation in Part IV.
- `axler_bounds` — function, lines 400–421; implementation in Part IV.
- `build_pi_features` — function, lines 490–496; implementation in Part IV.
- `build_pn_features` — function, lines 498–504; implementation in Part IV.
- `bias_pi_scalar` — function, lines 506–507; implementation in Part IV.
- `dbias_pi_dx` — function, lines 509–511; implementation in Part IV.
- `Fb_pi` — function, lines 513–516; implementation in Part IV.
- `Fz_pi` — function, lines 518–520; implementation in Part IV.
- `_lower_hull_core` — function, lines 536–560; implementation in Part IV.
- `lower_hull_tsn_base` — function, lines 562–576; implementation in Part IV.
- `_sieve_stream` — function, lines 630–663; implementation in Part IV.
- `c_emp_streaming` — function, lines 665–691; implementation in Part IV.
- `sieve_grid` — function, lines 693–708; implementation in Part IV.
- `_margin_floor_from_data` — function, lines 848–866; implementation in Part IV.
- `measure_margin_params` — function, lines 872–908; implementation in Part IV.
- `_self_margin` — function, lines 914–944; implementation in Part IV.
- `predict_c_emp` — function, lines 946–976; implementation in Part IV.
- `calibration_from_pattern` — function, lines 978–1095; implementation in Part IV.
- `_pc` — function, lines 1124–1137; implementation in Part IV.
- `nth_prime_vec` — function, lines 1139–1144; implementation in Part IV.
- `c_emp_via_primecount` — function, lines 1146–1165; implementation in Part IV.
- `inv_R_vec_fast` — function, lines 1172–1182; implementation in Part IV.
- `_mp_cemp` — function, lines 1225–1247; implementation in Part IV.
- `_cemp_pair` — function, lines 1249–1273; implementation in Part IV.
- `save_calibration` — function, lines 1285–1288; implementation in Part IV.
- `load_calibration` — function, lines 1290–1299; implementation in Part IV.
- `compute_calibration` — function, lines 1302–1526; implementation in Part IV.
- `thresh_pi_at` — function, lines 1540–1559; implementation in Part IV.
- `jump_walk_pi` — function, lines 1562–1735; implementation in Part IV.
- `make_pi_cloud_fn` — function, lines 1740–1773; implementation in Part IV.
- `find_boundary` — function, lines 1775–1811; implementation in Part IV.
- `pi_integer_components` — function, lines 1818–1847; implementation in Part IV.
- `_pi_components_seeded_by_pn` — function, lines 1849–1859; implementation in Part IV.
- `_merge_integer_components` — function, lines 1865–1880; implementation in Part IV.
- `_pn_point_components` — function, lines 1883–1904; implementation in Part IV.
- `_pi_raw_and_empirical_components` — function, lines 1907–1949; implementation in Part IV.
- `_joint_post_raw_geometry` — function, lines 1952–2103; implementation in Part IV.
- `jump_walk_pn` — function, lines 2108–2238; implementation in Part IV.
- `jump_predict` — function, lines 2243–2316; implementation in Part IV.
- `safe_jump_predict_lowmem` — function, lines 2321–2342; implementation in Part IV.
- `run_pn_cloud` — function, lines 2359–2460; implementation in Part IV.
- `run_pi_walk_with_snapshots` — function, lines 2464–2562; implementation in Part IV.
- `make_pi_cloud` — function, lines 2566–2591; implementation in Part IV.
- `find_boundary_continuous_with_gap` — function, lines 2595–2627; implementation in Part IV.
- `process_range` — function, lines 2631–2727; implementation in Part IV.
- `_legal_anchor` — function, lines 2769–2777; implementation in Part IV.
- `_residual_at` — function, lines 2782–2788; implementation in Part IV.
- `_legal_anchor` — function, lines 2792–2800; implementation in Part IV.
- `measure_residual_rho` — function, lines 2804–2833; implementation in Part IV.
- `rho_at` — function, lines 2836–2842; implementation in Part IV.
- `measure_anomaly_persistence` — function, lines 2845–2887; implementation in Part IV.
- `anomaly_prob_at` — function, lines 2890–2896; implementation in Part IV.
- `anchored_point_estimate` — function, lines 2899–2915; implementation in Part IV.
- `wave_point_estimate` — function, lines 2919–2951; implementation in Part IV.
- `online_predict` — function, lines 2954–2983; implementation in Part IV.
- `best_point_estimate` — function, lines 2987–3008; implementation in Part IV.
- `auto_calibrate` — function, lines 3092–3151; implementation in Part IV.
- `_require_exact_int` — function, lines 3158–3179; implementation in Part IV.
- `_finite_correlation` — function, lines 3182–3197; implementation in Part IV.
- `_Rprime_table` — function, lines 3200–3209; implementation in Part IV.
- `_merge_lower_hull_true_x` — function, lines 3212–3221; implementation in Part IV.
- `_merge_lower_hull_numba` — function, lines 3225–3235; implementation in Part IV.
- `certified_nth_prime_bounds` — function, lines 3240–3268; implementation in Part IV.
- `_legacy_math_v14_cli` — function, lines 3280–3403; implementation in Part IV.


# PART IV — COMPLETE FROZEN REFERENCE IMPLEMENTATION (VERBATIM)

The following is an archival reproduction of the unchanged reference implementation, not a proposed modification. Every symbol and coefficient not separately discussed above is defined by its exact occurrence here. A source listing is not a mathematical derivation of empirical choices.

```python
# PrimeCloud v18.2: PN candidate generation -> PI validation -> complete selected PI component bounds (no PN/PI set intersection).
# Base: user original PN->PI model, with v17 spectral-RMS PI threshold replacing fixed 0.3.
# NOTE: the Moebius tables (_MU_T / _MU_TERMS) here have 12 terms (k <= 19), as in Last_work;
#       R36 has the full 36-term [R-FULL] list.
import numpy as np
import math, time, gc, pickle
from scipy.special import expi
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from sklearn.linear_model import Ridge, RidgeCV
from functools import lru_cache

T_START = time.time()

# Original user setting retained for source-configuration compatibility.
# IMPORTANT: its old formula required C_EMP_PI inside RAW. That path is
# intentionally inapplicable after the explicitly requested C_EMP separation.
# PI threshold v17 uses spectral-RMS Ridge scale; PI_THRESH_FACTOR is ignored.
PI_THRESH_FACTOR = 1.2

# ══════════════════════════════════════════════════════════════
# 1. Segmented sieve
# ══════════════════════════════════════════════════════════════
def gen_primes(N, segment_size=8_000_000):
    if N <= 0: return np.array([], dtype=np.int64)
    upper_est = int(N * (math.log(N) + math.log(math.log(N)) + 2)) + 5000 if N > 6 else 15
    while True:
        limit_sqrt = int(math.isqrt(upper_est)) + 1
        base_sieve = np.ones(limit_sqrt + 1, dtype=bool); base_sieve[:2] = False
        for i in range(2, int(limit_sqrt**0.5) + 1):
            if base_sieve[i]: base_sieve[i*i::i] = False
        base_primes = np.flatnonzero(base_sieve)
        result = np.empty(N, dtype=np.int64)
        filled = 0; low = 2
        while low <= upper_est and filled < N:
            high = min(low + segment_size - 1, upper_est)
            seg = np.ones(high - low + 1, dtype=bool)
            for p in base_primes:
                if p * p > high: break
                start = max(p * p, ((low + p - 1) // p) * p)
                seg[start-low::p] = False
            if low <= 1: seg[:max(0, 2-low)] = False
            seg_primes = np.flatnonzero(seg) + low
            take = min(len(seg_primes), N - filled)
            result[filled:filled+take] = seg_primes[:take]
            filled += take; low = high + 1
            del seg, seg_primes
        if filled >= N: return result
        upper_est = int(upper_est * 1.5)



# ══════════════════════════════════════════════════════════════
# [COMPACT] Compressed prime storage -- 8x memory reduction
# ══════════════════════════════════════════════════════════════
# Problem: P_ALL as int64 takes 8 bytes per prime. At TRAIN_END=10^9 that
# is **8 GB** resident for the whole run -- the cause of out-of-memory
# crashes on Colab (12 GB total), before any computation starts.
#
# Solution: store **gaps**, not values. Every prime after 2 is odd, so
# gaps between them are always even => store gap/2 in a uint8 (one byte).
# Measured: max gap/2 over the first 20 million primes = 124, and the
# largest known gap below 2.3e10 is ~464 => gap/2 = 232 < 255. uint8
# suffices with margin.
#
# Retrieval: we keep anchor points every ANCHOR_STRIDE primes (negligible
# cost); any slice is reconstructed by cumsum from the nearest anchor.
# Verified: slices are **bit-identical** to the original array.
#
# Memory at TRAIN_END=10^9:  8.0 GB  ->  1.0 GB
ANCHOR_STRIDE = 1_000_000

class CompactPrimes:
    """Lossless compact storage for the first N primes with checked gap width."""
    def __init__(self,P):
        orig=np.asarray(P)
        if orig.ndim!=1 or not np.issubdtype(orig.dtype,np.integer):
            raise ValueError('prime data must be a 1-D integer sequence')
        if np.issubdtype(orig.dtype,np.unsignedinteger) and len(orig) and np.max(orig)>np.iinfo(np.int64).max:
            raise OverflowError('prime sequence exceeds signed 64-bit range')
        P=np.asarray(P,dtype=np.int64)
        if P.ndim!=1 or (len(P) and P[0]!=2) or (len(P)>1 and P[1]!=3):
            raise ValueError('expected consecutive primes beginning 2,3')
        self.n=len(P)
        podd=P[1:]
        dd=np.diff(podd)
        if len(dd) and (np.any(dd<=0) or np.any(dd%2)):
            raise ValueError('invalid gap sequence')
        d=dd//2
        max_gap=int(d.max()) if len(d) else 0
        typ=np.uint8 if max_gap<=255 else np.uint16 if max_gap<=65535 else np.uint32
        self.g=np.empty(len(podd),dtype=typ)
        if len(podd):self.g[0]=0
        if len(d):self.g[1:]=d
        self.odd0=3
        self.anchors=np.asarray(podd[::ANCHOR_STRIDE],dtype=np.int64).copy()
    def __len__(self):return self.n
    def _range(self,lo,hi):
        lo=max(0,int(lo));hi=min(self.n,int(hi))
        if hi<=lo:return np.empty(0,dtype=np.int64)
        head=np.array([2],dtype=np.int64) if lo==0 else np.empty(0,dtype=np.int64)
        j0=max(lo-1,0);j1=hi-1
        if j1<=j0:return head
        ai=j0//ANCHOR_STRIDE
        base=int(self.anchors[ai]);jstart=ai*ANCHOR_STRIDE
        d=self.g[jstart+1:j1].astype(np.int64)
        out=np.empty(j1-jstart,dtype=np.int64);out[0]=base
        if len(d):out[1:]=base+2*np.cumsum(d)
        tail=out[j0-jstart:]
        return np.concatenate((head,tail)) if len(head) else tail
    def __getitem__(self,k):
        if isinstance(k,slice):
            start,stop,step=k.indices(self.n)
            if step<0:return self[np.fromiter(range(start,stop,step),dtype=np.int64)]
            return self._range(start,stop)[::step]
        if isinstance(k,(int,np.integer)):
            idx=int(k);idx=idx+self.n if idx<0 else idx
            if not 0<=idx<self.n:raise IndexError('prime index out of bounds')
            return int(self._range(idx,idx+1)[0])
        idx=np.asarray(k)
        if idx.ndim!=1:raise IndexError('only 1-d fancy indexing is supported')
        if idx.dtype==bool:
            if len(idx)!=self.n:raise IndexError('boolean mask wrong length')
            idx=np.flatnonzero(idx)
        elif not np.issubdtype(idx.dtype,np.integer):
            raise IndexError('indices must be integer or boolean')
        idx=np.asarray(idx,dtype=np.int64)
        if np.any(idx>=self.n) or np.any(idx< -self.n):raise IndexError('out of bounds')
        idx=np.where(idx<0,idx+self.n,idx)
        out=np.empty(len(idx),dtype=np.int64)
        order=np.argsort(idx,kind='stable');seq=idx[order]
        i=0
        while i<len(seq):
            block=max(int(seq[i])-1,0)//ANCHOR_STRIDE
            j=i+1
            while j<len(seq) and max(int(seq[j])-1,0)//ANCHOR_STRIDE==block:j+=1
            lo=int(seq[i]);hi=int(seq[j-1])+1
            part=self._range(lo,hi)
            out[order[i:j]]=part[seq[i:j]-lo]
            i=j
        return out
    def astype(self,dt):return self[:].astype(dt)

def gen_primes_compact(N,segment_size=8_000_000):
    """Stream-first gap compression, safe even for segment_size=1.

    Promotion uint8->uint16->uint32 prevents a silent overflow if future
    prime gaps exceed the historical 510-unit cap.
    """
    N=_require_exact_int(N,'N',1)
    segment_size=_require_exact_int(segment_size,'segment_size',1)
    if N==1:return CompactPrimes(np.array([2],dtype=np.int64))
    g=np.empty(N-1,dtype=np.uint8)
    g[0]=0;anchors=[];nodd=0;previous=None
    for _,sp0 in _sieve_stream(N,segment_size):
        sp=sp0[sp0!=2]
        if not len(sp):continue
        take=len(sp)
        d=np.empty(take,dtype=np.int64)
        if previous is None:
            if int(sp[0])!=3:raise RuntimeError('expected 3 after 2')
            d[0]=0
        else:
            d[0]=(int(sp[0])-previous)//2
        if take>1:d[1:]=np.diff(sp)//2
        if np.any(d<0) or (previous is not None and (int(sp[0])-previous)%2):
            raise ValueError('invalid even prime-gap encoding')
        maximum=int(d.max())
        if maximum > np.iinfo(np.uint32).max:
            raise OverflowError('half-prime-gap exceeds uint32 representable range')
        if maximum>np.iinfo(g.dtype).max:
            typ=np.uint16 if maximum<=65535 else np.uint32
            g=g.astype(typ)
        g[nodd:nodd+take]=d
        k=((nodd+ANCHOR_STRIDE-1)//ANCHOR_STRIDE)*ANCHOR_STRIDE if nodd else 0
        while k<nodd+take:
            anchors.append(int(sp[k-nodd]));k+=ANCHOR_STRIDE
        previous=int(sp[-1]);nodd+=take
    if nodd!=N-1:raise RuntimeError('compressed sieve count mismatch')
    obj=CompactPrimes.__new__(CompactPrimes)
    obj.n=N;obj.odd0=3;obj.g=g;obj.anchors=np.asarray(anchors,dtype=np.int64)
    return obj


# ══════════════════════════════════════════════════════════════
# 2. li / inv_li
# ══════════════════════════════════════════════════════════════
_EI_LN2 = float(expi(math.log(2)))

def _NORM(a):
    """[NORM] Normalisation = ln^2(n), **kept as is**.
    Replacing it with sqrt(n)/ln(n) was tried (dispersion 1.08x vs 1.57x --
    more uniform *inside* the training range). But measurement on actual
    targets showed it **widens the interval 1.8x** (13,914 vs 7,846),
    because it grows faster under extrapolation beyond training.
    Lesson: lower dispersion **inside** training does not imply better
    extrapolation **outside** it. ln^2 is retained."""
    a = np.maximum(np.asarray(a, dtype=float), 3.0)
    return np.log(a)**2

def _NORM_S(v):
    v = max(float(v), 3.0)
    return math.log(v)**2

def li_val(x):
    if x <= 2: return 0.0
    return float(expi(math.log(x)) - _EI_LN2)

def li_vec(x_arr):
    x_arr = np.asarray(x_arr, dtype=float)
    out = np.zeros_like(x_arr)
    mask = x_arr > 2.0
    out[mask] = expi(np.log(x_arr[mask])) - _EI_LN2
    return out

def inv_li_vec(ns, iters=6):
    ns = np.asarray(ns, dtype=float)
    lns = np.log(np.maximum(ns, 2.))
    x = np.maximum(ns * (lns + np.log(np.maximum(lns, 1e-9))), 2.5)
    for _ in range(iters):
        lx = np.log(x)
        li_x = expi(lx) - _EI_LN2
        x = np.maximum(x - (li_x - ns) * lx, 2.5)
    return x


# ══════════════════════════════════════════════════════════════
# [UPGRADE-R] Riemann R function instead of li
# ══════════════════════════════════════════════════════════════
# R(x) = sum mu(k)/k * li(x^(1/k)) -- the most accurate known analytic
# approximation to pi(x).
# Why it beats li: li approximates J(x) = pi(x) + (1/2)pi(sqrt x) +
# (1/3)pi(cbrt x) + ..., i.e. it counts prime **powers** as well. Moebius
# inversion removes that excess.
# Cost: R needs ~7 expi calls instead of one. Solution: D(x)=R(x)-li(x) is
# a very smooth function -> tabulate it once on a fixed logarithmic grid
# and interpolate (error ~8e-6).
# The first Newton step is then **free**: if li(x0)=n exactly, then
# R(x0)-n = D(x0) directly.
_MU_T = [(2,-1),(3,-1),(5,-1),(6,1),(7,-1),(10,1),(11,-1),(13,-1),
         (14,1),(15,1),(17,-1),(19,-1)]
_D_LOGX = None; _D_VALS = None
_D_XMAX = 1e18; _D_NPTS = 200_000

def _D_exact(x):
    x = np.asarray(x, dtype=float); s = np.zeros_like(x)
    for k, mu in _MU_T:
        s += (mu/k) * li_vec(np.maximum(x, 2.0001)**(1.0/k))
    return s

def _ensure_D():
    global _D_LOGX, _D_VALS
    if _D_LOGX is None:
        lg = np.linspace(math.log(2.5), math.log(_D_XMAX), _D_NPTS)
        _D_LOGX = lg; _D_VALS = _D_exact(np.exp(lg))

def D_corr(x):
    x = np.asarray(x, dtype=float); _ensure_D()
    return np.interp(np.log(np.maximum(x, 2.5)), _D_LOGX, _D_VALS)

def _check_D_range(x, where):
    # [GUARD] D_corr's table ends at _D_XMAX and np.interp silently clamps beyond it, so R(x)
    # would be wrong with no warning. The check sits at the entry points that turn an index n
    # into an x (R^-1) or set the end of a pi walk -- NOT inside D_corr: the pi-side boundary
    # search (find_boundary_continuous_with_gap) legitimately probes absurdly large x
    # (up to ~1e88) where the clamp is harmless, and a raise there would break normal runs.
    x = np.asarray(x, dtype=float)
    if x.size and float(np.max(x)) > _D_XMAX * (1.0 + 1e-9):
        raise ValueError(f"{where}: x={float(np.max(x)):.6g} is beyond the D_corr table limit "
                         f"{_D_XMAX:.0e}; raise _D_XMAX (and re-validate R / R^-1 there) before using larger x.")

def R_vec(x):
    return li_vec(x) + D_corr(x)

def inv_R_vec(ns):
    """Solve R_truncated(x)=n with its actual piecewise interpolation derivative.

    Newton correction is damped to remain in the table. A scalar Brent solver
    handles any unconverged positions. Relative error is limited by float64 ULP.
    """
    from scipy.optimize import brentq
    n=np.asarray(ns,dtype=np.float64)
    if np.any(~np.isfinite(n)) or np.any(n<0.0):
        raise ValueError('inverse R needs finite nonnegative indices')
    original_shape=n.shape
    nf=np.atleast_1d(n).reshape(-1)
    if np.any(nf>R_vec(np.array([_D_XMAX]))[0]):
        raise ValueError('requested n exceeds R correction table range')
    x=np.maximum(np.asarray(inv_li_vec(nf),dtype=float),2.0)
    x=np.minimum(x,_D_XMAX)
    _ensure_D()
    for _ in range(12):
        lg=np.log(x)
        ids=np.clip(np.searchsorted(_D_LOGX,lg,side='right')-1,0,len(_D_LOGX)-2)
        d_c=(_D_VALS[ids+1]-_D_VALS[ids])/(_D_LOGX[ids+1]-_D_LOGX[ids])
        derivative=1.0/lg+d_c/x
        f=R_vec(x)-nf
        tol=np.maximum(1e-11,16*np.maximum(np.spacing(nf),np.abs(derivative)*np.spacing(x)))
        if np.all(np.abs(f)<=tol):break
        new=x-f/derivative
        new=np.where(np.isfinite(new)&(new>=2.0)&(new<=_D_XMAX),new,(x+np.where(f>0,2.0,_D_XMAX))/2.0)
        x=np.maximum(2.0,np.minimum(new,_D_XMAX))
    residual=np.abs(R_vec(x)-nf)
    tol=np.maximum(1e-9,64*np.maximum(np.spacing(nf), np.spacing(x)/np.log(np.maximum(x,2.1))))
    bad=np.flatnonzero(residual>tol)
    for i in bad:

        x[i]=brentq(lambda xx:float(R_vec(np.array([xx]))[0]-nf[i]),
                    2.0,_D_XMAX,xtol=1e-12,rtol=9e-16,maxiter=120)
    return x.reshape(original_shape) if original_shape else float(x[0])

def _li_asymptotic_val(x, K=10):
    """Asymptotic series approximation to li(x) -- cheap (no expi), but
    accurate only for large x."""
    lx = np.log(np.maximum(x, 2.))
    s = np.ones_like(x); term = np.ones_like(x)
    for k in range(1, K+1):
        term = term * k / lx
        s = s + term
    return (x/lx) * s


def inv_li_vec_fast_bulk(ns, small_thresh=100_000.0, cheap_iters=4, exact_iters=1):
    """
    [SPEED-FIX-2] A faster substitute for inv_li_vec, intended specifically
    for the full C_EMP scan over millions of points. Idea: the same correct
    Newton step (same derivative 1/ln(x)), but using the cheap asymptotic
    li(x) series (no expi) for the early iterations, followed by a single
    correction with the true expi to remove the series bias.
    **Warning**: the asymptotic series breaks down for small x (verified:
    error reaches 89% at n=10) -- so any n below small_thresh is computed
    by the ordinary exact method unchanged (cheap anyway, being few points).
    Verified: max error ~1e-12 across a mixed range (small+large), far more
    accurate than required, and ~2.2x faster than inv_li_vec(iters=4) on
    large inputs.
    """
    ns = np.asarray(ns, dtype=np.float64)
    out = np.empty_like(ns)
    small_mask = ns < small_thresh
    if small_mask.any():
        out[small_mask] = inv_li_vec(ns[small_mask], iters=6)
    big_mask = ~small_mask
    if big_mask.any():
        nb = ns[big_mask]
        lns = np.log(np.maximum(nb, 2.))
        x = np.maximum(nb * (lns + np.log(np.maximum(lns, 1e-9))), 2.5)
        for _ in range(cheap_iters):
            lx = np.log(x)
            li_x = _li_asymptotic_val(x, K=10)
            x = np.maximum(x - (li_x - nb) * lx, 2.5)
        for _ in range(exact_iters):
            lx = np.log(x)
            li_x = expi(lx) - _EI_LN2
            x = np.maximum(x - (li_x - nb) * lx, 2.5)
        out[big_mask] = x
    return out
def fast_ridge_correction(x_arr,G60,ridge_coef):
    """Float64 evaluation consistent with double-precision prediction.

    The old float32 argument and trigonometric phases lose resolution at
    x~10^16-10^18; measured low-scale errors cannot justify extrapolation.
    """
    x=np.atleast_1d(np.asarray(x_arr,dtype=np.float64))
    lx=np.log(np.maximum(x,2.0));base=np.sqrt(np.maximum(x,2.0))/lx
    coefficients=np.asarray(ridge_coef,dtype=np.float64)
    result=np.zeros(x.shape,dtype=np.float64)
    for k,g in enumerate(G60):
        angle=float(g)*lx
        result+=(base/float(g))*(coefficients[2*k]*np.cos(angle)+coefficients[2*k+1]*np.sin(angle))
    return result


# ══════════════════════════════════════════════════════════════════════
# [UPGRADE-R] Analytic base: Riemann R function instead of li
# ══════════════════════════════════════════════════════════════════════
# R(x) = sum_{k>=1} mu(k)/k * li(x^(1/k)) -- the most accurate known
# analytic approximation to pi(x).
# Measured against the true primes: R^-1 is closer to truth than li^-1 by
# a factor of 1.4x to 62x depending on scale. Using it as the base shrinks
# the residuals, hence shrinks C_EMP, hence narrows the final interval.
#
# Practical problem: R(x) needs ~7 expi calls per point instead of one
# (7x slower over 500 million points = unacceptable). Solution:
# D(x) = R(x) - li(x) is a very smooth, slowly varying function, so we
# tabulate it once on a logarithmic grid and use np.interp.
# Verified: interpolation error 1e-4 (negligible against residuals in the
# thousands), and 14x faster.
_MU_TERMS = [(2,-1),(3,-1),(5,-1),(6,1),(7,-1),(10,1),(11,-1),(13,-1),
             (14,1),(15,1),(17,-1),(19,-1)]
_D_LOGX = None; _D_VALS = None


def inv_R_scalar(n):
    return float(inv_R_vec(np.array([float(n)]))[0])


# ══════════════════════════════════════════════════════════════
# 3. Axler
# ══════════════════════════════════════════════════════════════
def axler_bounds(n):
    """[FIX-1] The constant in the upper bound was 0 instead of 10.667.
    Reference: Christian Axler, "New Estimates for the nth Prime Number",
    arXiv:1706.03651 (J. Integer Sequences 22, 2019), Theorem 1.1.
    Measured impact at n=10,000,000: width 221,484 -> 16,186 (13.7x narrower).
    The lower bound already matched Theorem 1.2 (n>=2, unconditional) -- unchanged.
    Three tiers according to each theorem's validity range:
      n >= 46,254,381 : Theorem 1.1 (tightest proven)
      688,383 <= n    : Dusart (1.10)
      below          : Rosser-Schoenfeld 1962"""
    y = math.log(n); w = math.log(y)
    tr = (w - 2.) / y
    if n >= 46_254_381:
        cr_u = (w*w - 6.*w + 10.667) / (2.*y*y)
        upper = n * (y + w - 1. + tr - cr_u)
    elif n >= 688_383:
        upper = n * (y + w - 1. + tr)
    else:
        upper = n * (y + w - 0.5)
    cr_l = (w*w - 6.*w + 11.508) / (2.*y*y)
    lower = n * (y + w - 1. + tr - cr_l)
    return lower, upper

# ══════════════════════════════════════════════════════════════
# 4. Correct zeta zeros
# ══════════════════════════════════════════════════════════════
G60_PI = np.array([
    14.1347251417346937904572519835625,
    21.0220396387715549926284795938969,
    25.0108575801456887632137909925628,
    30.4248761258595132103118975305840,
    32.9350615877391896906623689640747,
    37.5861781588256712572177634807053,
    40.9187190121474951873981269146334,
    43.3270732809149995194961221654068,
    48.0051508811671597279424727494277,
    49.7738324776723021819167846785638,
    52.9703214777144606441472966088808,
    56.4462476970633948043677594767060,
    59.3470440026023530796536486749922,
    60.8317785246098098442599018245241,
    65.1125440480816066608750542531836,
    67.0798105294941737144788288965221,
    69.5464017111739792529268575265547,
    72.0671576744819075825221079698261,
    75.7046906990839331683269167620305,
    77.1448400688748053726826648563047,
    79.3373750202493679227635928771161,
    82.9103808540860301831648374947706,
    84.7354929805170501057353112068276,
    87.4252746131252294065316678509191,
    88.8091112076344654236823480795095,
    92.4918992705584842962597252418105,
    94.6513440405198869665979258152080,
    95.8706342282453097587410292192467,
    98.8311942181936922333244201386224,
    101.3178510057313912287854479402924,
    103.7255380404783394163984081086952,
    105.4466230523260944936708324141119,
    107.1686111842764075151233519630860,
    111.0295355431696745246564503099445,
    111.8746591769926370856120787167707,
    114.3202209154527127658909372761910,
    116.2266803208575543821608043120647,
    118.7907828659762173229791397026999,
    121.3701250024206459189455329704998,
    122.9468292935525882008174603307700,
    124.2568185543457671847320079661301,
    127.5166838795964951242793237669060,
    129.5787041999560509857680339061800,
    131.0876885309326567235663724615015,
    133.4977372029975864501304920426407,
    134.7565097533738713313260641571699,
    138.1160420545334432001915551902824,
    139.7362089521213889504500465233824,
    141.1237074040211237619403538184753,
    143.1118458076206327394051238689139,
    146.0009824867655185474025075964246,
    147.4227653425596020495211850104316,
    150.0535204207848803514324672369594,
    150.9252576122414667618525246783058,
    153.0246938111988961982565442551854,
    156.1129092942378675697501893101691,
    157.5975918175940598875305031584988,
    158.8499881714204987241749947755403,
    161.1889641375960275194373441293695,
    163.0307096871819872433110390006880
])
G60_PN = G60_PI.copy()

def build_pi_features(x_arr):
    x_arr = np.atleast_1d(x_arr); lx = np.log(np.maximum(x_arr, 2.)); sq = np.sqrt(np.maximum(x_arr, 2.))
    X = np.zeros((len(x_arr), 2*len(G60_PI)), dtype=np.float32)
    for k, g in enumerate(G60_PI):
        A = sq/(lx*g); th = g*lx
        X[:, 2*k] = A*np.cos(th); X[:, 2*k+1] = A*np.sin(th)
    return X

def build_pn_features(li_arr):
    li_arr = np.atleast_1d(li_arr); lg = np.log(li_arr); sq = np.sqrt(li_arr)
    X = np.zeros((len(li_arr), 2*len(G60_PN)), dtype=np.float32)
    for k, g in enumerate(G60_PN):
        A = sq/(lg*g); th = g*lg
        X[:, 2*k] = A*np.cos(th); X[:, 2*k+1] = A*np.sin(th)
    return X

def bias_pi_scalar(x, A_BIAS_PI, B_BIAS_PI):
    return A_BIAS_PI*math.sqrt(x/math.log(max(x,3.))) + B_BIAS_PI

def dbias_pi_dx(x, A_BIAS_PI, B_BIAS_PI):
    h = max(x*1e-5, 1.0)
    return (bias_pi_scalar(x+h,A_BIAS_PI,B_BIAS_PI)-bias_pi_scalar(x-h,A_BIAS_PI,B_BIAS_PI))/(2*h)

def Fb_pi(x, A_BIAS_PI, B_BIAS_PI):
    # [UPGRADE-R] R is required here too, otherwise an inconsistency of size
    # D(x) arises against the calibration
    return li_val(x) + float(D_corr(np.array([float(x)]))[0]) + bias_pi_scalar(x, A_BIAS_PI, B_BIAS_PI)

def Fz_pi(x, A_BIAS_PI, B_BIAS_PI, ridge_coef):
    f = build_pi_features(np.array([x]))[0]
    return Fb_pi(x, A_BIAS_PI, B_BIAS_PI) + float(np.dot(f, ridge_coef))

# ══════════════════════════════════════════════════════════════
# 5. lower_hull_tsn_base
# ══════════════════════════════════════════════════════════════
try:
    from numba import njit
    _HAVE_NUMBA = True
except ImportError:
    _HAVE_NUMBA = False
    def njit(*a, **k):
        if len(a)==1 and callable(a[0]): return a[0]
        return lambda f: f

# Runtime fix: do not load stale Numba cache from dynamic Colab/RunPy modules.
@njit(cache=False)
def _lower_hull_core(p_arr_int):
    """[FIX-4] The identical monotone-chain algorithm, but:
    - n is generated internally (not stored as a full float64 array)
    - P_ALL is passed as is (int64, no float64 copy)
    - the hull grows dynamically instead of pre-allocating N elements
    The old version allocated ~20 GB at TRAIN_END=500M (an actual crash).
    The result is **bit-identical** (zero difference), and ~67x faster."""
    N = p_arr_int.shape[0]
    cap = 4096 if N > 4096 else max(N, 1)
    hx = np.empty(cap, dtype=np.float64); hy = np.empty(cap, dtype=np.float64)
    top = -1
    for i in range(N):
        x = np.float64(i + 1); y = np.float64(p_arr_int[i])
        while top >= 1:
            x1 = hx[top-1]; y1 = hy[top-1]; x2 = hx[top]; y2 = hy[top]
            if (x2-x1)*(y-y1) - (y2-y1)*(x-x1) <= 0: top -= 1
            else: break
        top += 1
        if top >= hx.shape[0]:
            nc = hx.shape[0]*2
            nhx = np.empty(nc, dtype=np.float64); nhy = np.empty(nc, dtype=np.float64)
            for j in range(hx.shape[0]): nhx[j] = hx[j]; nhy[j] = hy[j]
            hx = nhx; hy = nhy
        hx[top] = x; hy[top] = y
    return hx[:top+1], hy[:top+1]

def lower_hull_tsn_base(TRAIN_END,P_ALL,chunk=20_000_000):
    """Exact chunked monotone-chain lower hull using preserved global indices."""
    if TRAIN_END<2 or TRAIN_END>len(P_ALL): raise ValueError('invalid TRAIN_END')
    hx=np.empty(0,np.float64);hy=np.empty(0,np.float64)
    for cs in range(0,int(TRAIN_END),max(1,int(chunk))):
        ce=min(int(TRAIN_END),cs+int(chunk))
        bp=np.asarray(P_ALL[cs:ce],dtype=np.int64)
        local_x,local_y=_lower_hull_core(bp)
        local_x+=cs
        hx,hy=_merge_lower_hull_numba(hx,hy,local_x,local_y)
        del bp
    slopes=np.diff(hy)/np.diff(hx)
    fi=np.floor(slopes).astype(np.int64)
    hit=np.flatnonzero(np.diff(fi)>=2)
    return float(slopes[hit[0]] if len(hit) else slopes[-1])

# ══════════════════════════════════════════════════════════════
# 6. Calibration
# ══════════════════════════════════════════════════════════════
C_EMP_CHUNK = 2_000_000     # [PORT-1] C_EMP scan chunk size (memory bound)
C_EMP_MIN_ARG = 10_000      # critical C_EMP cutoff (see [FIX-2]) -- single source
BOUNDARY_SEEDS = 5          # [PORT-10] number of seeds for the boundary search.
# The pi(x) cloud has gaps and the search logic steps over them, so the
# final bound used to depend on the seed

# ══════════════════════════════════════════════════════════════
# [PRIMECOUNT] Alternative source for p(n) -- no sieving, no storage
# ══════════════════════════════════════════════════════════════
# Calibration needs p(n) at points inside the training range. The current
# approach sieves TRAIN_END numbers and stores them -- O(n) memory, which
# is the model's ceiling (~10^9 on 12 GB).
#
# primecount implements Meissel-Lehmer/Gourdon: it computes p(n) by
# combinatorial decomposition with O(x^(1/3)) memory and no storage.
# Actual measurement:
#     nth_prime(10^9) = 22,801,763,489  in 0.002s
#     100,000 calls   ~ 60 seconds
#     agreement with sieve: 100% (tested at 3 well-separated points)
#
# And C_EMP does not need every point -- measured on a sample:
#      30,000 points -> 95.4% of the full value
#     100,000 points -> 98.3%   (and at TRAIN_END=60M: 98.26% -- stable)
# The shortfall is covered by PC_MARGIN.
#
# Impact at TRAIN_END=10^9:
#     memory          : 1 GB    ->  ~0
#     calibration time: 3,712s  ->  ~60s   (60x faster)
#     ceiling         : 10^9    ->  primecount's own limit

# ══════════════════════════════════════════════════════════════
# [STREAM] Full C_EMP scan without storage -- full coverage, bounded memory
# ══════════════════════════════════════════════════════════════
# Problem: the only way to guarantee C_EMP is to inspect **every** point.
# A full sieve does that, but stores P_ALL (1-8 GB at 10^9), which is the
# model's ceiling. Sampling breaks the ceiling but loses the guarantee
# (a 2-5% shortfall covered by an estimated margin).
#
# Solution: sieve **in chunks** and compute C_EMP on each chunk as soon as
# it is generated, then discard it. Every point is genuinely inspected =>
# **full coverage, zero shortfall, no estimated margin**.
# Memory equals the chunk size only (~8 MB) regardless of TRAIN_END.
#
# Measured at TRAIN_END=20M: a full sieve pass takes 2.3s with an 8 MB
# peak per chunk (versus 2.9s + 160 MB for full storage).
#
# Two passes are needed: the first collects the calibration grid (to fit
# the coefficients), the second computes C_EMP -- because C_EMP requires
# ridge_coef and median to be already known.
def _sieve_stream(TRAIN_END,segment_size=8_000_000):
    """Yield each of the first N primes ONCE in consecutive, bounded chunks.

    The original code restarted from index zero when the estimated ceiling
    failed, which would silently duplicate already yielded data. We reject
    invalid upper estimates rather than silently replaying yielded indices.
    """
    N=int(TRAIN_END);segment_size=int(segment_size)
    if N<1 or segment_size<1:raise ValueError('N and segment size must be positive')
    up=(int(math.ceil(N*(math.log(N)+math.log(math.log(N))+2)))+5000) if N>6 else 15
    if up>np.iinfo(np.int64).max:
        raise OverflowError('prime bound exceeds int64; use arbitrary precision method')
    rt=math.isqrt(up)+1
    base=np.ones(rt+1,dtype=bool);base[:2]=False
    for i in range(2,math.isqrt(rt)+1):
        if base[i]:base[i*i::i]=False
    primes=np.flatnonzero(base)
    cnt=0
    for low in range(2,up+1,segment_size):
        if cnt>=N:break
        high=min(low+segment_size-1,up)
        seg=np.ones(high-low+1,dtype=bool)
        for p in primes:
            p=int(p)
            if p*p>high:break
            st=max(p*p,((low+p-1)//p)*p)
            seg[st-low::p]=False
        part=np.flatnonzero(seg).astype(np.int64)+low
        k=min(len(part),N-cnt)
        if k:
            yield cnt,part[:k]
            cnt+=k
    if cnt!=N:
        raise RuntimeError(f'prime ceiling insufficient: yielded {cnt} of {N}')

def c_emp_streaming(TRAIN_END, ridge_coef, median_corr, side='pn',
                    A_BIAS=None, B_BIAS=None, segment_size=8_000_000):
    """C_EMP over **every** point, without storage. Full coverage, bounded
    memory."""
    best = 0.0
    for i0, sp in _sieve_stream(TRAIN_END, segment_size):
        nf = np.arange(i0+1, i0+1+len(sp), dtype=np.float64)
        if side == 'pn':
            bx = inv_R_vec_fast(nf)
            F = bx + fast_ridge_correction(bx, G60_PN, ridge_coef) + median_corr
            k = nf >= C_EMP_MIN_ARG
            if k.any():
                m = float(np.max(np.abs(sp[k].astype(np.float64) - F[k]) / _NORM(nf[k])))
                if m > best: best = m
            del bx, F, k
        else:
            xv = sp.astype(np.float64)
            F = (R_vec(xv) + (A_BIAS*np.sqrt(xv/np.log(np.maximum(xv,3.)))+B_BIAS)
                 + fast_ridge_correction(xv, G60_PI, ridge_coef) + median_corr)
            k = xv >= C_EMP_MIN_ARG
            if k.any():
                m = float(np.max(np.abs(nf[k] - F[k]) / _NORM(xv[k])))
                if m > best: best = m
            del xv, F, k
        del nf, sp
    gc.collect()
    return best

def sieve_grid(TRAIN_END, grid_n, segment_size=8_000_000):
    """Collects a calibration grid (spread indices plus their values) in a
    single pass, without full storage."""
    idx = np.unique(np.linspace(0, TRAIN_END-1, grid_n).astype(np.int64))
    out = np.empty(len(idx), dtype=np.int64); w = 0
    for i0, sp in _sieve_stream(TRAIN_END, segment_size):
        hi = i0 + len(sp)
        j = w
        while j < len(idx) and idx[j] < hi: j += 1
        if j > w:
            out[w:j] = sp[idx[w:j] - i0]
            w = j
        del sp
        if w >= len(idx): break
    gc.collect()
    return idx, out



# ══════════════════════════════════════════════════════════════
# [SCALES] Reference calibration series -- embedded in the code
# ══════════════════════════════════════════════════════════════
# C_EMP values measured directly at 13 scales (doublings from 5K to
# 20.48M). Used to predict larger scales (mode F) when direct calibration
# is infeasible.
# Embedded here rather than in an external file: this is a fixed part of
# the method, and its absence made mode F fail completely at huge scales
# -- the very case it exists for.
# Add your own new direct calibrations here to improve prediction
# accuracy.
BUILTIN_SCALES = [5_000, 10_000, 20_000, 40_000, 80_000, 160_000, 320_000,
                  640_000, 1_280_000, 2_560_000, 5_120_000, 10_240_000, 20_480_000,
                  500_000_000_000]
BUILTIN_C_EMP_PN = [0.0000, 0.7246, 1.2808, 1.8939, 2.6032, 3.3351, 3.8619,
                    4.2948, 7.3161, 10.1835, 11.2533, 14.5728, 18.9255,
                    1495.315]
BUILTIN_C_EMP_PI = [0.0464, 0.0679, 0.0693, 0.0993, 0.1215, 0.1491, 0.1820,
                    0.1969, 0.3005, 0.4053, 0.4422, 0.5354, 0.7463,
                    38.3170]
# The last point (5x10^11) was measured directly with mode C on Colab:
#   C_EMP_PN=1644.847 and C_EMP_PI=42.1487 were printed **after** applying
#   PC_MARGIN=1.10, so the raw values stored here are printed / 1.10.
#   Containment verified: p(10^12)=29,996,224,275,833 lies inside
#   [29,996,221,520,875 , 29,996,226,242,554] with a 41.7% edge margin,
#   and 116.65x narrower than Axler.
#   The exponent between it and 20.48M is 0.4325 -- consistent with the
#   0.377 seen within the series.


# ══════════════════════════════════════════════════════════════
# [PATTERN] Predicting C_EMP from a prior series of scales -- no full
# calibration
# ══════════════════════════════════════════════════════════════
# Problem: C_EMP requires p(n) at thousands of points. The per-call cost
# through primecount grows with n:
#     n=10^9 : 5.8ms   => 100K-point calibration = 10 minutes
#     n=10^12: 151ms   => 4 hours
#     n=10^13: 622ms   => 17 hours
# So the practical ceiling is ~10^12-10^13.
#
# Idea: calibrate on a series of **doublings** (5K, 10K, 20K, ...), infer
# the pattern, then predict the next scale instead of calibrating it.
#
# What measurement over 12 scales (10K to 20.48M) showed:
#   - C_EMP updates in **clusters**, not gradually: only 20 clusters over
#     20 million points (451 grouped updates), with millions of points
#     between them showing no change at all.
#   - coefficient of variation of the gaps = 1.46 (above pure Poisson 1.0)
#     => the clusters are **clustered**, not regular, so a jump is
#     possible at any time.
#   - but as the scale grows:
#       * clusters become rarer  (position-gap correlation r=+0.52)
#       * increments shrink      (position-increment correlation r=-0.45)
#     first 5 increments: 31.8% 2.5% 3.8% 37.1% 3.1%
#     last  5 increments:  1.7% 15.4% 0.2%  2.6% 3.5%
#
# Local prediction accuracy (from the last two scales, over 10
# predictions):
#     mean error +4.8%, worst undershoot -34.7%, worst overshoot +26.0%
#     and the last two predictions: -14.7% and -0.3%
#     => accuracy improves with scale
#
# Hence the margin **decreases with scale**: it covers the worst measured
# undershoot at small scales, and tightens at large ones where the
# prediction proved more accurate.
#
# **Warning:** this is an estimate, not a guarantee. C_EMP is a
# maximum-value statistic and the clusters are clustered, so a jump larger
# than any measured is possible at any scale.
GRID_TIME_BUDGET = 900.0    # [AUTO-GRID] time budget for building the grid (seconds)
GRID_MIN, GRID_MAX = 500, 20_000   # grid-size limits

# [MARGIN] The margin is derived statistically from the distribution of
# prediction error, not chosen by hand.
#
# **The model's actual range is a 2x jump only.** Training on [1, N]
# yields a valid containment on [N+1, 2N] and no further. So the
# distribution is measured at jumps <= 2x exclusively -- measuring it over
# longer jumps mixes in cases outside the operating range and spuriously
# inflates sigma (0.54 when mixing all jumps versus 0.147 when correctly
# restricted).
#
# Measured distribution of log(predicted / actual) at jumps <= 2x:
#     mu = +0.0172   sigma = 0.1473   (23 cases, scales 10K-20M)
#
# Margin = exp(z*sigma - mu):
#     95%   confidence -> 1.253
#     99%   confidence -> 1.385
#     99.7% confidence -> 1.529
#
# **Cross-region validation** (what distinguishes this from hand-tuning):
#   derived from 10K-20M and tested on an independent series 25M-62M that
#   never entered the derivation: **zero failures out of 9 cases** at all
#   three confidence levels.
#   The margin actually required there was 1.179, so 1.385 covered it with
#   17% to spare.
#
# The two former values (1.55 and 1.25) were hand-picked after seeing the
# results, with no derivation and no cross-region validation.
MARGIN_Z       = 3.00       # nominal confidence (3.0=99.7%, 2.33=99%, 1.65=95%)
# [CORRECTION] It was 2.33, giving a margin of 1.68. The only available
# paired measurement (at 2^40, where modes C and F were both run at the
# same scale) showed the prediction **undershot by 43%**: predicted
# 2,102.6 versus measured 3,681.8, so ~1.75 was required while the margin
# was 1.68 -- i.e. **4% short**.
# Containment still succeeded, but thanks to internal slack in the
# envelope (+/-2*C_EMP*ln^2), not to margin accuracy. Raising it to z=3.0
# makes it cover the available measurement.
#
# **Warning:** the confidence level is nominal, not real. The formula
# exp(z*sigma - mu) assumes a normal error distribution, which is
# **unverified** -- indeed C_EMP is a maximum-value statistic, and such
# distributions belong to a different, heavier-tailed family (Gumbel). The
# figure 99.7% therefore describes the formula, not reality.
# -- Paired measurements: prediction (F) versus direct measurement at the
#    **same scale** --
# Each row: (scale, measured C_EMP, raw predicted C_EMP)
# The measurement comes strictly from inside [1, N] -- by full sieve or by
# a primecount sample.
# The prediction uses only scales **smaller than N** (no data at N itself).
# The ratio (measured / predicted) is the margin actually required at that
# scale.
#
# This is the sole source of MARGIN_FLOOR: not a chosen number, but the
# largest ratio observed. Add any new paired measurement here and the
# margin updates automatically.
MARGIN_CALIBRATION = [
    (2_000_000,          8.009,     10.309),   # ratio 0.777 (prediction overshot)
    (4_000_000,          9.696,     12.600),   # 0.770 (overshot)
    (8_000_000,         14.313,     12.001),   # 1.193
    (12_000_000,        15.237,     15.461),   # 0.986
    (20_000_000,        19.224,     18.705),   # 1.028
    (2**40,           3681.800,   2102.565),   # 1.751  <- worst observed
]

def _margin_floor_from_data(safety=1.00):
    """
    The largest margin observed to be actually required, from
    MARGIN_CALIBRATION.

    safety=1.00 means no manual inflation whatsoever -- the value is the
    worst measured, nothing more. (An earlier version used 1.85 = 1.751
    plus a 6% guess; that increment had no source.)

    The six current measurements span ratios 0.770 - 1.751, mean 1.084.
    In four of them the prediction already overshot (ratio < 1), so no
    margin was needed there.

    **Caveat:** five of the six lie at scales 2M-20M, and only one at
    2^40 -- which is the worst. The possibility that the required margin
    grows with scale remains open, and one point cannot settle it.
    """
    if not MARGIN_CALIBRATION: return 1.20
    return max(m/p for _, m, p in MARGIN_CALIBRATION if p > 0) * safety

MARGIN_FLOOR   = _margin_floor_from_data()
MARGIN_MAX_JMP = 2.5        # largest jump admitted (the model's operating range)
MARGIN_MIN_N   = 8          # minimum number of cases to accept a measurement

def measure_margin_params(scales, values):
    """
    Computes (mu, sigma, n) of the prediction-error distribution **from
    the given series** -- not from constants.

    Method: for every pair within the series whose jump is
    <= MARGIN_MAX_JMP, predict the second point from the first using the
    same algorithm as predict_c_emp, and record log(predicted / actual).
    mu and sigma are the statistics of that distribution.

    Why measured on every call rather than frozen? Because the series
    grows: each new direct calibration added changes the distribution.
    Freezing (mu, sigma) would make the rule operate on stale inputs and
    ignore every later measurement -- losing its purpose and becoming a
    fixed number dressed as a formula.

    Restricting to jumps <= 2.5x is deliberate: the model's operating
    range is [N+1, 2N], and longer jumps are not use cases. Including them
    spuriously inflates sigma (0.54 when mixed versus 0.147 when properly
    restricted).
    Returns None if there are too few cases, in which case the margin
    falls back to MARGIN_FLOOR.
    """
    s = np.asarray(scales, dtype=float); v = np.asarray(values, dtype=float)
    ok = v > 0; s, v = s[ok], v[ok]
    if len(s) < 3: return None
    rs = []
    for i in range(2, len(s)):
        if s[i-1]/s[i-2] > MARGIN_MAX_JMP: continue
        k = float(np.clip(math.log(v[i-1]/v[i-2]) / math.log(s[i-1]/s[i-2]), 0.0, 1.0))
        for j in range(i, len(s)):
            if s[j]/s[i-1] > MARGIN_MAX_JMP: continue
            raw = v[i-1] * (s[j]/s[i-1])**k
            if raw > 0 and v[j] > 0: rs.append(math.log(raw/v[j]))
    if len(rs) < MARGIN_MIN_N: return None
    a = np.asarray(rs)
    return float(a.mean()), float(a.std()), len(rs)

PATTERN_MARGIN_SMALL = 1.55   # (legacy -- unused, kept for compatibility)
PATTERN_MARGIN_LARGE = 1.25
PATTERN_LARGE_SCALE  = 10_000_000

def _self_margin(scales, values):
    """
    [SELF-MARGIN] Measures the margin from the accuracy of earlier
    predictions **within the same series**: predict each scale from the
    one before it, record the ratio (predicted / actual), then take the
    worst undershoot and compensate for it.

    Why not used alone? Because it learns only from the past, so the first
    large jump catches it out. Measured on the actual series: coverage 8/9
    -- it failed at TRAIN_END=1,280,000 where the error was -34.7% while
    the learned margin was only 1.10 (because every prior error had been
    positive).
    It is therefore combined with a fixed floor: max(self, floor).
    """
    s = np.asarray(scales, dtype=float); v = np.asarray(values, dtype=float)
    ok = v > 0; s, v = s[ok], v[ok]
    if len(s) < 3: return None
    worst = 1.0; used = 0
    for i in range(2, len(s)):
        # [SCOPE-FIX] Restrict the measurement to jumps <= 2.5x -- the
        # model's actual operating range (training on [1,N] contains
        # [N+1,2N]).
        # Longer jumps within the series (e.g. 20,480,000 -> 5x10^11 =
        # 24,414x) are not operating cases, and including them spuriously
        # inflated the margin (1.93 instead of 1.39).
        if s[i]/s[i-1] > 2.5: continue
        k = math.log(v[i-1]/v[i-2]) / math.log(s[i-1]/s[i-2])
        raw = v[i-1] * (s[i]/s[i-1])**float(np.clip(k, 0.0, 1.0))
        worst = min(worst, raw/v[i]); used += 1
    if used < 3 or worst <= 0: return None
    return (1.0/worst) * 1.10

def predict_c_emp(scales, values, TRAIN_END):
    """
    Predicts C_EMP at TRAIN_END from a prior series of (scale, value).
    Method: a local exponent from the last two scales (more accurate than
    a global fit: mean error 16.1% versus 22.9%), multiplied by a margin.

    Margin = max(self-measured from the series, fixed floor).
    The floor guards against the first jump (where the self-measured value
    alone fails), while the self-measured value widens automatically if a
    jump larger than any observed appears. As the series lengthens and
    accuracy improves, the margin approaches its floor -- i.e. it
    **tightens automatically**.
    """
    s = np.asarray(scales, dtype=float); v = np.asarray(values, dtype=float)
    ok = v > 0; s, v = s[ok], v[ok]
    if len(s) < 2:
        raise ValueError("at least two actual scales are required")
    k = math.log(v[-1]/v[-2]) / math.log(s[-1]/s[-2])
    k = float(np.clip(k, 0.0, 1.0))
    raw = v[-1] * (float(TRAIN_END)/s[-1])**k
    # [MARGIN] The statistically derived margin, widened automatically if
    # the self-measured margin reveals an error larger than any seen
    # during derivation.
    # [MARGIN] mu and sigma are measured from the series itself on every
    # call, so any newly added direct calibration updates them
    # automatically.
    _mp = measure_margin_params(scales, values)
    derived = math.exp(MARGIN_Z * _mp[1] - _mp[0]) if _mp else MARGIN_FLOOR
    self_m = _self_margin(s, v)
    margin = max(MARGIN_FLOOR, derived, self_m or 0.0)
    return raw * margin, raw, k, margin

def calibration_from_pattern(TRAIN_END, P_ALL, scales, pn_values, pi_values,
                             verbose=True, grid=None, measure_c_emp=False):
    """
    Full calibration **without any sieving**: coefficients from a grid via
    primecount, and C_EMP by pattern prediction.

    [FIX] The first version called compute_calibration, which needs
    P_ALL -- so it failed exactly at the scales that cannot be sieved,
    which is the case mode F exists for. Now: P_ALL is used if available,
    otherwise the grid values are fetched from primecount (a few hundred
    thousand calls at most).
    """
    # [PERF] Grid size when P_ALL is absent: each point is a primecount
    # call, whose cost grows with TRAIN_END. At 2^50 a call takes about a
    # second, so 200K points would take days.
    # 20K points suffice for the coefficients (mean-type statistics --
    # verified that a grid is enough for them).
    if grid is None:
        if P_ALL is not None:
            grid = min(CALIB_GRID, 200_000)
        else:
            # [AUTO-GRID] Each point is a primecount call, whose cost
            # explodes with n:
            #   10¹² -> 222ms | 10¹³ -> 633ms | 10¹⁴ -> 2.1s | 2⁵⁰ -> 9.9s
            # A fixed size (20,000) meant 55 hours at 2^50 with no
            # warning at all.
            # Now: measure the cost of one call directly, and compute the
            # size that fits GRID_TIME_BUDGET -- so this cannot recur at
            # any scale.
            #
            # The floor of 500 is justified by measurement: the effect of
            # grid size on estimation error (TRAIN_END=5M) was 878 at 500
            # versus 1,040 at 200,000 -- i.e. the small grid is not worse
            # (the coefficients are mean-type statistics).
            _t0 = time.time()
            _m = _pc()
            if _m is not None:
                _m.nth_prime(int(TRAIN_END))
            _per = max(time.time() - _t0, 1e-6)
            grid = int(np.clip(GRID_TIME_BUDGET / _per, GRID_MIN, GRID_MAX))
            if verbose:
                _est = grid * _per
                print(f"  [auto-grid] one call={_per*1000:.0f}ms => "
                      f"grid {grid:,} points (~{_est/60:.0f} min)")
        grid = min(grid, TRAIN_END)
    idx = np.unique(np.linspace(0, TRAIN_END-1, min(grid, TRAIN_END)).astype(np.int64))
    n_g = (idx + 1).astype(np.float64)
    if P_ALL is not None:
        p_g = P_ALL[idx].astype(np.float64)
    else:
        m = _pc()
        if m is None:
            raise RuntimeError(
                "\n" + "="*62 + "\n"
                "primecount is not installed -- required for scales beyond the sieve limit.\n\n"
                "Install (one cell before running):\n"
                "    !pip install primecountpy\n\n"
                "then re-run this cell.\n" + "="*62)
        if verbose:
            print(f"  [pattern] fetching {len(idx):,} grid points via primecount...")
        p_g = np.array([m.nth_prime(int(v)) for v in n_g], dtype=np.float64)

    li_inv_g = inv_R_vec(n_g)
    X = build_pn_features(li_inv_g)
    y = p_g - li_inv_g
    si = np.random.RandomState(42).choice(len(idx), min(15000, len(idx)), replace=False)
    al = RidgeCV(alphas=np.logspace(-1,3,30), fit_intercept=False, cv=3).fit(X[si], y[si]).alpha_
    Rp = Ridge(alpha=al, fit_intercept=False, solver='cholesky').fit(X, y)
    median_pn = float(np.median(y - Rp.predict(X)))
    del X; gc.collect()

    li_g = R_vec(p_g)
    r0b = n_g - li_g
    sm = savgol_filter(r0b, min(51, len(r0b)-(1-len(r0b)%2)), 3)
    mk = (np.abs(r0b-sm) <= np.percentile(np.abs(r0b-sm), 70.4)) & (n_g >= 10)
    def _bf(x,a,b): return a*np.sqrt(x/np.log(np.maximum(x,3.)))+b
    (A_B, B_B), _ = curve_fit(_bf, p_g[mk], r0b[mk], maxfev=10000)
    Fb_pi_g = li_g + _bf(p_g, A_B, B_B)
    Xp = build_pi_features(p_g)
    pi_residual = n_g - Fb_pi_g
    alpha_pi = RidgeCV(alphas=np.logspace(-1,3,30), fit_intercept=False, cv=3).fit(Xp[si], pi_residual[si]).alpha_
    Ri = Ridge(alpha=alpha_pi, fit_intercept=False, solver='cholesky').fit(Xp, pi_residual)
    median_pi = float(np.median((n_g - Fb_pi_g) - Ri.predict(Xp)))
    del Xp; gc.collect()

    # base slope from the grid (substitute for the convex hull when P_ALL
    # is unavailable)
    if P_ALL is not None:
        TSN = lower_hull_tsn_base(TRAIN_END, P_ALL)
    else:
        TSN = float(np.min(np.diff(p_g)/np.diff(n_g)))

    if measure_c_emp:
        # [C] Measure C_EMP by primecount sampling instead of predicting
        # -- more accurate but slower
        pn_pred = c_emp_via_primecount(TRAIN_END, Rp.coef_, median_pn, 'pn')
        pi_pred = c_emp_via_primecount(TRAIN_END, Ri.coef_, median_pi, 'pi',
                                       A_BIAS=A_B, B_BIAS=B_B)
        pn_raw, k_pn, m_ = pn_pred/PC_MARGIN, float('nan'), PC_MARGIN
    else:
        pn_pred, pn_raw, k_pn, m_ = predict_c_emp(scales, pn_values, TRAIN_END)
        pi_pred, pi_raw, k_pi, _  = predict_c_emp(scales, pi_values, TRAIN_END)

    def c_vec(a_): return np.zeros_like(np.asarray(a_, dtype=np.float64))
    cal_pn = dict(TRAIN_END=TRAIN_END, c_vec=c_vec, dc_dn=lambda n: 0.0,
                  ridge_coef_pn=Rp.coef_, DECAY=0.999, C_EMP_PN=pn_pred,
                  TSN_BASE=TSN, THRESH_PN=2.28*(TRAIN_END**0.486)*0.84,
                  median_corr_pn=median_pn, RHO=None, ANOM=None,
                  C_EMP_PREDICTED=True)
    cal_pi = dict(TRAIN_END=TRAIN_END, A_BIAS_PI=A_B, B_BIAS_PI=B_B,
                  ridge_coef=Ri.coef_, C_EMP_PI=pi_pred, DECAY_PI=0.999,
                  CP_PI=5, THRESH_PI=0.3, x_min=float(p_g[0]),
                  median_corr_pi=median_pi, C_EMP_PREDICTED=True)
    if verbose:
        _lbl = "measured" if measure_c_emp else "pattern"
        print(f"  [{_lbl} TRAIN_END={TRAIN_END:,}] C_EMP_PN={pn_pred:.3f} "
              f"(raw {pn_raw:.3f}, k={k_pn:.3f}) | C_EMP_PI={pi_pred:.4f} | margin={m_:.2f}")
    return cal_pn, cal_pi


STREAM_CEMP = False   # True = full C_EMP scan without storage (coverage + bounded memory)


USE_PRIMECOUNT = False      # True = obtain p(n) via primecount instead of sieving
PC_SAMPLE      = 100_000    # number of sample points for C_EMP
PC_MARGIN      = 1.10       # margin compensating the sampling shortfall -- justified statistically:
# Measured shortfall (C_EMP from a 100K sample / full C_EMP) across 5
# scales from 2M to 20M:
#     0.9780  0.9829  0.9532  0.9835  0.9917
#     min=0.9532   mean=0.9779   sd=0.0131
# Required margins:
#     to cover the worst measured : 1/0.9532 = 1.049
#     to cover (mean - 3 sigma)   : 1/0.9386 = 1.065
# We chose 1.10 -- wider than both with a safety gap, costing only 10% in
# width.
#
# **Warning: this is not a guarantee.** A method for computing the margin
# automatically from the sample itself was tried (a saturation indicator,
# C_EMP(50K)/C_EMP(100K)) and failed: its correlation with the actual
# shortfall was r=-0.32, with inverted sign. The shortfall is therefore
# **estimated statistically from prior measurements, not computed from
# the current run**.
# And C_EMP is a maximum-value statistic in any case, so it carries no
# guarantee -- the only proven component in the model is the Axler bound.

_PC = None
def _pc():
    """Loads primecount on first use; returns None if unavailable."""
    global _PC
    if _PC is None:
        try:
            from primecountpy import primecount as _m
            _PC = _m
        except ImportError:
            try:
                import primecount as _m
                _PC = _m
            except ImportError:
                _PC = False
    return _PC if _PC else None

def nth_prime_vec(ns):
    """p(n) for an array of indices, via primecount."""
    m = _pc()
    if m is None:
        raise RuntimeError("primecount unavailable -- install it: pip install primecountpy")
    return np.array([m.nth_prime(int(v)) for v in np.atleast_1d(ns)], dtype=np.float64)

def c_emp_via_primecount(TRAIN_END, ridge_coef, median_corr, side='pn',
                         A_BIAS=None, B_BIAS=None, K=None):
    """
    Computes C_EMP over a spread sample, fetching p(n) from primecount --
    no sieve. The result is multiplied by PC_MARGIN to compensate for what
    the sample misses.
    """
    K = K or PC_SAMPLE
    ns = np.unique(np.linspace(C_EMP_MIN_ARG, TRAIN_END, int(min(K, TRAIN_END))).astype(np.int64))
    nf = ns.astype(np.float64)
    pv = nth_prime_vec(ns)
    if side == 'pn':
        bx = inv_R_vec(nf)
        F = bx + fast_ridge_correction(bx, G60_PN, ridge_coef) + median_corr
        ratio = np.abs(pv - F) / _NORM(nf)
    else:
        F = (R_vec(pv) + (A_BIAS*np.sqrt(pv/np.log(np.maximum(pv,3.)))+B_BIAS)
             + fast_ridge_correction(pv, G60_PI, ridge_coef) + median_corr)
        ratio = np.abs(nf - F) / _NORM(pv)
    return float(np.max(ratio)) * PC_MARGIN


# depending on the seed (measured: 3 of 4 targets give different bounds,
# differing by 1-3%). Taking the union over several spread seeds makes the
# result **independent of any particular choice**.

def inv_R_vec_fast(ns):
    """[PORT-5] Same as inv_R_vec but with the fast li seed -- used
    **inside the C_EMP loop only** (the one place that traverses millions
    of points). The walk and jump paths keep the original inv_R_vec
    unchanged."""
    ns = np.asarray(ns, dtype=float)
    x = inv_li_vec_fast_bulk(ns)
    _check_D_range(x, "inv_R_vec_fast")            # [GUARD]
    x = x - D_corr(x)*np.log(x)
    lx = np.log(x)
    return x - ((li_vec(x) + D_corr(x)) - ns)*lx

ANOMALY_OFF = False         # [TEST-SWITCH] True = disable the anomaly mechanism (testing only)

# [EXTRACT] How the interval is extracted from the candidate cloud.
# 0.0 = min/max over all candidates (original behaviour, **default**)
# >0  = percentile: trim EXTRACT_PCTL% from each tail before taking the
#       bounds.
#
# Why might it help? With min/max the result is set by the two extreme
# candidates alone, so any improvement in the cloud's **density** (from
# the anomaly mechanism or from branching) evaporates and never shows.
# With a percentile, density is what determines the result.
#
# Preliminary measurement (18 targets, TRAIN_END=2M):
#     min/max        : width 7,365  contained 18/18  min margin 21.7%
#     percentile 0.5%: width 7,032  contained 18/18  min margin 22.0%
#     percentile 2%  : width 6,578  contained 18/18  min margin 18.8%
#     percentile 5%  : width 6,003  contained 18/18  min margin 16.0%
#
# **Warning:** the percentile abandons coverage in principle (it trims
# tails deliberately), and 18 targets is a very small sample. Do not adopt
# it before measuring on hundreds of targets across at least 3 scales.
# The default 0.0 leaves behaviour unchanged.
EXTRACT_PCTL = 0.0

SAMPLE_CAP = 500_000
CALIB_GRID = 600_000
CP_PN = 40


# ══════════════════════════════════════════════════════════════
# [PORT-3] Parallel execution of the two C_EMP loops
# ══════════════════════════════════════════════════════════════
# C_EMP_PN and C_EMP_PI are entirely independent (each needs only P_ALL
# and its own coefficients), and they are the two most expensive items in
# calibration. Running them together in two processes saves about half.
# P_ALL is inherited by fork without copying (Linux/Colab). Under spawn we
# fall back to sequential execution.
# Note: **numbers** are passed, not closures -- closures are not
# picklable.
_MP_P = None

def _mp_cemp(args):
    side, TE, prm = args
    P_ALL = _MP_P
    best = 0.0
    for cs in range(0, TE, C_EMP_CHUNK):
        ce = min(cs + C_EMP_CHUNK, TE)
        if side == 'pn':
            arg = np.arange(cs+1, ce+1, dtype=np.float64)
            b = inv_R_vec_fast(arg)
            F = b + fast_ridge_correction(b, G60_PN, prm['rc']) + prm['med']
            tgt = P_ALL[cs:ce]; den = arg
        else:
            arg = P_ALL[cs:ce].astype(np.float64)
            F = (R_vec(arg) + (prm['A']*np.sqrt(arg/np.log(np.maximum(arg,3.)))+prm['B'])
                 + fast_ridge_correction(arg, G60_PI, prm['rc']) + prm['med'])
            tgt = np.arange(cs+1, ce+1, dtype=np.float64); den = arg
        k = den >= 10_000
        if k.any():
            m = float(np.max(np.abs(tgt[k] - F[k]) / _NORM(den[k])))
            if m > best: best = m
        del arg, F, tgt, den, k
        if cs % (C_EMP_CHUNK*10) == 0: gc.collect()
    return best

def _cemp_pair(TRAIN_END, P_ALL, pn_prm, pi_prm, parallel):
    global _MP_P
    _ensure_D()   # [FIX] built once in the parent before fork, so both
                  # children inherit it without rebuilding
    if parallel:
        import multiprocessing as _mp
        if _mp.get_start_method(allow_none=True) in ('fork', None):
            try:
                import concurrent.futures as cf
                print("  [C_EMP: running in parallel on two processes]")
                _MP_P = P_ALL
                with cf.ProcessPoolExecutor(max_workers=2) as ex:
                    f1 = ex.submit(_mp_cemp, ('pn', TRAIN_END, pn_prm))
                    f2 = ex.submit(_mp_cemp, ('pi', TRAIN_END, pi_prm))
                    r = (f1.result(), f2.result())
                _MP_P = None
                return r
            except Exception as _e:
                print(f"  [parallel failed: {_e} -- falling back to sequential]")
                _MP_P = None
    if parallel: print("  [parallel unavailable (spawn) -- sequential]")
    _MP_P = P_ALL
    r = (_mp_cemp(('pn', TRAIN_END, pn_prm)), _mp_cemp(('pi', TRAIN_END, pi_prm)))
    _MP_P = None
    return r



# ══════════════════════════════════════════════════════════════
# [CACHE] Saving and reusing a calibration
# ══════════════════════════════════════════════════════════════
# Calibration costs 2,388s at TRAIN_END=500M, whereas a single jump takes
# 11s. Saving it once and reloading makes any later target from the same
# TRAIN_END almost instantaneous.
# Note: c_vec and dc_dn are functions (not picklable), so we store their
# parameters and rebuild them on load -- the result is **bit-identical**.
def save_calibration(path, cal_pn, cal_pi):
    d_pn = {k: v for k, v in cal_pn.items() if k not in ('c_vec', 'dc_dn')}
    with open(path, 'wb') as f:
        pickle.dump({'pn': d_pn, 'pi': cal_pi, 'ver': 'original-architecture-math-only-v14-12mobius'}, f)

def load_calibration(path):
    with open(path, 'rb') as f:
        d = pickle.load(f)
    if d.get('ver') != 'original-architecture-math-only-v14-12mobius':
        raise ValueError('Old calibration cache version is incompatible with corrected equations')
    cal_pn = d['pn']
    # c(n) is removed in this version => both functions are zero (see [C-REMOVED])
    cal_pn['c_vec'] = lambda n_arr: np.zeros_like(np.asarray(n_arr, dtype=np.float64))
    cal_pn['dc_dn'] = lambda n: 0.0
    return cal_pn, d['pi']


def compute_calibration(TRAIN_END, P_ALL, verbose=True, parallel=False):
    t_start = time.time()
    rng = np.random.RandomState(42)
    grid_n = min(CALIB_GRID, TRAIN_END)
    idx_grid = np.unique(np.linspace(0, TRAIN_END-1, grid_n).astype(np.int64))
    n_g = (idx_grid + 1).astype(np.float64)
    p_g = P_ALL[idx_grid].astype(np.float64)
    li_inv_g = inv_R_vec(n_g)          # [UPGRADE-R]

    r0 = p_g - li_inv_g
    r0_sm = savgol_filter(r0, min(51, len(r0)-(1-len(r0)%2)), 3)
    mask = (np.abs(r0-r0_sm) <= np.percentile(np.abs(r0-r0_sm), 70.4)) & (n_g>=10)
    # [C-REMOVED] c(n)=a*sqrt(n*ln n)+b corrected the bias of li^-1. After
    # upgrading the base to R^-1 that bias vanished, so it began fitting
    # noise.
    # Measured on 7,997 out-of-training targets (with full Ridge refit):
    #     with c(n)   : mean error 599   median 532
    #     without c(n): mean error 554   median 480   (8.15% better, p=1e-47)
    # Its purpose (centring the model on the sample) is served more
    # accurately by median_corr -- which is computed over **all** points
    # after Ridge, not on a partial 70.4% mask.
    a_c, b_c = 0.0, 0.0
    def c_vec(n_arr): return np.zeros_like(np.asarray(n_arr, dtype=np.float64))
    def dc_dn(n): return 0.0
    del r0, r0_sm, mask; gc.collect()

    samp_idx = rng.choice(len(idx_grid), min(SAMPLE_CAP, len(idx_grid)), replace=False)
    Fb_pn_samp = li_inv_g[samp_idx] + c_vec(n_g[samp_idx])
    y_pn = p_g[samp_idx] - Fb_pn_samp
    Xs_pn = build_pn_features(li_inv_g[samp_idx])
    si2 = rng.choice(len(samp_idx), min(15000,len(samp_idx)), replace=False)
    alpha_pn = RidgeCV(alphas=np.logspace(-1,3,30), fit_intercept=False, cv=3).fit(Xs_pn[si2], y_pn[si2]).alpha_
    Rdg_pn = Ridge(alpha=alpha_pn, fit_intercept=False, solver='cholesky').fit(Xs_pn, y_pn)
    ridge_coef_pn = Rdg_pn.coef_
    del Xs_pn, y_pn; gc.collect()

    Fb_pn_grid = li_inv_g + c_vec(n_g)
    Xg_pn = build_pn_features(li_inv_g)
    base_pred_pn = Rdg_pn.predict(Xg_pn)
    median_pn = float(np.median((p_g - Fb_pn_grid) - base_pred_pn))
    del Fb_pn_grid, Xg_pn, base_pred_pn
    # [FIX-2] It was n_g >= 10. Measured at TRAIN_END=2,000,000:
    #   cutoff n>=10     -> C_EMP=64.250, worst point at n=11 (!)
    #   cutoff n>=100    -> C_EMP=16.800, worst point at n=111
    #   cutoff n>=10,000 -> C_EMP= 7.553, worst point at n=1,887,313
    #                       (a genuine gap)
    # Reason: ln^2(11)=5.75 is a trivially small denominator, so dividing
    # by it inflates the ratio spuriously.
    # The value is entirely stable from n>=10,000 upward. Impact: C_EMP
    # inflated 8.5x.
    # [NORM] Normalisation: ln^2(n) is borrowed from Cramer's conjecture
    # on the largest prime gap, but C_EMP measures the **model's** residual,
    # not gaps. Measured across 4 bins (dispersion = max/min; ideal 1.0):
    #     ln^2(n) [current] : 1.57x
    #     √n·ln(n)        : 1.21x
    #     n^0.4           : 1.11x
    #     √n              : 1.10x
    #     sqrt(n)/ln(n)     : 1.08x   <- most uniform
    # sqrt(n)/ln(n) is also the form the Riemann hypothesis suggests for
    # the size of pi(x) fluctuations, so it is not a purely empirical
    # choice. It was nevertheless rejected: see the [NORM] note above --
    # it widens the interval 1.8x under extrapolation.
    # [PORT-1] C_EMP over the **entire** training range (every point, in
    # chunks) instead of a fixed 600K grid. C_EMP is a **maximum-value**
    # statistic, so a sample can miss the true extreme. Chunking is
    # necessary for memory (a single batch at 50M would OOM).
    del Fb_pn_samp; gc.collect()
    if STREAM_CEMP:
        # [STREAM] full scan without storage -- full coverage, bounded
        # memory, no margin needed
        C_EMP_PN = c_emp_streaming(TRAIN_END, ridge_coef_pn, median_pn, 'pn')
    elif USE_PRIMECOUNT and _pc() is not None:
        # [PRIMECOUNT] sampling without a sieve -- higher ceiling, but
        # with an estimated margin
        C_EMP_PN = c_emp_via_primecount(TRAIN_END, ridge_coef_pn, median_pn, 'pn')
    _PN_PRM = dict(rc=ridge_coef_pn, med=median_pn)   # [PORT-3] computed later, in parallel with PI
    if not (STREAM_CEMP or (USE_PRIMECOUNT and _pc() is not None)):
        C_EMP_PN = None

    # [FIX-3] It used to take a 3000-point sample of the checkpoints and
    # then compute the correlation between them. But after subsampling,
    # the spacing between adjacent points became (range/3000), not CP_PN.
    # At TRAIN_END=50M this measures a lag of ~16,680 instead of 40 -- an
    # entirely wrong lag.
    # Proved with a synthetic AR(1) series of known correlation 0.85:
    #     old method    -> -0.0355  (error 0.886)
    #     block method  -> +0.849   (error 0.001)
    # Fix: short blocks of **consecutive** points at the true CP spacing,
    # distributed over the range, preserving lag=CP while still covering
    # the whole range.
    _nb, _bl = 300, min(40, max(3, TRAIN_END // (CP_PN * 6)))
    _span_end = max(CP_PN + 1.0, TRAIN_END - CP_PN*_bl - 1)
    _anchors = np.linspace(CP_PN, _span_end, min(_nb, max(2, int(_span_end/CP_PN))))
    _pts = np.empty(len(_anchors)*_bl, dtype=np.float64); _i = 0
    for _a in _anchors:
        _c = float(_a)
        for _ in range(_bl):
            _c += CP_PN; _pts[_i] = _c; _i += 1
    # MATH-FIX: both prime truth and full fitted model at SAME integer rank.
    _rank = np.rint(_pts).astype(np.int64)
    if np.any(_rank < 1) or np.any(_rank > TRAIN_END):
        raise ArithmeticError('DECAY_PN training checkpoints outside range')
    _true = P_ALL[_rank-1].astype(np.float64)
    _rank_float = _rank.astype(np.float64)
    _base = inv_R_vec(_rank_float)
    _model = (_base + c_vec(_rank_float)
              + fast_ridge_correction(_base, G60_PN, ridge_coef_pn) + median_pn)
    _rb = (_true - _model).reshape(len(_anchors), _bl)
    _a1 = _rb[:, :-1].ravel(); _a2 = _rb[:, 1:].ravel()
    DECAY_pn, _pn_corr_defined = _finite_correlation(_a1, _a2)
    del _pts, _rank, _rank_float, _base, _true, _model, _rb, _a1, _a2; gc.collect()

    # [ONLINE] online-mode measurements -- cheap (a few thousand points)
    _cl_on = dict(c_vec=c_vec, ridge_coef_pn=ridge_coef_pn, median_corr_pn=median_pn)
    try:
        RHO  = measure_residual_rho(TRAIN_END, P_ALL, _cl_on)
        ANOM = measure_anomaly_persistence(TRAIN_END, P_ALL, _cl_on)
    except Exception as _e:
        print(f'  [ONLINE] RHO/ANOM measurement failed: {_e}')
        RHO, ANOM = None, None
    TSN_BASE = lower_hull_tsn_base(TRAIN_END, P_ALL)
    THRESH_PN = 2.28*(TRAIN_END**0.486)*0.84
    cal_pn = dict(TRAIN_END=TRAIN_END, c_vec=c_vec, dc_dn=dc_dn, ridge_coef_pn=ridge_coef_pn,
                  DECAY=DECAY_pn, C_EMP_PN=C_EMP_PN, TSN_BASE=TSN_BASE, THRESH_PN=THRESH_PN,
                  median_corr_pn=median_pn, RHO=RHO, ANOM=ANOM, DECAY_PN_CORRELATION_DEFINED=_pn_corr_defined)

    li_g = R_vec(p_g)                 # [UPGRADE-R]
    r0b = n_g - li_g
    r0b_sm = savgol_filter(r0b, min(51, len(r0b)-(1-len(r0b)%2)), 3)
    maskb = (np.abs(r0b-r0b_sm) <= np.percentile(np.abs(r0b-r0b_sm), 70.4)) & (n_g>=10)
    def _bf(x,a,b): return a*np.sqrt(x/np.log(np.maximum(x,3.)))+b
    (A_BIAS_PI,B_BIAS_PI),_ = curve_fit(_bf, p_g[maskb], r0b[maskb], maxfev=10000)
    Fb_pi_g = li_g + (A_BIAS_PI*np.sqrt(p_g/np.log(np.maximum(p_g,3.)))+B_BIAS_PI)
    del r0b, r0b_sm, maskb; gc.collect()

    samp_idx2 = rng.choice(len(idx_grid), min(SAMPLE_CAP, len(idx_grid)), replace=False)
    Xs_pi = build_pi_features(p_g[samp_idx2])
    y_pi = (n_g - Fb_pi_g)[samp_idx2]
    si3 = rng.choice(len(samp_idx2), min(15000,len(samp_idx2)), replace=False)
    alpha_pi = RidgeCV(alphas=np.logspace(-1,3,30), fit_intercept=False, cv=3).fit(Xs_pi[si3], y_pi[si3]).alpha_
    Rdg_pi = Ridge(alpha=alpha_pi, fit_intercept=False, solver='cholesky').fit(Xs_pi, y_pi)
    del Xs_pi, y_pi; gc.collect()

    Xg_pi = build_pi_features(p_g)
    base_pred_pi = Rdg_pi.predict(Xg_pi)
    median_pi = float(np.median((n_g - Fb_pi_g) - base_pred_pi))
    del Xg_pi, base_pred_pi, Fb_pi_g, li_g, idx_grid, n_g, p_g, li_inv_g; gc.collect()
    # [PORT-1] same principle for the pi(x) side
    _PI_PRM = dict(rc=Rdg_pi.coef_, med=median_pi, A=A_BIAS_PI, B=B_BIAS_PI)
    if STREAM_CEMP:
        C_EMP_PI = c_emp_streaming(TRAIN_END, Rdg_pi.coef_, median_pi, 'pi',
                                   A_BIAS=A_BIAS_PI, B_BIAS=B_BIAS_PI)
    elif USE_PRIMECOUNT and _pc() is not None:
        C_EMP_PI = c_emp_via_primecount(TRAIN_END, Rdg_pi.coef_, median_pi, 'pi',
                                        A_BIAS=A_BIAS_PI, B_BIAS=B_BIAS_PI)
    else:
        C_EMP_PN, C_EMP_PI = _cemp_pair(TRAIN_END, P_ALL, _PN_PRM, _PI_PRM, parallel)
    cal_pn['C_EMP_PN'] = C_EMP_PN
    gc.collect()

    CP_PI = 5
    # [FIX-CONSISTENCY] It used to walk **every** checkpoint from the
    # first prime to the end of training -- 10 million points at
    # TRAIN_END=50M, with a full build_pi_features for them. It was also
    # **inconsistent** with DECAY_PN, which had been fixed to the block
    # method (since a spread-out sample measures the wrong lag).
    # Both sides now use the same method: short blocks of consecutive
    # points at the true CP spacing, distributed over the range.
    # Measured at TRAIN_END=5M:
    #     full walk : DECAY_PI=0.99846695  in 3.67s (1,000,138 points)
    #     blocks    : DECAY_PI=0.99849368  in 0.04s (12,000 points)
    #     difference 2.7e-05, 84x faster
    _nbp, _blp = 300, 40
    _xtrain_max = float(P_ALL[TRAIN_END-1])
    _max_anchor = _xtrain_max - (_blp + 1) * CP_PI * math.log(_xtrain_max)
    if _max_anchor <= 2.5:
        raise ValueError('insufficient training domain for DECAY_PI checkpoints')
    _anch_pi = np.minimum(inv_R_vec(np.linspace(CP_PN, TRAIN_END-1, _nbp)), _max_anchor)
    _pts_pi = np.empty(_nbp*_blp, dtype=np.float64); _j = 0
    for _aa in _anch_pi:
        _c = float(_aa)
        for _ in range(_blp):
            _c += CP_PI*math.log(max(_c, 3.)); _pts_pi[_j] = _c; _j += 1
    if np.any(_pts_pi > _xtrain_max) or np.any(_pts_pi < 2.0):
        raise ArithmeticError('DECAY_PI checkpoints escaped training; no count clipping permitted')
    # [COMPACT-FIX] It used searchsorted over P_ALL[:TRAIN_END], which
    # materialises the whole array (8 GB at 10^9). Now each point is
    # searched only within the anchor block it falls in, so only small
    # slices are reconstructed.
    if isinstance(P_ALL, CompactPrimes):
        _tk = np.empty(len(_pts_pi), dtype=np.float64)
        _ord = np.argsort(_pts_pi); _sp = _pts_pi[_ord]
        _ai = np.searchsorted(P_ALL.anchors, _sp, side='right') - 1
        _ai = np.clip(_ai, 0, len(P_ALL.anchors)-1)
        _i = 0
        while _i < len(_sp):
            _blk = _ai[_i]; _j = _i
            while _j < len(_sp) and _ai[_j] == _blk: _j += 1
            _lo = 1 + _blk*ANCHOR_STRIDE
            _hi = min(_lo + ANCHOR_STRIDE + 1, TRAIN_END)
            _seg = P_ALL[_lo:_hi]
            # [EDGE-FIX] a point beyond the last prime in range used to
            # return TRAIN_END+1; we clip it at TRAIN_END to match the
            # original searchsorted.
            _tk[_ord[_i:_j]] = _lo + np.searchsorted(_seg, _sp[_i:_j], side='right')
            _i = _j; del _seg
    else:
        _tk = np.searchsorted(P_ALL[:TRAIN_END], _pts_pi, side='right').astype(np.float64)
    if np.any(_tk > TRAIN_END):
        raise ArithmeticError('DECAY_PI true count exceeds training')
    _fz = ((A_BIAS_PI*np.sqrt(_pts_pi/np.log(np.maximum(_pts_pi,3.)))+B_BIAS_PI)
           + R_vec(_pts_pi)
           + fast_ridge_correction(_pts_pi, G60_PI, Rdg_pi.coef_) + median_pi)
    _rbp = (_tk - _fz).reshape(_nbp, _blp)
    DECAY_PI, _pi_corr_defined = _finite_correlation(_rbp[:,:-1].ravel(), _rbp[:,1:].ravel())
    del _pts_pi, _tk, _fz, _rbp; gc.collect()

    cal_pi = dict(TRAIN_END=TRAIN_END, A_BIAS_PI=A_BIAS_PI, B_BIAS_PI=B_BIAS_PI,
                  ridge_coef=Rdg_pi.coef_, C_EMP_PI=C_EMP_PI, DECAY_PI=DECAY_PI,
                  CP_PI=CP_PI, THRESH_PI=0.3, x_min=float(P_ALL[0]),
                  median_corr_pi=median_pi, DECAY_PI_CORRELATION_DEFINED=_pi_corr_defined)
    if verbose:
        print(f"  [calib TRAIN_END={TRAIN_END:,}] C_EMP_PN={C_EMP_PN:.3f} DECAY_PN={DECAY_pn:.8f} | "
              f"C_EMP_PI={C_EMP_PI:.4f} DECAY_PI={DECAY_PI:.8f}  ({time.time()-t_start:.1f}s)")
    return cal_pn, cal_pi

# ══════════════════════════════════════════════════════════════
# 7. jump_walk_pi -- independent memory
# ══════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════
# [PI-THRESH v17] Spectral RMS threshold derived from the calibrated PI
# Ridge correction. User-approved replacement of the legacy fixed 0.3.
# Theta_pi(x)^2 = median_corr_pi^2 + x/(2 log(x)^2)
#                 * sum_k (a_k^2 + b_k^2) / gamma_k^2.
# Independent of C_EMP, with units of prime counts. This is a spectral
# RMS model scale, NOT a proven bound for pi(x) error or cloud containment.
# The PN -> PI architecture and branch transition rules stay unchanged.
# ══════════════════════════════════════════════════════════════
def thresh_pi_at(cal_pi, x):
    x = float(x)
    if not math.isfinite(x) or x <= 1.0:
        raise ValueError('PI threshold requires finite x > 1')
    coef = np.asarray(cal_pi['ridge_coef'], dtype=np.float64)
    gamma = np.asarray(G60_PI, dtype=np.float64)
    if coef.ndim != 1 or len(coef) != 2*len(gamma):
        raise ValueError('PI Ridge coefficient count must equal 2 * zeta-zero count')
    if (not np.all(np.isfinite(coef)) or not np.all(np.isfinite(gamma))
            or np.any(gamma <= 0.0)):
        raise ValueError('PI spectral threshold requires finite coefficients and positive zeta zeros')
    median = float(cal_pi['median_corr_pi'])
    if not math.isfinite(median):
        raise ValueError('PI median correction must be finite')
    power = float(np.sum((coef[0::2]**2 + coef[1::2]**2) / (gamma**2),
                         dtype=np.float64))
    threshold_sq = median*median + x*power/(2.0*math.log(x)**2)
    if not math.isfinite(threshold_sq) or threshold_sq < 0.0:
        raise ArithmeticError('PI spectral threshold is not finite')
    return math.sqrt(threshold_sq)


def jump_walk_pi(x_start, x_end_target, cal_pi, off0=None, sf0=None, wo0=None,
                  MX_PI=1024, SEL_PI=512, PS_SCALE_PI=7,
                  CONF_THRESH_PI=None, PHASE1_CPS_PI=6, PHASE3_CPS_PI=6, ANOM_DECAY_PI=0.4):
    _check_D_range(x_end_target, "jump_walk_pi")   # [GUARD]
    if CONF_THRESH_PI is None: CONF_THRESH_PI = 1e9 if ANOMALY_OFF else 0.08   # [TEST-SWITCH]
    A_BIAS_PI=cal_pi['A_BIAS_PI']; B_BIAS_PI=cal_pi['B_BIAS_PI']; ridge_coef=cal_pi['ridge_coef']
    CP_PI=cal_pi['CP_PI']; DECAY_PI=cal_pi['DECAY_PI']
    # THRESH_PI is computed at every checkpoint by thresh_pi_at()   # [PI-THRESH]
    median_corr_pi = cal_pi.get('median_corr_pi', 0.0)

    x_cur = x_start
    off = np.zeros(1, float) if off0 is None else off0.copy()
    sf  = np.zeros(1, float) if sf0  is None else sf0.copy()
    cb  = np.zeros(1, int); ch = np.zeros(1, int)

    if wo0 is None:
        wo = np.full(len(off), (Fz_pi(x_start, A_BIAS_PI, B_BIAS_PI, ridge_coef)
                                - Fb_pi(x_start, A_BIAS_PI, B_BIAS_PI)
                                + median_corr_pi))
    else:
        wo = wo0.copy()

    anom_phase = np.zeros(len(off), dtype=int)
    anom_dir   = np.zeros(len(off), dtype=int)
    anom_amp   = np.zeros(len(off), dtype=float)
    anom_ph_cps= np.zeros(len(off), dtype=int)

    snapshots = []
    step = 0
    _max_steps = int(3 * max(x_end_target - x_cur, 1) / (CP_PI * math.log(max(x_cur, 3.)))) + 1000

    while x_cur < x_end_target and step < _max_steps:
        step += 1
        ln_x = math.log(x_cur)
        x_new = x_cur + CP_PI * ln_x
        if x_new > x_end_target: x_new = x_end_target
        dx_actual = x_new - x_cur
        if dx_actual <= 0.0: raise ArithmeticError('PI walk did not advance')
        fraction = dx_actual / (CP_PI * ln_x)

        Fb = Fb_pi(x_new, A_BIAS_PI, B_BIAS_PI)
        w = DECAY_PI * wo * math.exp(-fraction)
        vals = Fb + w + off + sf * (CP_PI * fraction)
        ae = np.where(anom_phase == 1, anom_dir * anom_amp,
                      np.where(anom_phase == 3, -anom_dir * anom_amp * ANOM_DECAY_PI, 0.0))

        snapshots.append((
            x_cur, x_new,
            float(np.min(off)), float(np.max(off)),
            float(np.min(sf)),  float(np.max(sf)),
            float(np.min(wo)),  float(np.max(wo)),
            float(np.min(ae)),  float(np.max(ae)),
            bool(np.max(np.abs(ae)) > 0.5),
            # MATH-FIX: preserve correlations inside each *retained* branch.
            # This does not restore branches discarded by original pruning.
            (off.copy(), sf.copy(), wo.copy(), ae.copy())
        ))

        tsn_pi = 1.0 / ln_x
        dFb_dx = _Rprime_table(x_new) + dbias_pi_dx(x_new, A_BIAS_PI, B_BIAS_PI)
        dWn = DECAY_PI * wo * (-1.0 / (CP_PI * ln_x)) * math.exp(-fraction)
        dfb = dFb_dx + sf / ln_x + dWn
        dev = np.abs(dfb - tsn_pi) / max(abs(tsn_pi), 1e-12)

        th  = np.full(len(off), 0.02)
        cf_ = np.full(len(off), 0.60)
        m3 = cb >= 3
        if m3.any():
            mv = np.minimum(cb[m3], 5)
            th[m3] = np.minimum(0.05, 0.02 * (1.5 ** mv))
            cf_[m3] = np.maximum(0.30, 0.60 * (0.9 ** cb[m3]))
        m5 = ch >= 5
        if m5.any(): th[m5] = np.maximum(0.005, th[m5] * (0.85 ** ch[m5]))

        ff = np.clip((tsn_pi - dfb) * cf_ * CP_PI, -1.5, 1.5)
        ps = np.abs(ff) * PS_SCALE_PI
        sn = np.clip(sf + ff, -1.5, 1.5)
        sr = np.clip(-ff, -1.5, 1.5)

        is_hold = dev <= th
        THRESH_PI = thresh_pi_at(cal_pi, x_cur)          # [PI-THRESH]
        is_brake = ~is_hold & (vals > Fb + THRESH_PI)
        is_jump = ~is_hold & ~is_brake & (vals < Fb - THRESH_PI)

        no, ns_, nc, nh = [], [], [], []
        nw = []
        na_phase, na_dir, na_amp, na_ph = [], [], [], []

        for j in range(len(off)):
            o, s_, c_, h_ = off[j], sf[j], int(cb[j]), int(ch[j])
            f_, p_, sn_, sr_ = ff[j], ps[j], sn[j], sr[j]
            wj = wo[j]
            ap, ad, aa, apc = anom_phase[j], anom_dir[j], anom_amp[j], anom_ph_cps[j]

            if is_hold[j]:
                no.append(o); ns_.append(f_ * 0.3); nc.append(0); nh.append(h_ + 1)
                nw.append(wj)
                na_phase.append(ap); na_dir.append(ad); na_amp.append(aa); na_ph.append(apc)
            elif is_brake[j]:
                no.append(o - p_); ns_.append(sn_); nc.append(c_ + 1); nh.append(0)
                nw.append(wj)
                na_phase.append(ap); na_dir.append(ad); na_amp.append(aa); na_ph.append(apc)
            elif is_jump[j]:
                no.append(o + p_); ns_.append(sn_); nc.append(c_ + 1); nh.append(0)
                nw.append(wj)
                na_phase.append(ap); na_dir.append(ad); na_amp.append(aa); na_ph.append(apc)
            else:
                o3 = (o + p_) if dfb[j] > tsn_pi else (o - p_)
                no.extend([o + p_, o - p_, o3])
                ns_.extend([sn_, sn_, sr_])
                nc.extend([c_ + 1] * 3)
                nh.extend([0] * 3)
                nw.extend([wj, wj, wj])
                na_phase.extend([ap, ap, ap])
                na_dir.extend([ad, ad, ad])
                na_amp.extend([aa, aa, aa])
                na_ph.extend([apc, apc, apc])

        off = np.array(no, float)
        sf  = np.array(ns_, float)
        cb  = np.array(nc, int)
        ch  = np.array(nh, int)
        wo  = np.array(nw, float)
        anom_phase = np.array(na_phase, int)
        anom_dir   = np.array(na_dir, int)
        anom_amp   = np.array(na_amp, float)
        anom_ph_cps= np.array(na_ph, int)

        if len(off) > MX_PI:
            idx = np.argsort(off, kind='stable')
            keep = np.round(np.linspace(0, len(idx) - 1, SEL_PI)).astype(int)
            off = off[idx[keep]]; sf = sf[idx[keep]]
            cb  = cb[idx[keep]];  ch = ch[idx[keep]]
            wo  = wo[idx[keep]]
            anom_phase = anom_phase[idx[keep]]
            anom_dir   = anom_dir[idx[keep]]
            anom_amp   = anom_amp[idx[keep]]
            anom_ph_cps= anom_ph_cps[idx[keep]]

        if len(off):
            off -= off.mean()

        target_wo = (Fz_pi(x_new, A_BIAS_PI, B_BIAS_PI, ridge_coef) - Fb
                     + median_corr_pi)
        wo = DECAY_PI * wo + (1.0 - DECAY_PI) * target_wo

        confidence = np.abs(wo) / max(THRESH_PI, 1e-12)
        for j in range(len(off)):
            if anom_phase[j] == 0:
                if confidence[j] >= CONF_THRESH_PI:
                    anom_phase[j] = 1
                    anom_dir[j] = 1 if wo[j] > 0 else -1
                    anom_amp[j] = abs(wo[j])
                    anom_ph_cps[j] = 1
            elif anom_phase[j] == 1:
                anom_ph_cps[j] += 1
                if anom_ph_cps[j] > PHASE1_CPS_PI:
                    anom_phase[j] = 2
                    anom_ph_cps[j] = 1
            elif anom_phase[j] == 2:
                anom_phase[j] = 3
                anom_ph_cps[j] = 1
            elif anom_phase[j] == 3:
                anom_ph_cps[j] += 1
                if anom_ph_cps[j] > PHASE3_CPS_PI:
                    anom_phase[j] = 0
                    anom_dir[j] = 0
                    anom_amp[j] = 0.0
                    anom_ph_cps[j] = 0

        x_cur = x_new
    if x_cur < x_end_target:
        raise ArithmeticError('PI walk exhausted its step budget before target endpoint')
    return snapshots

# ══════════════════════════════════════════════════════════════
# 8. make_pi_cloud_fn
# ══════════════════════════════════════════════════════════════
def make_pi_cloud_fn(snaps, cal_pi):
    if len(snaps)==0: return lambda x: (float('inf'), -float('inf'))
    PI_SNAP_X = np.array([s[1] for s in snaps])
    A_BIAS_PI=cal_pi['A_BIAS_PI']; B_BIAS_PI=cal_pi['B_BIAS_PI']
    DECAY_PI=cal_pi['DECAY_PI']; CP_PI=cal_pi['CP_PI']
    @lru_cache(maxsize=200000)
    def _cached(x_query):
        if x_query < snaps[0][0] or x_query > PI_SNAP_X[-1]:
            return (float('inf'), -float('inf'))
        idx=int(np.searchsorted(PI_SNAP_X,x_query)); idx=min(max(idx,0),len(snaps)-1)
        snap = snaps[idx]
        (x_prev, x_new, off_min, off_max, sf_min, sf_max,
         wo_min, wo_max, ae_min, ae_max, ae_present) = snap[:11]
        if x_query < x_prev: return (float('inf'), -float('inf'))
        frac=(x_query-x_prev)/max(x_new-x_prev,1e-9); frac=min(max(frac,0.0),1.0)
        delta_k=frac*CP_PI
        Fb_q=Fb_pi(x_query,A_BIAS_PI,B_BIAS_PI)
        decay = DECAY_PI*math.exp(-delta_k/CP_PI)
        if len(snap) >= 12:
            off_b, sf_b, wo_b, ae_b = snap[11]
            values = Fb_q + decay*wo_b + off_b + sf_b*delta_k
            if ae_present:
                values = np.concatenate((values,values+ae_b))
            # Exact over retained branches, not a rectangle made of unrelated extremes.
            return float(np.min(values)), float(np.max(values))
        # Compatibility for original 11-field snapshots. The interval-hull is
        # conservative for the retained branches but not necessarily tight.
        vals_min = Fb_q + decay*wo_min + off_min + sf_min*delta_k
        vals_max = Fb_q + decay*wo_max + off_max + sf_max*delta_k
        if ae_present:
            vals_min = min(vals_min, vals_min + ae_min)
            vals_max = max(vals_max, vals_max + ae_max)
        return float(vals_min), float(vals_max)
    return _cached

def find_boundary(n, x_center, direction, cloud_fn, x_step_init=50, max_expand=50,
                   max_extra_steps=500, x_limit=None):
    x_out=None; step=x_step_init; expand_count=0; last_inside=x_center
    def _past_limit(x):
        return x_limit is not None and ((direction>0 and x>x_limit) or (direction<0 and x<x_limit))
    while expand_count<max_expand:
        x_try=x_center+direction*step
        if _past_limit(x_try): return last_inside, False
        if x_try<=2.0: break
        cmin,cmax=cloud_fn(x_try)
        if cmin<=n<=cmax: last_inside=x_try; step*=1.5; expand_count+=1
        else:
            if x_out is None: x_out=x_try
            step*=1.5; expand_count+=1; break
    if x_out is None: return last_inside, False
    extra_step=step; x_current=x_out; found_again=False
    for _ in range(max_extra_steps):
        x_try=x_current+direction*extra_step
        if _past_limit(x_try): break
        if x_try<=2.0: break
        cmin,cmax=cloud_fn(x_try)
        if cmin<=n<=cmax: x_after_gap=x_try; found_again=True; break
        x_current=x_try; extra_step*=1.5
    if found_again:
        x1,x2=x_out,x_after_gap
        for _ in range(40):
            mid=(x1+x2)/2; cmin,cmax=cloud_fn(mid)
            if cmin<=n<=cmax: x2=mid
            else: x1=mid
        return x2,True
    else:
        x1,x2=last_inside,x_out
        for _ in range(40):
            mid=(x1+x2)/2; cmin,cmax=cloud_fn(mid)
            if cmin<=n<=cmax: x1=mid
            else: x2=mid
        return x1,False


# ══════════════════════════════════════════════════════════════
# MATHEMATICAL ENDPOINT REPAIR: all admissible INTEGER positions in the
# theorem interval and the actually evaluated PI snapshot domain.
# ══════════════════════════════════════════════════════════════
def pi_integer_components(n_target, pi_cloud_at, lower, upper, *, domain=None):
    """Enumerate EVERY integer x in [lower,upper] with pi_lo(x)<=n<=pi_hi(x).

    This proves *enumeration completeness on the checked integer grid*
    relative to pi_cloud_at. It does NOT certify the empirical PI band,
    recover pruned branches, or prove a true-prime containment guarantee.
    No seed/root heuristics; disconnected components are retained.
    """
    n_target=_require_exact_int(n_target,'n_target',1)
    low=int(math.ceil(lower)); high=int(math.floor(upper))
    if domain is not None:
        low=max(low,int(math.ceil(domain[0])))
        high=min(high,int(math.floor(domain[1])))
    if high<low:return []
    if high > (1 << 53):
        raise ValueError('PI integer grid exceeds exact float64 integer representation; cannot certify enumeration')
    components=[]
    running_start=None
    prev=None
    for x in range(low,high+1):
        l,r=pi_cloud_at(float(x))
        included=(math.isfinite(l) and math.isfinite(r) and l<=n_target<=r)
        if included:
            if running_start is None:running_start=x
            prev=x
        elif running_start is not None:
            components.append((running_start,prev))
            running_start=None
    if running_start is not None:components.append((running_start,prev))
    return components

def _pi_components_seeded_by_pn(pi_components, accepted_pn_candidates):
    """Retain ONLY integer PI components touching an accepted PN location.

    A PN location is a float: its outward integer bracket [floor(v),ceil(v)]
    is used solely for discrete-grid membership.  We do not treat PN seeds as
    endpoints and we do not include any unseeded PI components.
    """
    brackets=[(int(math.floor(v)),int(math.ceil(v)))
              for v in accepted_pn_candidates if math.isfinite(float(v))]
    return [(a,b) for a,b in pi_components
            if any(lo <= b and hi >= a for lo,hi in brackets)]


# ══════════════════════════════════════════════════════════════
# POST-RAW PN->PI: PN selects PI components; final accepted set is their contiguous hull
# ══════════════════════════════════════════════════════════════
def _merge_integer_components(components, lower, upper):
    """Clip, sort and merge overlapping/adjacent closed integer intervals."""
    items=[]
    lo_bound=int(lower); hi_bound=int(upper)
    for a,b in components:
        a=max(lo_bound,int(a)); b=min(hi_bound,int(b))
        if a<=b: items.append((a,b))
    if not items: return []
    items.sort()
    out=[list(items[0])]
    for a,b in items[1:]:
        if a<=out[-1][1]+1:
            if b>out[-1][1]: out[-1][1]=b
        else:
            out.append([a,b])
    return [(a,b) for a,b in out]


def _pn_point_components(points, radius, lower, upper, *, raw_brackets=False):
    """Represent PN branch points on the integer x-grid without filling gaps.

    Raw diagnostics use the outward [floor(v), ceil(v)] bracket for each
    floating branch point. Post-margin support uses the exact integer points
    in [v-radius, v+radius], rounded outward by one float ULP before ceil/floor.
    """
    intervals=[]
    radius=float(radius)
    for value in np.asarray(points,dtype=np.float64).reshape(-1):
        v=float(value)
        if not math.isfinite(v):
            continue
        if raw_brackets:
            a=int(math.floor(v)); b=int(math.ceil(v))
        else:
            left=math.nextafter(v-radius,-math.inf)
            right=math.nextafter(v+radius,math.inf)
            a=int(math.ceil(left)); b=int(math.floor(right))
        if a<=b:
            intervals.append((a,b))
    return _merge_integer_components(intervals,lower,upper)


def _pi_raw_and_empirical_components(n_target, pi_cloud_at, lower, upper,
                                     c_emp_pi, multiplier=1.0):
    """Scan the same integer domain for RAW PI and post-RAW PI-margin sets.

    RAW membership is evaluated first and is exactly
        pi_raw_low(x) <= n_target <= pi_raw_high(x).
    Only the separate empirical predicate uses C_EMP_PI:
        pi_raw_low(x)-M_PI(x) <= n_target <= pi_raw_high(x)+M_PI(x),
        M_PI(x)=multiplier*C_EMP_PI*log(x)^2.
    The returned component lists preserve every disconnected run.
    """
    n_target=_require_exact_int(n_target,'n_target',1)
    low=int(math.ceil(lower)); high=int(math.floor(upper))
    if high<low: return [],[]
    if high>(1<<53):
        raise ValueError('PI integer grid exceeds exact float64 integer representation')
    c=float(c_emp_pi); mult=float(multiplier)
    if not (math.isfinite(c) and c>=0.0 and math.isfinite(mult) and mult>=0.0):
        raise ValueError('C_EMP_PI and margin multiplier must be finite and nonnegative')
    raw_runs=[]; emp_runs=[]
    raw_start=None; emp_start=None
    raw_prev=None; emp_prev=None
    for x in range(low,high+1):
        pi_lo,pi_hi=pi_cloud_at(float(x))
        finite=(math.isfinite(pi_lo) and math.isfinite(pi_hi) and pi_lo<=pi_hi)
        # RAW predicate: deliberately contains no C_EMP term.
        raw_ok=bool(finite and pi_lo<=n_target<=pi_hi)
        if raw_ok:
            if raw_start is None: raw_start=x
            raw_prev=x
        elif raw_start is not None:
            raw_runs.append((raw_start,raw_prev)); raw_start=None; raw_prev=None
        # Post-RAW PI dilation in pi-count units. This is heuristic, not a theorem.
        delta=mult*c*_NORM_S(float(x)) if finite else 0.0
        emp_ok=bool(finite and pi_lo-delta<=n_target<=pi_hi+delta)
        if emp_ok:
            if emp_start is None: emp_start=x
            emp_prev=x
        elif emp_start is not None:
            emp_runs.append((emp_start,emp_prev)); emp_start=None; emp_prev=None
    if raw_start is not None: raw_runs.append((raw_start,raw_prev))
    if emp_start is not None: emp_runs.append((emp_start,emp_prev))
    return raw_runs,emp_runs


def _joint_post_raw_geometry(n_target, pn_raw_candidates, pi_raw, theorem_bounds,
                             snapshot_domain, cal_pn, cal_pi,
                             margin_multiplier=1.0):
    """Use raw PN candidates as seeds to select complete PI integer components.

    First pass each raw PN branch candidate x_i through the *raw* PI cloud.
    The raw candidates accepted by PI are retained in diagnostics. In the
    post-RAW phase, evaluate the PI empirical band at the SAME raw candidate
    locations (no PN support dilation), and keep every full integer component
    of PI's empirical admissible set containing an accepted candidate's
    outward [floor(x_i),ceil(x_i)] integer bracket.

    PI empirical membership is
      pi_raw_low(x)-C_EMP_PI*log(x)^2 <= n <= pi_raw_high(x)+C_EMP_PI*log(x)^2.
    There is NO set intersection or endpoint clipping using a PN envelope.
    The theorem domain limits the PI scan independently. The FINAL
    admissible set is one contiguous closed integer interval from the
    smallest selected PI lower bound to the largest selected PI upper bound.
    EVERY integer in between, including gaps between selected PI components,
    belongs to the final containment set. Original disjoint PI components
    remain available strictly as diagnostics. All containment claims are empirical.
    """
    n_target=_require_exact_int(n_target,'n_target',1)
    Ln,Un=map(int,theorem_bounds)
    if snapshot_domain is None:
        raise ValueError('PI snapshot domain is empty; cannot construct D_n')
    domain_lo=max(int(math.ceil(Ln)),int(math.ceil(snapshot_domain[0])),3)
    domain_hi=min(int(math.floor(Un)),int(math.floor(snapshot_domain[1])))
    if domain_hi<domain_lo:
        return dict(pi_integer_search_domain=(domain_lo,domain_hi),
                    raw_pn_integer_bracket_components=[],raw_pi_integer_components=[],
                    raw_integer_bracket_intersection_components=[],raw_intersection_components=[],
                    pn_empirical_integer_components=[],pi_empirical_integer_components=[],
                    empirical_integer_components=[],selected_pi_integer_components=[],
                    final_integer_components=[],final_interval=None,
                    final_included_integer_count=0,final_fills_internal_gaps=True,
                    raw_compatible_pn_candidates=[],
                    raw_compatible=False,pi_integer_scan_complete=True,
                    fallback_used=True,fallback_reason='EMPTY_THEOREM_SNAPSHOT_DOMAIN')

    points=np.asarray(pn_raw_candidates,dtype=np.float64).reshape(-1)
    if points.size==0 or not np.all(np.isfinite(points)):
        raise ArithmeticError('PN raw candidate cloud is empty or non-finite')

    # RAW PN/PI compatibility is tested before reading either C_EMP value.
    raw_compatible_points=[]
    for v in points:
        if domain_lo<=v<=domain_hi:
            lo,hi=pi_raw(float(v))
            if math.isfinite(lo) and math.isfinite(hi) and lo<=n_target<=hi:
                raw_compatible_points.append(float(v))
    raw_compatible_points=sorted(set(raw_compatible_points))
    raw_accepted_pn_point_components=[(v,v) for v in raw_compatible_points]
    raw_interval=((raw_compatible_points[0],raw_compatible_points[-1])
                  if raw_compatible_points else None)

    # PN point brackets are diagnostic only, never used to cut PI boundaries.
    raw_pn_integer_brackets=_pn_point_components(
        points,0.0,domain_lo,domain_hi,raw_brackets=True)

    # Post-RAW: PI margin only. PN's empirical multiplier is deliberately
    # NOT used to widen candidates or to clip PI-component endpoints.
    c_pi=cal_pi.get('C_EMP_PI')
    if c_pi is None:
        raise ValueError('C_EMP_PI is required for the PI empirical band')
    c_pi=float(c_pi); mult=float(margin_multiplier)
    if not (math.isfinite(c_pi) and c_pi>=0.0
            and math.isfinite(mult) and mult>=0.0):
        raise ValueError('C_EMP_PI and margin_multiplier must be finite and nonnegative')

    # One integer scan extracts all (possibly disconnected) PI components.
    raw_pi_components,pi_empirical_components=_pi_raw_and_empirical_components(
        n_target,pi_raw,domain_lo,domain_hi,c_pi,mult)

    # PN->PI: use raw PN float candidates for *both* tests, not PN dilations.
    # The PI margin acts in prime-count units, evaluated at candidate x_i.
    empirical_compatible_points=[]
    for v in points:
        v=float(v)
        if domain_lo<=v<=domain_hi:
            lo,hi=pi_raw(v)
            if math.isfinite(lo) and math.isfinite(hi) and lo<=hi:
                delta=mult*c_pi*_NORM_S(v)
                if lo-delta<=n_target<=hi+delta:
                    empirical_compatible_points.append(v)
    empirical_compatible_points=sorted(set(empirical_compatible_points))

    # Select WHOLE PI components touched by admissible PN seeds.
    # PN candidates are selectors, not independent geometric endpoints.
    raw_seeded_pi_components=_pi_components_seeded_by_pn(
        raw_pi_components,raw_compatible_points)
    empirical_components=_pi_components_seeded_by_pn(
        pi_empirical_components,empirical_compatible_points)
    accepted_pn=len(empirical_compatible_points)

    # FINAL DOMAIN (not just two display numbers): close all internal gaps.
    # Each integer L <= x <= U is included, even if it fails PI membership
    # between two selected disjoint PI components. The original components
    # are retained separately ONLY for diagnostic inspection.
    empirical_interval=((min(a for a,b in empirical_components),
                         max(b for a,b in empirical_components))
                        if empirical_components else None)
    final_integer_components=([empirical_interval] if empirical_interval is not None else [])
    final_included_integer_count=(empirical_interval[1]-empirical_interval[0]+1
                                  if empirical_interval is not None else 0)
    return dict(
        pi_integer_search_domain=(domain_lo,domain_hi),
        pi_integer_scan_complete=True,
        raw_pn_candidates=points.tolist(),
        raw_pn_integer_bracket_components=raw_pn_integer_brackets,
        raw_pi_integer_components=raw_pi_components,
        raw_integer_bracket_intersection_components=raw_seeded_pi_components,  # legacy diagnostic alias
        raw_seeded_pi_integer_components=raw_seeded_pi_components,
        raw_intersection_components=raw_accepted_pn_point_components,  # legacy diagnostic alias
        raw_accepted_pn_point_components=raw_accepted_pn_point_components,
        raw_interval=raw_interval,
        raw_compatible_pn_candidates=raw_compatible_points,
        raw_compatible=bool(raw_compatible_points),
        pn_margin_in_x_units=0.0,
        pn_margin_used_for_bounds=False,
        pi_margin_formula='C_EMP_PI * log(x)^2 in pi-count units',
        pn_margin_formula='UNUSED: raw PN candidates are seeds, not bound constraints',
        pn_empirical_integer_components=[],  # legacy field; no PN dilation
        empirical_compatible_pn_candidates=empirical_compatible_points,
        pi_empirical_integer_components=pi_empirical_components,
        selected_pi_integer_components=empirical_components,  # diagnostics only
        empirical_integer_components=final_integer_components,  # authoritative continuous set
        final_integer_components=final_integer_components,
        empirical_interval=empirical_interval,
        final_interval=empirical_interval,
        final_included_integer_count=final_included_integer_count,
        final_fills_internal_gaps=True,
        n_pi_integer_components_raw=len(raw_pi_components),
        n_pi_integer_components_empirical=len(pi_empirical_components),
        n_raw_intersection_components=len(raw_accepted_pn_point_components),  # legacy alias
        n_raw_integer_bracket_intersection_components=len(raw_seeded_pi_components),  # legacy alias
        n_empirical_intersection_components=len(empirical_components),  # legacy alias
        n_raw_seeded_pi_components=len(raw_seeded_pi_components),
        n_selected_pi_components=len(empirical_components),
        n_raw_valid_after_pi_filter=len(raw_compatible_points),
        n_valid_after_pi_filter=accepted_pn,
        n_empirical_pn_candidates=int(len(points)),
        raw_cloud_uses_c_emp=False,
        c_emp_applied_only_after_raw_extraction=True,
        all_disconnected_components_preserved=True,  # diagnostic PI components only
        selected_pi_components_preserved_for_diagnostics=True,
        returned_interval_is_outer_envelope=True,
        margin_multiplier=float(mult),
        fallback_used=not bool(empirical_components),
        fallback_reason=None if empirical_components else 'NO_PI_COMPONENT_SELECTED_BY_PN_CANDIDATES',
        endpoint_source='FULL_PI_COMPONENTS_SELECTED_BY_PN_CANDIDATES',
    )

# ══════════════════════════════════════════════════════════════
# 9. jump_walk_pn -- independent memory + dimension fix
# ══════════════════════════════════════════════════════════════
def jump_walk_pn(n_start, n_end_target, cal_pn, off0=None, sf0=None):
    CP, MX, SEL = 40, 128, 64
    c_vec=cal_pn['c_vec']; DECAY=cal_pn['DECAY']
    THRESH_PN=cal_pn['THRESH_PN']; TSN_BASE=cal_pn['TSN_BASE']; TRAIN_END=cal_pn['TRAIN_END']
    ridge_coef_pn = cal_pn['ridge_coef_pn']; dc_dn_fn = cal_pn['dc_dn']
    median_corr_pn = cal_pn.get('median_corr_pn', 0.0)

    off=np.zeros(1,float) if off0 is None else off0.copy()
    sf=np.zeros(1,float) if sf0 is None else sf0.copy()
    cb=np.zeros(1,int); ch=np.zeros(1,int)

    li_start = inv_R_scalar(float(n_start))  # [UPGRADE-R]
    feat_start = build_pn_features(np.array([li_start]))[0]
    wo = np.full(len(off), float(np.dot(feat_start, ridge_coef_pn)) + median_corr_pn)

    anom_phase = np.zeros(len(off), dtype=int)
    anom_dir   = np.zeros(len(off), dtype=int)
    anom_amp   = np.zeros(len(off), dtype=float)
    anom_ph_cps= np.zeros(len(off), dtype=int)

    lcp=n_start
    # [TEST-SWITCH] ANOMALY_OFF=True disables the anomaly mechanism to
    # measure its effect. Default False -- i.e. model behaviour
    # unchanged.
    PHASE1_CPS,PHASE3_CPS,CONF_THRESH = 2,2,(1e9 if ANOMALY_OFF else 0.30)
    av_hist=[]
    n_range = np.arange(n_start+1, n_end_target+1, dtype=np.float64)
    li_range = inv_R_vec(n_range)            # [UPGRADE-R]
    c_range = c_vec(n_range)
    n = n_start
    for _i in range(len(n_range)):
        n += 1
        li_n = li_range[_i]
        Fb = li_n + c_range[_i]
        delta=float(n-lcp)
        w = DECAY * wo * math.exp(-delta/CP)
        vals = Fb + w + off + sf*delta

        ln_n=math.log(n)
        ae = np.where(anom_phase==1, anom_dir*anom_amp,
                      np.where(anom_phase==3, -anom_dir*anom_amp*0.4, 0.0))

        candidates = []
        for j in range(len(off)):
            base = Fb + w[j] + off[j] + sf[j]*delta
            candidates.append(base)
            if abs(ae[j]) > 1.0:
                candidates.append(base + ae[j])
                candidates.append(base - ae[j])
        av = np.array(candidates, dtype=float)
        av_hist.append((n, av.copy()))

        if n % CP == 0:
            tsn=TSN_BASE+(n-TRAIN_END)/TRAIN_END
            dWn=DECAY*wo*(-1.0/CP)*math.exp(-delta/CP)
            dfb = (1.0 / _Rprime_table(li_n)) + dc_dn_fn(n) + sf + dWn
            dev=np.abs(dfb-tsn)/max(abs(tsn),1e-9)
            th=np.full(len(off),0.02); cf_=np.full(len(off),0.60)
            m3=cb>=3
            if m3.any():
                mv=np.minimum(cb[m3],5); th[m3]=np.minimum(0.05,0.02*(1.5**mv)); cf_[m3]=np.maximum(0.30,0.60*(0.9**cb[m3]))
            m5=ch>=5
            if m5.any(): th[m5]=np.maximum(0.005,th[m5]*(0.85**ch[m5]))
            ff=np.clip((tsn-dfb)*cf_,-1.5,1.5); ps=np.abs(ff)*50
            sn=np.clip(sf+ff,-1.5,1.5); sr=np.clip(-ff,-1.5,1.5)
            is_hold=(dev<=th); is_brake=~is_hold&(vals>Fb+THRESH_PN); is_jump=~is_hold&~is_brake&(vals<Fb-THRESH_PN)

            no,ns_,nc,nh=[],[],[],[]
            nw=[]
            na_phase,na_dir,na_amp,na_ph=[],[],[],[]
            for j in range(len(off)):
                o,s_,c_,h_=off[j],sf[j],int(cb[j]),int(ch[j])
                f_,p_,sn_,sr_=ff[j],ps[j],sn[j],sr[j]
                wj=wo[j]
                ap,ad,aa,apc=anom_phase[j],anom_dir[j],anom_amp[j],anom_ph_cps[j]
                if is_hold[j]:
                    no.append(o); ns_.append(f_*0.3); nc.append(0); nh.append(h_+1)
                    nw.append(wj)
                    na_phase.append(ap); na_dir.append(ad); na_amp.append(aa); na_ph.append(apc)
                elif is_brake[j]:
                    no.append(o-p_); ns_.append(sn_); nc.append(c_+1); nh.append(0)
                    nw.append(wj)
                    na_phase.append(ap); na_dir.append(ad); na_amp.append(aa); na_ph.append(apc)
                elif is_jump[j]:
                    no.append(o+p_); ns_.append(sn_); nc.append(c_+1); nh.append(0)
                    nw.append(wj)
                    na_phase.append(ap); na_dir.append(ad); na_amp.append(aa); na_ph.append(apc)
                else:
                    o3=(o+p_) if dfb[j]>tsn else (o-p_)
                    no.extend([o+p_,o-p_,o3]); ns_.extend([sn_,sn_,sr_])
                    nc.extend([c_+1]*3); nh.extend([0]*3)
                    nw.extend([wj,wj,wj])
                    na_phase.extend([ap,ap,ap]); na_dir.extend([ad,ad,ad])
                    na_amp.extend([aa,aa,aa]); na_ph.extend([apc,apc,apc])

            off=np.array(no,float); sf=np.array(ns_,float)
            cb=np.array(nc,int); ch=np.array(nh,int)
            wo=np.array(nw,float)
            anom_phase=np.array(na_phase,int); anom_dir=np.array(na_dir,int)
            anom_amp=np.array(na_amp,float); anom_ph_cps=np.array(na_ph,int)

            if len(off)>MX:
                idx=np.argsort(off); ki=np.round(np.linspace(0,len(idx)-1,SEL)).astype(int)
                off=off[idx[ki]]; sf=sf[idx[ki]]; cb=cb[idx[ki]]; ch=ch[idx[ki]]
                wo=wo[idx[ki]]
                anom_phase=anom_phase[idx[ki]]; anom_dir=anom_dir[idx[ki]]
                anom_amp=anom_amp[idx[ki]]; anom_ph_cps=anom_ph_cps[idx[ki]]

            if len(off): off -= off.mean()

            ridge_pred = float(np.dot(build_pn_features(np.array([li_n]))[0], ridge_coef_pn))
            target_wo = ridge_pred + median_corr_pn
            wo = DECAY * wo + (1.0 - DECAY) * target_wo

            confidence=abs(wo)/(max(abs(THRESH_PN),1e-9))
            for j in range(len(off)):
                if anom_phase[j]==0:
                    if confidence[j]>=CONF_THRESH:
                        anom_phase[j]=1; anom_dir[j]=1 if wo[j]>0 else -1
                        anom_amp[j]=abs(wo[j]); anom_ph_cps[j]=1
                elif anom_phase[j]==1:
                    anom_ph_cps[j]+=1
                    if anom_ph_cps[j]>PHASE1_CPS: anom_phase[j]=2; anom_ph_cps[j]=1
                elif anom_phase[j]==2: anom_phase[j]=3; anom_ph_cps[j]=1
                elif anom_phase[j]==3:
                    anom_ph_cps[j]+=1
                    if anom_ph_cps[j]>PHASE3_CPS:
                        anom_phase[j]=0; anom_dir[j]=0; anom_amp[j]=0.0; anom_ph_cps[j]=0

            lcp=n
    return av_hist

# ══════════════════════════════════════════════════════════════
# 10. jump_predict
# ══════════════════════════════════════════════════════════════
def jump_predict(n_target, cal_pi, cal_pn, settle_steps=800, margin_multiplier=1.0):
    """PN raw candidates -> PI tests -> full selected PI component endpoints.

    PI alone defines the empirical-component boundaries. PN candidates seed
    the components to select; they never intersect or clip the PI components.
    All disconnected components are preserved in diagnostics. This remains
    a heuristic interval, not a proven guarantee of prime containment.
    """
    n_target=_require_exact_int(n_target,'n_target',1)
    settle_steps=_require_exact_int(settle_steps,'settle_steps',1)
    multiplier=float(margin_multiplier)
    if not math.isfinite(multiplier) or multiplier<0.0:
        raise ValueError('margin_multiplier must be finite and nonnegative')

    x_target_raw=inv_R_scalar(float(n_target))
    x_target=x_target_raw+cal_pn['c_vec'](np.array([float(n_target)]))[0]
    n_jump_start_pn=max(cal_pn['TRAIN_END']+1,n_target-settle_steps*40)
    n_jump_start_pn=min(n_target,((n_jump_start_pn+39)//40)*40)
    av_hist_pn=jump_walk_pn(n_jump_start_pn,n_target,cal_pn)
    pn_raw=av_hist_pn[-1][1] if av_hist_pn else np.array([x_target],dtype=np.float64)
    pn_raw=np.asarray(pn_raw,dtype=np.float64)
    if pn_raw.size==0 or not np.all(np.isfinite(pn_raw)):
        raise ArithmeticError('PN raw candidate cloud is empty or non-finite')

    ln_x=math.log(max(x_target,3.0))
    # PI raw walk domain is independent of C_EMP.
    walk_radius=settle_steps*cal_pi['CP_PI']*ln_x
    x_jump_start=max(x_target-walk_radius,cal_pi.get('x_min',1000.0))
    x_end=x_target+walk_radius
    snaps=jump_walk_pi(x_jump_start,x_end,cal_pi)
    if not snaps:
        raise ArithmeticError('PI raw walk returned no snapshots')
    pi_raw=make_pi_cloud_fn(snaps,cal_pi)
    Ln,Un=certified_nth_prime_bounds(n_target)
    snapshot_domain=(snaps[0][0],snaps[-1][1])

    diag=dict(n_candidates_pn=int(len(pn_raw)),pn_raw=(float(np.min(pn_raw)),float(np.max(pn_raw))),
              pn_raw_candidates=pn_raw.tolist(),pi_raw_uses_c_emp=False,pn_raw_uses_c_emp=False,
              axler_lower=Ln,axler_upper=Un,raw_pi_snapshots=(x_jump_start,x_end),
              raw_contains_no_certification=True)
    geom=_joint_post_raw_geometry(n_target,pn_raw,pi_raw,(Ln,Un),snapshot_domain,
                                  cal_pn,cal_pi,margin_multiplier=multiplier)
    diag.update(geom)
    diag['raw_interval']=geom.get('raw_interval')
    diag['raw_integer_components']=geom.get('raw_seeded_pi_integer_components',[])
    diag['raw_accepted_pn_point_components']=geom.get('raw_accepted_pn_point_components',[])
    diag['pi_integer_components_before_pn_filter']=geom.get('n_pi_integer_components_empirical',0)
    diag['n_valid_after_pi_filter']=geom.get('n_valid_after_pi_filter',0)
    diag['n_empirical_pn_candidates']=geom.get('n_empirical_pn_candidates',0)
    diag['c_emp_pn_applied_after_raw']=False
    diag['c_emp_pi_applied_after_raw']=True
    diag['c_emp_applied_only_after_raw_extraction']=True
    diag['extract_pctl_ignored_for_component_completeness']=True
    diag['pi_integer_scan_complete']=geom.get('pi_integer_scan_complete',False)

    components=geom.get('empirical_integer_components',[])
    if not components:
        diag['fallback_used']=True
        diag['fallback_reason']=geom.get('fallback_reason') or 'NO_PI_COMPONENT_SELECTED_BY_PN_CANDIDATES'
        diag['n_valid_after_pi_filter']=0
        # The theorem fallback is also ONE contiguous final containment set.
        diag['final_integer_components']=[(Ln,Un)]
        diag['final_interval']=(Ln,Un)
        diag['final_included_integer_count']=Un-Ln+1
        diag['final_interval_source']='THEOREM_FALLBACK'
        return Ln,Un,diag

    x1,x2=geom['empirical_interval']
    diag['fallback_used']=False
    diag['final_interval_source']='CONTIGUOUS_HULL_OF_PI_COMPONENTS_SELECTED_BY_PN'
    diag.update(x_center=float(np.median(pn_raw)),extract_pctl=0.0,
                endpoint_source='CONTIGUOUS_HULL_OF_PI_COMPONENTS_SELECTED_BY_PN',
                cloud_lower=x1,cloud_upper=x2,clipped_low=False,clipped_high=False)
    return int(x1),int(x2),diag

# ══════════════════════════════════════════════════════════════
# 11. safe_jump_predict_lowmem
# ══════════════════════════════════════════════════════════════
def safe_jump_predict_lowmem(n_target, cal_pn, cal_pi, P_ALL, ratio=2.0, settle_steps=800, verbose=True):
    ratio = float(ratio)
    if not math.isfinite(ratio) or ratio <= 1.0:
        raise ValueError('ratio must be finite and > 1 to ensure forward progress')
    stages = []
    cur_pn, cur_pi, cur_P = cal_pn, cal_pi, P_ALL
    while n_target > ratio * cur_pn['TRAIN_END']:
        new_train_end = int(ratio * cur_pn['TRAIN_END'])
        if new_train_end <= cur_pn['TRAIN_END']:
            raise ArithmeticError('recalibration failed to advance TRAIN_END')
        if verbose: print(f"  target {n_target:,} is beyond {ratio}x the current TRAIN_END ({cur_pn['TRAIN_END']:,}) "
                           f"-- recalibrating at {new_train_end:,}...")
        t0=time.time()
        del cur_P
        gc.collect()
        cur_P = gen_primes(new_train_end + 5000)
        cur_pn, cur_pi = compute_calibration(new_train_end, cur_P, verbose=verbose)
        stages.append((new_train_end, time.time()-t0))
    x1, x2, diag = jump_predict(n_target, cur_pi, cur_pn, settle_steps=settle_steps)
    diag['recalibration_stages'] = stages
    diag['final_TRAIN_END'] = cur_pn['TRAIN_END']
    return x1, x2, diag

# ══════════════════════════════════════════════════════════════
# 12. Direct execution
# ══════════════════════════════════════════════════════════════


# ══════════════════════════════════════════════════════════════
# [PORT-2] Full walk path -- with calibration **shared** with the jump
# ══════════════════════════════════════════════════════════════
# It used to live in a separate file with its own calibration, so every
# fix had to be applied twice -- the source of most silent divergences.
# Both now call the same compute_calibration, so the walk inherits
# automatically: the R^-1 base, the full C_EMP scan, parallelism, and all
# speedups.
PER_BRANCH_EXTREMES = False   # enabling it multiplies memory ~8-9x over millions of points

def run_pn_cloud(cal, TRAIN_END, TEST_START, TEST_END, P_ALL, return_cloud=True):
    CP, MX, SEL = 40, 128, 64
    # MATH-FIX: propagate over EVERY rank after training, but output
    # predictions only for TEST_START..TEST_END; never skip hidden states.
    n_walk = np.arange(TRAIN_END+1, TEST_END+1, dtype=np.float64)
    li_te = inv_R_vec(n_walk)
    Fb_te = li_te + cal['c_vec'](n_walk)

    wo = np.full(1, float(np.dot(build_pn_features(np.array([inv_R_scalar(float(TRAIN_END))]))[0], cal['ridge_coef_pn'])) + cal['median_corr_pn'])
    lcp = TRAIN_END
    off = np.zeros(1, float); sf = np.zeros(1, float)
    cb = np.zeros(1, int); ch = np.zeros(1, int)
    N = TEST_END - TRAIN_END
    av_list = []
    # [WAVE-CLEANUP] anomaly mechanism removed -- its effect on the walk
    # was exactly 0.00%
    TSN_BASE = cal['TSN_BASE']; THRESH_PN = cal['THRESH_PN']
    ridge_coef = cal['ridge_coef_pn']; median_corr = cal['median_corr_pn']; DECAY = cal['DECAY']

    for i in range(N):
        n = TRAIN_END + 1 + i; Fb = float(Fb_te[i]); delta = float(n - lcp)
        w = DECAY * wo * math.exp(-delta / CP)          # now an array (was scalar)
        vals = Fb + w + off + sf * delta
        ln_n = math.log(n)

        # the global anomaly state remains a shared scalar -- it was not
        # the source of the off.mean() problem, so it is left as is
        ae = 0.0

        if PER_BRANCH_EXTREMES:
            candidates = []
            for j in range(len(off)):
                base = float(vals[j])
                candidates.append(base)
                if abs(ae) > 1.0:
                    candidates.append(base + ae); candidates.append(base - ae)
            av = np.array(candidates, dtype=float)
        else:
            # a fixed set of extremes (6 values) around a representative w
            # = the mean of wo across branches, instead of repeating them
            # per branch (keeps av_list the same size as the original)
            w_rep = DECAY * float(wo.mean()) * math.exp(-delta / CP)
            extreme_vals = np.array([], dtype=float)
            if abs(ae) > 1.:
                av = np.concatenate([vals, vals+ae, vals-ae, extreme_vals])
            else:
                av = np.concatenate([vals, extreme_vals])

        if return_cloud and n >= TEST_START: av_list.append(av.copy())

        if n % CP == 0:
            tsn = TSN_BASE + (n - TRAIN_END)/TRAIN_END
            dWn = DECAY * wo * (-1.0/CP) * math.exp(-delta/CP)
            li_n = li_te[i]
            dfb = (1.0 / _Rprime_table(li_n)) + cal['dc_dn'](n) + sf + dWn
            dev = np.abs(dfb - tsn) / max(abs(tsn), 1e-9)
            th = np.full(len(off), 0.02); cf_ = np.full(len(off), 0.60)
            m3 = cb >= 3
            if m3.any():
                mv = np.minimum(cb[m3], 5)
                th[m3] = np.minimum(0.05, 0.02*(1.5**mv)); cf_[m3] = np.maximum(0.30, 0.60*(0.9**cb[m3]))
            m5 = ch >= 5
            if m5.any(): th[m5] = np.maximum(0.005, th[m5]*(0.85**ch[m5]))
            ff = np.clip((tsn-dfb)*cf_, -1.5, 1.5); ps = np.abs(ff)*50
            sn = np.clip(sf+ff, -1.5, 1.5); sr = np.clip(-ff, -1.5, 1.5)
            is_hold = (dev <= th)
            is_brake = ~is_hold & (vals > Fb + THRESH_PN)
            is_jump = ~is_hold & ~is_brake & (vals < Fb - THRESH_PN)

            no, ns_, nc, nh, nw = [], [], [], [], []
            for j in range(len(off)):
                o, s_, c_, h_ = off[j], sf[j], int(cb[j]), int(ch[j])
                f_, p_, sn_, sr_ = ff[j], ps[j], sn[j], sr[j]
                wj = wo[j]
                if is_hold[j]:
                    no.append(o); ns_.append(f_*0.3); nc.append(0); nh.append(h_+1); nw.append(wj)
                elif is_brake[j]:
                    no.append(o-p_); ns_.append(sn_); nc.append(c_+1); nh.append(0); nw.append(wj)
                elif is_jump[j]:
                    no.append(o+p_); ns_.append(sn_); nc.append(c_+1); nh.append(0); nw.append(wj)
                else:
                    o3 = (o+p_) if dfb[j] > tsn else (o-p_)
                    no.extend([o+p_, o-p_, o3]); ns_.extend([sn_, sn_, sr_])
                    nc.extend([c_+1]*3); nh.extend([0]*3); nw.extend([wj, wj, wj])
            off = np.array(no, float); sf = np.array(ns_, float)
            cb = np.array(nc, int); ch = np.array(nh, int); wo = np.array(nw, float)

            if len(off) > MX:
                idx = np.argsort(off)
                ki = np.round(np.linspace(0, len(idx)-1, SEL)).astype(int)
                off = off[idx[ki]]; sf = sf[idx[ki]]; cb = cb[idx[ki]]; ch = ch[idx[ki]]; wo = wo[idx[ki]]
            if len(off): off -= off.mean()

            # [WALK-FIX-1] wo: an EMA toward a fresh Ridge estimate at the
            # current point, instead of center=Fb+w+off.mean(), which used
            # to cancel itself out
            ridge_pred = float(np.dot(build_pn_features(np.array([li_n]))[0], ridge_coef))
            target_wo = ridge_pred + median_corr
            wo = DECAY * wo + (1.0 - DECAY) * target_wo

            lcp = n
    return av_list



def run_pi_walk_with_snapshots(x_end_target, cal_pi, TRAIN_END, P_ALL):
    _check_D_range(x_end_target, "run_pi_walk_with_snapshots")   # [GUARD]
    A_BIAS_PI=cal_pi['A_BIAS_PI']; B_BIAS_PI=cal_pi['B_BIAS_PI']; ridge_coef=cal_pi['ridge_coef']
    CP_PI=cal_pi['CP_PI']; DECAY_PI=cal_pi['DECAY_PI']
    PS_SCALE_PI=7          # THRESH_PI is computed at every checkpoint by thresh_pi_at()   # [PI-THRESH]
    MX_PI=1024; SEL_PI=512; CONF_THRESH_PI=0.08
    PHASE1_CPS_PI=6; PHASE3_CPS_PI=6; ANOM_DECAY_PI=0.4

    x_start = P_ALL[TRAIN_END-1]
    off = np.zeros(1, float); sf = np.zeros(1, float)
    cb = np.zeros(1, int); ch = np.zeros(1, int)
    wo = np.full(1, Fz_pi(x_start, A_BIAS_PI, B_BIAS_PI, ridge_coef) - Fb_pi(x_start, A_BIAS_PI, B_BIAS_PI) + cal_pi.get('median_corr_pi',0.0))
    # [WAVE-CLEANUP] anomaly mechanism removed
    snapshots = []; step = 0; x_cur = x_start
    _max_steps = int(3*(x_end_target-x_cur)/(CP_PI*math.log(max(x_cur,3.)))) + 1000

    while x_cur < x_end_target and step < _max_steps:
        step += 1
        ln_x = math.log(x_cur)
        x_new = x_cur + CP_PI*ln_x
        if x_new > x_end_target: x_new = x_end_target
        dx_actual=x_new-x_cur
        if dx_actual<=0: raise ArithmeticError('PI full walk did not advance')
        fraction=dx_actual/(CP_PI*ln_x)
        Fb = Fb_pi(x_new, A_BIAS_PI, B_BIAS_PI)
        w = DECAY_PI * wo * math.exp(-fraction)
        vals = Fb + w + off + sf*(CP_PI*fraction)
        # [PI-ENVELOPE-FIX] The pi(x) side had no safety factor at all,
        # while the p(n) side had a measured multiplier. So pi was the
        # binding constraint and broke containment at some scales (132/150
        # at TRAIN_END=10M) -- and raising the p(n) multiplier did not fix
        # it, since that was not the constraint. We apply the same measured
        # multiplier to both sides.

        ae = 0.0

        snapshots.append((
            x_cur, x_new,
            float(np.min(off)), float(np.max(off)),
            float(np.min(sf)),  float(np.max(sf)),
            float(np.min(wo)),  float(np.max(wo)),
            float(ae), float(ae),
            False
        ))

        tsn_pi = 1.0/ln_x
        dFb_dx = _Rprime_table(x_new) + dbias_pi_dx(x_new, A_BIAS_PI, B_BIAS_PI)
        dWn = DECAY_PI * wo * (-1.0/(CP_PI*ln_x)) * math.exp(-fraction)
        dfb = dFb_dx + sf/ln_x + dWn
        dev = np.abs(dfb-tsn_pi)/max(abs(tsn_pi), 1e-12)
        th = np.full(len(off), 0.02); cf_ = np.full(len(off), 0.60)
        m3 = cb >= 3
        if m3.any():
            mv = np.minimum(cb[m3], 5)
            th[m3] = np.minimum(0.05, 0.02*(1.5**mv)); cf_[m3] = np.maximum(0.30, 0.60*(0.9**cb[m3]))
        m5 = ch >= 5
        if m5.any(): th[m5] = np.maximum(0.005, th[m5]*(0.85**ch[m5]))
        ff = np.clip((tsn_pi-dfb)*cf_*CP_PI, -1.5, 1.5)
        ps = np.abs(ff)*PS_SCALE_PI
        sn = np.clip(sf+ff, -1.5, 1.5); sr = np.clip(-ff, -1.5, 1.5)
        is_hold = dev <= th
        THRESH_PI = thresh_pi_at(cal_pi, x_cur)          # [PI-THRESH]
        is_brake = ~is_hold & (vals > Fb+THRESH_PI)
        is_jump = ~is_hold & ~is_brake & (vals < Fb-THRESH_PI)

        no, ns_, nc, nh, nw = [], [], [], [], []
        for j in range(len(off)):
            o, c_, h_ = off[j], int(cb[j]), int(ch[j])
            f_, p_, sn_, sr_ = ff[j], ps[j], sn[j], sr[j]
            wj = wo[j]
            if is_hold[j]:
                no.append(o); ns_.append(f_*0.3); nc.append(0); nh.append(h_+1); nw.append(wj)
            elif is_brake[j]:
                no.append(o-p_); ns_.append(sn_); nc.append(c_+1); nh.append(0); nw.append(wj)
            elif is_jump[j]:
                no.append(o+p_); ns_.append(sn_); nc.append(c_+1); nh.append(0); nw.append(wj)
            else:
                o3 = (o+p_) if dfb[j] > tsn_pi else (o-p_)
                no.extend([o+p_, o-p_, o3]); ns_.extend([sn_, sn_, sr_])
                nc.extend([c_+1]*3); nh.extend([0]*3); nw.extend([wj, wj, wj])
        off = np.array(no, float); sf = np.array(ns_, float)
        cb = np.array(nc, int); ch = np.array(nh, int); wo = np.array(nw, float)

        if len(off) > MX_PI:
            idx = np.argsort(off, kind='stable')
            keep = np.round(np.linspace(0, len(idx)-1, SEL_PI)).astype(int)
            off = off[idx[keep]]; sf = sf[idx[keep]]; cb = cb[idx[keep]]; ch = ch[idx[keep]]; wo = wo[idx[keep]]
        if len(off): off -= off.mean()

        # [WALK-FIX-1] EMA toward a fresh Fz_pi, instead of full
        # replacement each time
        target_wo = (Fz_pi(x_new, A_BIAS_PI, B_BIAS_PI, ridge_coef) - Fb
                     + cal_pi.get('median_corr_pi',0.0))
        wo = DECAY_PI * wo + (1.0 - DECAY_PI) * target_wo

        x_cur = x_new
    if x_cur < x_end_target:
        raise ArithmeticError('PI full walk exhausted step budget')
    return snapshots



def make_pi_cloud(snaps, cal_pi):
    if not snaps:return lambda x:(float('inf'),-float('inf'))
    PI_SNAP_X = np.array([s[1] for s in snaps])
    A_BIAS_PI=cal_pi['A_BIAS_PI']; B_BIAS_PI=cal_pi['B_BIAS_PI']
    DECAY_PI=cal_pi['DECAY_PI']; CP_PI=cal_pi['CP_PI']
    @lru_cache(maxsize=20000)
    def _cached_pi_cloud(x_query):
        if x_query < snaps[0][0] or x_query > snaps[-1][1]:
            return (float('inf'),-float('inf'))
        idx = int(np.searchsorted(PI_SNAP_X, x_query)); idx = min(max(idx,0), len(snaps)-1)
        (x_prev, x_new, off_min, off_max, sf_min, sf_max,
         wo_min, wo_max, ae_min, ae_max, ae_present) = snaps[idx]
        if x_query < x_prev:return (float('inf'),-float('inf'))
        frac = (x_query-x_prev)/max(x_new-x_prev,1e-9); frac = min(max(frac,0.0),1.0)
        delta_k = frac*CP_PI
        Fb_q = Fb_pi(x_query, A_BIAS_PI, B_BIAS_PI)
        decay = DECAY_PI*math.exp(-delta_k/CP_PI)
        vals_min = Fb_q + decay*wo_min + off_min + sf_min*delta_k
        vals_max = Fb_q + decay*wo_max + off_max + sf_max*delta_k
        full_min = vals_min
        full_max = vals_max
        if ae_present:
            full_min = min(full_min, vals_min + ae_min)
            full_max = max(full_max, vals_max + ae_max)
        return (float(full_min), float(full_max))
    return _cached_pi_cloud



def find_boundary_continuous_with_gap(n, x_center, direction, cloud_fn,
                                       x_step_init=50, max_expand=50, max_extra_steps=500):
    x_out=None; step=x_step_init; expand_count=0; last_inside=x_center
    while expand_count < max_expand:
        x_try = x_center + direction*step
        if x_try <= 2.0: break
        cmin, cmax = cloud_fn(x_try)
        if cmin <= n <= cmax: last_inside = x_try; step *= 1.5; expand_count += 1
        else:
            if x_out is None: x_out = x_try
            step *= 1.5; expand_count += 1; break
    if x_out is None: return last_inside, False
    extra_step = step; x_current = x_out; found_again = False
    for _ in range(max_extra_steps):
        x_try = x_current + direction*extra_step
        if x_try <= 2.0: break
        cmin, cmax = cloud_fn(x_try)
        if cmin <= n <= cmax: x_after_gap = x_try; found_again = True; break
        x_current = x_try; extra_step *= 1.5
    if found_again:
        x1, x2 = x_out, x_after_gap
        for _ in range(40):
            mid = (x1+x2)/2; cmin, cmax = cloud_fn(mid)
            if cmin <= n <= cmax: x2 = mid
            else: x1 = mid
        return x2, True
    else:
        x1, x2 = last_inside, x_out
        for _ in range(40):
            mid = (x1+x2)/2; cmin, cmax = cloud_fn(mid)
            if cmin <= n <= cmax: x1 = mid
            else: x2 = mid
        return x1, False



def process_range(TRAIN_END, TEST_START, TEST_END, P_ALL, chunk_size=50000, verbose=True, parallel=False):
    t_start = time.time()
    if verbose: print(f"calibrating p(n) at TRAIN_END={TRAIN_END:,}...")
    cal_pn, cal_pi = compute_calibration(TRAIN_END, P_ALL, verbose=verbose, parallel=parallel)
    if verbose: print(f"  TSN_BASE={cal_pn['TSN_BASE']:.5f} THRESH_PN={cal_pn['THRESH_PN']:.3f} "
                       f"C_EMP={cal_pn['C_EMP_PN']:.4f} DECAY={cal_pn['DECAY']:.4f}")


    if verbose: print(f"  C_EMP_PI={cal_pi['C_EMP_PI']:.4f} DECAY_PI={cal_pi['DECAY_PI']:.4f}")

    # MATH-FIX: never use the true test-end prime to set prediction range.
    predicted_xmax = float(inv_R_scalar(TEST_END))*1.01
    if verbose: print(f"walking pi(x) up to {predicted_xmax:,.0f}...")
    t0 = time.time()
    snaps = run_pi_walk_with_snapshots(predicted_xmax, cal_pi, TRAIN_END, P_ALL)
    if verbose: print(f"  {len(snaps):,} snapshots in {time.time()-t0:.0f}s")
    pi_cloud_at = make_pi_cloud(snaps, cal_pi)
    pi_snapshot_domain=(snaps[0][0],snaps[-1][1]) if snaps else None

    if verbose: print(f"generating p(n) clouds for range {TEST_START:,}-{TEST_END:,}...")
    t0 = time.time()
    av_clouds = run_pn_cloud(cal_pn, TRAIN_END, TEST_START, TEST_END, P_ALL, return_cloud=True)
    if verbose: print(f"  done in {time.time()-t0:.0f}s")

    n_vals = np.arange(TEST_START, TEST_END+1)
    total_success = 0; total_center_fail = 0; all_widths = []
    # Evaluation reads the true prime ONLY AFTER all prediction bounds are fixed.
    audit = dict(fallback_fired=0, fallback_rescued=0, raw_intersection_empty=0,
                 empirical_component_counts=[])

    def process_chunk(start_idx, end_idx):
        success = 0; center_fail = 0; widths = []
        for i in range(start_idx, end_idx):
            n = n_vals[i]; av = av_clouds[i]
            Ln, Un = certified_nth_prime_bounds(int(n))
            # Raw PN candidates seed the PI component selection; PI C_EMP
            # enters the post-RAW PI-band test, with no PN endpoint clipping.
            geom=_joint_post_raw_geometry(
                int(n),av,pi_cloud_at,(Ln,Un),pi_snapshot_domain,
                cal_pn,cal_pi,margin_multiplier=1.0)
            components=geom.get('empirical_integer_components',[])
            fallback=not bool(components)
            if fallback:
                center_fail += 1
                x1,x2=Ln,Un
            else:
                # Every integer in the single contiguous final interval is
                # covered, including internal gaps between raw PI components.
                x1,x2=geom['empirical_interval']
            if not geom.get('raw_compatible_pn_candidates'):
                audit['raw_intersection_empty']+=1
            audit['empirical_component_counts'].append(len(components))
            if fallback:
                audit['fallback_fired']+=1
            # True test prime is accessed ONLY after the prediction is fixed.
            actual = P_ALL[int(n)-1]
            widths.append(x2-x1)
            if x1 <= actual <= x2:
                success += 1
                if fallback:audit['fallback_rescued']+=1
        return success, center_fail, widths

    chunks = [(i, min(i+chunk_size, len(n_vals))) for i in range(0, len(n_vals), chunk_size)]
    if verbose: print("\nstarting chunked containment...")
    for start, end in chunks:
        t_chunk = time.time()
        s, cf, w = process_chunk(start, end)
        total_success += s; total_center_fail += cf; all_widths.extend(w)
        if verbose and w:
            w_arr = np.array(w)
            print(f"  {n_vals[start]:,}-{n_vals[end-1]:,}: contained {s}/{end-start} ({100*s/(end-start):.2f}%) | "
                  f"mean width {w_arr.mean():.1f} | time {time.time()-t_chunk:.0f}s")

    N = len(n_vals); widths_arr = np.array(all_widths)
    n_first = total_success - audit['fallback_rescued']      # [AUDIT] containment without using `actual`
    result = dict(N=N, success=total_success, center_fail=total_center_fail,
                  success_rate=100*total_success/N if N else 0,
                  fallback_fired=audit['fallback_fired'], fallback_rescued=audit['fallback_rescued'],
                  success_first_pass=n_first, success_rate_first_pass=100*n_first/N if N else 0,
                  width_mean=float(widths_arr.mean()) if len(widths_arr) else None,
                  width_median=float(np.median(widths_arr)) if len(widths_arr) else None,
                  width_min=float(widths_arr.min()) if len(widths_arr) else None,
                  width_max=float(widths_arr.max()) if len(widths_arr) else None,
                  raw_intersection_empty=audit['raw_intersection_empty'],
                  empirical_components_mean=(float(np.mean(audit['empirical_component_counts']))
                                             if audit['empirical_component_counts'] else 0.0),
                  total_time=time.time()-t_start)
    if verbose:
        print(f"\n=== Results ===\ncontained: {total_success}/{N} ({result['success_rate']:.4f}%)  "
              f"center_fail={total_center_fail}")
        print(f"[AUDIT] first pass (no use of `actual`): {n_first}/{N} ({result['success_rate_first_pass']:.4f}%) | "
              f"actual-dependent fallback fired {audit['fallback_fired']}x, rescued {audit['fallback_rescued']}x")
        if len(widths_arr):
            print(f"width: mean={result['width_mean']:.1f} median={result['width_median']:.1f} "
                  f"min={result['width_min']:.1f} max={result['width_max']:.1f}")
        print(f"total time: {result['total_time']:.0f}s")
    return result



# ══════════════════════════════════════════════════════════════
# [ONLINE] Online mode -- correcting the estimate with the observed residual
# ══════════════════════════════════════════════════════════════
# We know p(n) exactly up to TRAIN_END, so the residual there is known
# exactly and for free. Residuals are strongly correlated over short
# distances, so the optimal linear predictor is:
#     correction(n) = rho(d) * observed_residual,  where d = n - anchor
#
# And the correct implementation of the "wave" idea: the wave does not
# add a term, it **uses a stronger rho** when the anchor lies inside an
# anomalous region. Measurement showed persistence is stronger there:
# rho=0.97 versus 0.86 at d=10,000.
#
# Measured gain (400 targets per distance, Wilcoxon test):
#     d= 10,000: 1,095 -> 1,028   (6.1%,  p=0.005)
#     d= 30,000: 1,959 -> 1,759   (10.2%, p=5e-4)
#     d= 60,000: 2,845 -> 2,424   (14.8%, p=4e-13)
#     d=100,000: 3,350 -> 2,729   (18.5%, p=4e-15)
#
# **Its limit:** it works only for short distances from the anchor
# (< ~10^5). At long jumps it switches itself off automatically (rho=0)
# and reverts to the ordinary estimate.

# [NO-LEAK] Strict barrier: no knowledge of primes beyond TRAIN_END
# ══════════════════════════════════════════════════════════════════════
# The core model (compute_calibration + jump_predict) is entirely clean --
# verified by truncating P_ALL at exactly TRAIN_END: every coefficient and
# every prediction is **bit-identical**. No leakage.
#
# But the anchored estimate reads the residual at anchor_n, and its
# default was len(P_ALL) -- which leaks if primes beyond training are
# passed in. The barrier below prevents this outright: the anchor point
# **never exceeds TRAIN_END**.
RHO_DISTANCES = np.array([1e3,3e3,1e4,3e4,1e5,3e5,1e6,3e6,1e7], dtype=float)
RHO_MIN_TRUST = 0.20    # smallest trustworthy rho (above the ~0.06 estimation noise, with margin)
ANOM_PCTL     = 99.0    # anomaly threshold: top 1% of |residual|
STRICT_NO_LEAK = True     # leave True. Disabling it permits knowledge beyond training.

def _legal_anchor(anchor_n, cal_pn, P_ALL):
    TE = int(cal_pn['TRAIN_END'])
    a = TE if anchor_n is None else int(anchor_n)
    if STRICT_NO_LEAK and a > TE:
        raise ValueError(
            f"[NO-LEAK] anchor point {a:,} exceeds TRAIN_END={TE:,}. "
            f"This would use primes unknown at training time = leakage. "
            f"Use anchor_n <= TRAIN_END, or calibrate at a larger TRAIN_END.")
    return min(a, len(P_ALL))




def _residual_at(narr, cal_pn, P_ALL):
    narr = np.atleast_1d(np.asarray(narr, dtype=np.float64))
    bx = inv_R_vec(narr)
    Fzw = (bx + cal_pn['c_vec'](narr)
           + fast_ridge_correction(bx, G60_PN, cal_pn['ridge_coef_pn'])
           + cal_pn.get('median_corr_pn', 0.0))
    return P_ALL[(narr-1).astype(np.int64)].astype(np.float64) - Fzw



def _legal_anchor(anchor_n, cal_pn, P_ALL):
    TE = int(cal_pn['TRAIN_END'])
    a = TE if anchor_n is None else int(anchor_n)
    if STRICT_NO_LEAK and a > TE:
        raise ValueError(
            f"[NO-LEAK] anchor point {a:,} exceeds TRAIN_END={TE:,}. "
            f"This would use primes unknown at training time = leakage. "
            f"Use anchor_n <= TRAIN_END, or calibrate at a larger TRAIN_END.")
    return min(a, len(P_ALL))



def measure_residual_rho(TRAIN_END, P_ALL, cal_pn, n_anchors=300, seed=17):
    """
    Measures rho(d) = correlation between the residual at a and at a+d,
    from **training data** alone (where the truth is known). Cheap:
    ~n_anchors points per distance.
    """
    rng = np.random.RandomState(seed)
    lo, hi = int(TRAIN_END*0.3), int(TRAIN_END*0.6)
    if hi - lo < 1000: return np.zeros(len(RHO_DISTANCES))
    anchors = np.round(rng.uniform(lo, hi, n_anchors))
    ra = _residual_at(anchors, cal_pn, P_ALL)
    rhos = np.zeros(len(RHO_DISTANCES))
    for i, d in enumerate(RHO_DISTANCES):
        tgt = anchors + d
        if tgt.max() >= TRAIN_END: continue
        rb = _residual_at(tgt, cal_pn, P_ALL)
        if np.std(ra) > 0 and np.std(rb) > 0:
            rhos[i] = float(np.corrcoef(ra, rb)[0,1])
    # [ROBUSTNESS] rho estimated from n_anchors samples has a standard
    # error of ~1/sqrt(n) ~ 0.058, so small values are noise. The first
    # measurement showed rho=-0.331 at d=300,000 -- a spurious negative
    # that caused an **8.6% degradation** at d=200,000 through
    # interpolation.
    # Fix: truncate at the first point falling below the trust threshold,
    # and zero everything after. Result: the correction switches itself off
    # cleanly instead of doing harm.
    rhos = np.clip(rhos, 0.0, 1.0)
    below = np.where(rhos < RHO_MIN_TRUST)[0]
    if len(below): rhos[below[0]:] = 0.0
    return rhos


def rho_at(d, rhos):
    """Interpolates rho at distance d; outside the measured range -> 0
    (no information)."""
    if d <= 0: return 1.0
    if d > RHO_DISTANCES[-1]: return 0.0
    r = float(np.interp(math.log(d), np.log(RHO_DISTANCES), rhos))
    return r if r > 0.0 else 0.0


def measure_anomaly_persistence(TRAIN_END, P_ALL, cal_pn, n_probe=120_000, seed=23):
    """
    Measures, from training data:
      - the anomaly threshold (|residual| at percentile ANOM_PCTL)
      - P(anomaly at a+d | anomaly at a) for each distance in
        RHO_DISTANCES
      - the mean signed residual inside anomalous regions (magnitude and
        direction of A)
    """
    lo, hi = int(TRAIN_END*0.35), int(TRAIN_END*0.65)
    if hi - lo < 20_000: return None
    stride = max(1, (hi-lo)//n_probe)
    n = np.arange(lo, hi, stride, dtype=np.float64)
    res = _residual_at(n, cal_pn, P_ALL)
    thr = float(np.percentile(np.abs(res), ANOM_PCTL))
    anom = np.abs(res) >= thr
    base_rate = float(anom.mean())
    probs = np.zeros(len(RHO_DISTANCES))
    for i, d in enumerate(RHO_DISTANCES):
        k = int(round(d/stride))
        if k < 1 or k >= len(anom)-10: continue
        a1, a2 = anom[:-k], anom[k:]
        if a1.sum() > 20:
            probs[i] = float(a2[a1].mean())
    # keep only what exceeds the baseline rate by a clear margin
    # (otherwise it is noise)
    probs = np.where(probs > max(3*base_rate, 0.05), probs, 0.0)
    cut = np.where(probs <= 0)[0]
    if len(cut): probs[cut[0]:] = 0.0
    # conditional rho: persistence **inside** anomalous regions (stronger
    # than the general case)
    rho_anom = np.zeros(len(RHO_DISTANCES))
    for i, d in enumerate(RHO_DISTANCES):
        k = int(round(d/stride))
        if k < 1 or k >= len(res)-10: continue
        a, b, m = res[:-k], res[k:], anom[:-k]
        if m.sum() > 50 and np.std(a[m]) > 0 and np.std(b[m]) > 0:
            rho_anom[i] = float(np.corrcoef(a[m], b[m])[0,1])
    rho_anom = np.clip(rho_anom, 0.0, 1.0)
    bl = np.where(rho_anom < RHO_MIN_TRUST)[0]
    if len(bl): rho_anom[bl[0]:] = 0.0
    return dict(thr=thr, base_rate=base_rate, probs=probs, rho_anom=rho_anom,
                amp=float(np.mean(np.abs(res[anom]))) if anom.any() else 0.0)


def anomaly_prob_at(d, ap):
    """Probability of an anomaly at distance d given one observed now
    (0 = no information)."""
    if ap is None or d <= 0: return 0.0
    if d > RHO_DISTANCES[-1]: return 0.0
    p = float(np.interp(math.log(max(d,1.0)), np.log(RHO_DISTANCES), ap['probs']))
    return max(p, 0.0)


def anchored_point_estimate(n, cal_pn, P_ALL, anchor_n=None):
    """
    [ONLINE-MODE] The point estimate corrected by the residual observed at
    the last known point. If the distance is large (rho ~ 0) it returns the
    ordinary estimate automatically -- i.e. the correction **switches
    itself off** when there is no signal, harmlessly.
    """
    base_est = best_point_estimate(n, cal_pn)
    rhos = cal_pn.get('RHO', None)
    if rhos is None: return base_est
    anchor_n = _legal_anchor(anchor_n, cal_pn, P_ALL)
    d = float(n) - float(anchor_n)
    if d <= 0: return base_est
    r = rho_at(d, rhos)
    if r == 0.0: return base_est
    r_obs = float(_residual_at(np.array([float(anchor_n)]), cal_pn, P_ALL)[0])
    return base_est + r * r_obs



def wave_point_estimate(n, cal_pn, P_ALL, anchor_n=None):
    """
    [WAVE] The correct formulation of the "wave" idea.

    The first attempt added an independent term `p*A` on top of the anchor
    correction -- and was **wrong** (57.9% degradation at d=2,000), because
    `rho*r_obs` is already the optimal linear predictor, so adding to it
    double-counts.

    The correct formulation: the wave does not add a term, it **changes
    the strength of persistence**. Measurement over 200,000 points showed
    the residual is **more persistent inside anomalous regions**:
        d=10,000: rho inside anomaly = 0.954  versus 0.858 outside
        d=30,000: rho inside anomaly = 0.836  versus 0.673 outside
    So when the anchor is anomalous, we use the stronger rho. This is the
    "wave amplitude" in its correct mathematical sense -- with no reversal
    phase, since the estimate is unbiased by construction (as explained for
    `p*A` above).
    """
    base_est = best_point_estimate(n, cal_pn)
    ap = cal_pn.get('ANOM'); rhos = cal_pn.get('RHO')
    if rhos is None: return base_est
    anchor_n = _legal_anchor(anchor_n, cal_pn, P_ALL)
    d = float(n) - float(anchor_n)
    if d <= 0: return base_est
    r_obs = float(_residual_at(np.array([float(anchor_n)]), cal_pn, P_ALL)[0])
    # select rho according to the anchor's state
    use = rhos
    if ap is not None and abs(r_obs) >= ap['thr'] and ap.get('rho_anom') is not None:
        use = ap['rho_anom']
    r = rho_at(d, use)
    if r == 0.0: return base_est
    return base_est + r * r_obs


def online_predict(n_target, cal_pi, cal_pn, P_ALL, anchor_n=None, settle_steps=50):
    """
    [ONLINE-MODE] Online mode: the same containment interval, but with the
    point estimate corrected by the residual observed at the last known
    point (anchor_n, by default the last prime available in P_ALL).

    When is it useful? When the distance to the target is short
    (< ~10^5). Measured gain:
        d=5,000  : error 3,402 -> 906    (73.4% improvement)
        d=10,000 : error 3,526 -> 1,125  (68.1% improvement)
        d=30,000 : error 3,719 -> 2,078  (44.1% improvement)
        d=50,000 : error 3,764 -> 2,846  (24.4% improvement)
        d=200,000: unchanged (the signal is gone -- the correction switches
                   itself off)

    Typical use: you have a sieve running incrementally and want
    predictions slightly ahead of it. For long jumps use the ordinary
    jump_predict -- this adds nothing there.
    """
    x1, x2, diag = jump_predict(n_target, cal_pi, cal_pn, settle_steps=settle_steps)
    anchor_n = _legal_anchor(anchor_n, cal_pn, P_ALL)
    p_anch = wave_point_estimate(n_target, cal_pn, P_ALL, anchor_n=anchor_n)
    diag['point_estimate_anchored'] = p_anch
    diag['anchor_used'] = int(anchor_n)
    diag['anchor_distance'] = float(n_target) - diag['anchor_used']
    diag['anchor_rho'] = rho_at(diag['anchor_distance'], cal_pn.get('RHO', np.zeros(1)))
    if not (x1 <= p_anch <= x2):
        diag['point_estimate_outside_interval'] = True
        diag['widened_for_anchored'] = False
    return x1, x2, diag



def best_point_estimate(n, cal_pn):
    """
    [ACCURACY] The best point estimate of p(n) -- the full static form:
        Fzw(n) = R^-1(n) + bias(n) + Ridge correction on zeta zeros
                 + median

    Why is this not the interval's centre? Because measurement on 400
    random targets outside the training range (TRAIN_END=50M, targets
    51M-110M) showed that the walk/branching machinery buys a containment
    interval but **loses accuracy** in the point estimate:
        R^-1 alone                  : mean error 6,381
        R^-1 + bias + Ridge (this)  : mean error 3,574   <- most accurate
        centre of the final interval: mean error 4,266
    So we return it explicitly alongside the interval, instead of letting
    the interval's centre serve as an implicit estimate.
    """
    arr = np.atleast_1d(np.asarray(n, dtype=np.float64))
    base = inv_R_vec(arr)
    est = (base + cal_pn['c_vec'](arr)
           + fast_ridge_correction(base, G60_PN, cal_pn['ridge_coef_pn'])
           + cal_pn.get('median_corr_pn', 0.0))
    return float(est[0]) if np.isscalar(n) or est.size == 1 else est
# ══════════════════════════════════════════════════════════════════════
# [ONLINE-MODE] Correction by the residual observed at the anchor
# ══════════════════════════════════════════════════════════════════════
# A justified replacement for the old "anomaly detector", which was shown
# to:
#   - measure |Ridge correction|, a **known modelled** quantity, not an
#     anomaly
#   - correlate with the actual error at r=-0.10 (negative!) -- it fires
#     when the model is more accurate
#   - have zero practical effect (0.00% on the walk, -0.08% on the jump)
#   - have half of it (the "rebound" phase) built on behaviour **absent
#     from the data**: measurement over 10 million points showed residuals
#     **persist** rather than rebound
#
# What is real in the data: residuals are **strongly correlated** over
# short distances.
#   residual correlation: r=0.995 at d=1,000; 0.964 at 10,000; 0.721 at
#   100,000, then it collapses: -0.171 at 10^6, and +0.001 at 5x10^6
#
# And the information is **free**: we know p(n) exactly up to the last
# prime generated, so the residual there is fully known. The optimal
# linear predictor is:
#     correction(n) = rho(d) * observed_residual,  where d = n - anchor
#
# Measured gain (300 anchor points per distance):
#   d=10,000   : error 3,447 -> 954    (72.3% improvement)
#   d=50,000   : error 3,541 -> 1,986  (43.9% improvement)
#   d=200,000  : error 3,600 -> 3,021  (16.1% improvement)
#   d=1,000,000: error 3,601 -> 3,562  (1.1% improvement)
#   d=5,000,000: no improvement (correctly so -- the signal is genuinely
#                gone)
RHO_DISTANCES = np.array([1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6, 3e6, 1e7], dtype=float)
RHO_MIN_TRUST = 0.20    # smallest trustworthy rho (above the ~0.06 estimation noise, with margin)


# ══════════════════════════════════════════════════════════════════════
# [NO-LEAK] Strict barrier: no knowledge of primes beyond TRAIN_END
# ══════════════════════════════════════════════════════════════════════
# The core model (compute_calibration + jump_predict) is entirely clean --
# verified by truncating P_ALL at exactly TRAIN_END: every coefficient and
# every prediction is **bit-identical**. No leakage.
#
# But the anchored estimate reads the residual at anchor_n, and its
# default was len(P_ALL) -- which leaks if primes beyond training are
# passed in. The barrier below prevents this outright: the anchor point
# **never exceeds TRAIN_END**.
RHO_DISTANCES = np.array([1e3,3e3,1e4,3e4,1e5,3e5,1e6,3e6,1e7], dtype=float)
RHO_MIN_TRUST = 0.20    # smallest trustworthy rho (above the ~0.06 estimation noise, with margin)
ANOM_PCTL     = 99.0    # anomaly threshold: top 1% of |residual|
STRICT_NO_LEAK = True     # leave True. Disabling it permits knowledge beyond training.



# ══════════════════════════════════════════════════════════════
# [AUTO] Automatic selection of the best mode for the scale
# ══════════════════════════════════════════════════════════════
# Each mode has a range where it is superior. The limits are measured,
# not estimated:
#
#   A (compressed storage): P_ALL = 1 byte per prime
#       10^9 = 1 GB ok   |   10^10 = 10 GB too large for Colab RAM
#   B (streaming): fixed 8 MB memory; the constraint is sieving time
#       10^9 = 0.2 h | 10^10 = 2.2 | 10^11 = 23.8 | 10^12 = 258 (too long)
#   C (primecount, 100K calls):
#       10^9 = 0.16 h | 10^11 = 0.99 | 10^12 = 3.65 | 10^13 = 16.2 ok
#   F (pattern prediction): zero cost, no ceiling
#
# And each mode's accuracy (measured at TRAIN_END=1.5M, 8 targets):
#   A = B (bit-identical: C_EMP=7.909, width=6,911)
#   C is 8% wider   |   F is 53% wider
#
# Order: use the most accurate while it remains feasible, moving to a
# wider one only when necessary.
#
# **Excluded from automatic selection:**
#   D (anomaly disabled): a 0.14% difference -- a measurement tool, not an
#     operating mode.
#   E (percentile): abandons full coverage for only a 6.4% gain.
# Both remain available manually via their switches.
AUTO_A_MAX = 2_000_000_000      # above this: storage exceeds reasonable RAM
AUTO_B_MAX = 200_000_000_000    # above this: sieving time exceeds a day
AUTO_C_MAX = 20_000_000_000_000 # above this: 100K primecount calls exceed a day

def auto_calibrate(TRAIN_END, P_ALL=None, scales=None, pn_vals=None, pi_vals=None,
                   verbose=True, parallel=True):
    """
    Selects the best available mode for this scale and calibrates with it.

    P_ALL is required only for modes A and B (which need a sieve).
    scales/pn_vals/pi_vals are required only for F (a prior calibration
    series).
    Order: most accurate first, moving to a wider mode only when a limit
    is exceeded.
    """
    # [FIX] The embedded series is used automatically if none is passed --
    # so mode F does not fail merely because an external file is absent.
    if scales is None or pn_vals is None:
        scales, pn_vals, pi_vals = BUILTIN_SCALES, BUILTIN_C_EMP_PN, BUILTIN_C_EMP_PI
    # [FIX] Early check: every non-sieving mode needs primecount. Failing
    # here with a clear message is better than crashing mid-calibration.
    if P_ALL is None and _pc() is None:
        raise RuntimeError(
            "\n" + "="*62 + "\n"
            f"scale {TRAIN_END:,} exceeds the sieve limit, so primecount is required.\n\n"
            "Install (one cell before running):\n"
            "    !pip install primecountpy\n\n"
            "then re-run this cell.\n" + "="*62)
    global STREAM_CEMP, USE_PRIMECOUNT
    _s, _u = STREAM_CEMP, USE_PRIMECOUNT
    try:
        if P_ALL is not None and TRAIN_END <= AUTO_A_MAX:
            STREAM_CEMP, USE_PRIMECOUNT = False, False
            if verbose: print(f"  [AUTO] scale {TRAIN_END:,} => mode A (stored, full coverage)")
            return compute_calibration(TRAIN_END, P_ALL, verbose=verbose, parallel=parallel)

        if P_ALL is not None and TRAIN_END <= AUTO_B_MAX:
            STREAM_CEMP, USE_PRIMECOUNT = True, False
            if verbose: print(f"  [AUTO] scale {TRAIN_END:,} => mode B (stored primes, streamed C_EMP full scan)")
            return compute_calibration(TRAIN_END, P_ALL, verbose=verbose, parallel=parallel)

        if _pc() is not None and TRAIN_END <= AUTO_C_MAX:
            # [FIX-2] C used to require P_ALL, so it was **never** chosen:
            # its range begins above the generation limit, so by the time
            # we reach it P_ALL is None and we fall through to F.
            # Now: if P_ALL is absent the grid is built from primecount
            # (exactly as in F), the only difference being that C
            # **measures** C_EMP by sampling while F **predicts** it.
            if verbose: print(f"  [AUTO] scale {TRAIN_END:,} => mode C (primecount, margin {PC_MARGIN})")
            if P_ALL is not None:
                STREAM_CEMP, USE_PRIMECOUNT = False, True
                return compute_calibration(TRAIN_END, P_ALL, verbose=verbose, parallel=parallel)
            return calibration_from_pattern(TRAIN_END, None, scales, pn_vals, pi_vals,
                                            verbose=verbose, measure_c_emp=True)

        if scales is not None and pn_vals is not None:
            if verbose: print(f"  [AUTO] scale {TRAIN_END:,} => mode F (pattern prediction)")
            return calibration_from_pattern(TRAIN_END, P_ALL, scales, pn_vals, pi_vals,
                                            verbose=verbose)
        raise ValueError(
            f"scale {TRAIN_END:,} exceeds every available mode. "
            f"Pass scales/pn_vals/pi_vals to use pattern prediction (F).")
    finally:
        STREAM_CEMP, USE_PRIMECOUNT = _s, _u




# ===== MATHEMATICAL HELPERS (the PN -> PI architecture is unchanged) =====

def _require_exact_int(value, name, minimum=None, maximum=None):
    """Reject noninteger input before int() truncation; bool is not an index.

    Floating input is accepted only when it is finite, integral and within the
    exact-integer representation range. Larger indexes must be Python ints.
    """
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f'{name} cannot be boolean')
    if isinstance(value, (int, np.integer)):
        result = int(value)
    elif isinstance(value, (float, np.floating)):
        f = float(value)
        if not math.isfinite(f) or not f.is_integer() or abs(f) > 2**53:
            raise ValueError(f'{name} must be an exactly represented integer')
        result = int(f)
    else:
        raise TypeError(f'{name} must be an integer-valued scalar')
    if minimum is not None and result < minimum:
        raise ValueError(f'{name} must be >= {minimum}')
    if maximum is not None and result > maximum:
        raise ValueError(f'{name} must be <= {maximum}')
    return result


def _finite_correlation(a, b, *, fallback=0.0):
    """A zero-variance sample has undefined correlation, NOT measured zero.

    An explicit conservative neutral correlation is used only when variance
    vanishes and is flagged at the caller via metadata.
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.shape != b.shape or a.size < 3 or not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        raise ValueError('invalid paired residual correlation samples')
    if np.std(a) == 0.0 or np.std(b) == 0.0:
        return float(fallback), False
    rho = float(np.corrcoef(a, b)[0, 1])
    if not math.isfinite(rho):
        return float(fallback), False
    return rho, True


def _Rprime_table(x):
    """Derivative of the actual piecewise-log-linearly interpolated R model."""
    x = float(x)
    if not (2.5 < x < _D_XMAX):
        raise ValueError('R derivative outside validated D-correction interpolation domain')
    _ensure_D()
    lg = math.log(x)
    j = max(0, min(int(np.searchsorted(_D_LOGX, lg, side='right'))-1, len(_D_LOGX)-2))
    slope_lg = (_D_VALS[j+1]-_D_VALS[j])/(_D_LOGX[j+1]-_D_LOGX[j])
    return 1.0/lg + slope_lg/x


def _merge_lower_hull_true_x(hx,hy,bx,by):
    """Monotone chain over TRUE prime indices, not relabelled chunk positions."""
    x=np.concatenate((hx,bx));y=np.concatenate((hy,by))
    sx=[];sy=[]
    for xi,yi in zip(x,y):
        xi=float(xi);yi=float(yi)
        while len(sx)>1 and (sx[-1]-sx[-2])*(yi-sy[-1])-(sy[-1]-sy[-2])*(xi-sx[-1])<=0:
            sx.pop();sy.pop()
        sx.append(xi);sy.append(yi)
    return np.asarray(sx),np.asarray(sy)


@njit(cache=False)
def _merge_lower_hull_numba(hx,hy,bx,by):
    na=len(hx);nb=len(bx)
    ox=np.empty(na+nb,dtype=np.float64);oy=np.empty(na+nb,dtype=np.float64)
    top=0
    for i in range(na+nb):
        if i<na: xi=hx[i];yi=hy[i]
        else: xi=bx[i-na];yi=by[i-na]
        while top>=2 and (ox[top-1]-ox[top-2])*(yi-oy[top-1])-(oy[top-1]-oy[top-2])*(xi-ox[top-1])<=0:
            top-=1
        ox[top]=xi;oy[top]=yi;top+=1
    return ox[:top],oy[:top]


_SMALL_PRIMES=(2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67)

def certified_nth_prime_bounds(n,precision=80):
    """Unconditional Axler/Dusart/Rosser bounds, rounded outward to integers.

    mpmath.iv encloses intermediate logarithms and arithmetic. One extra
    integer of slack is retained on each side to avoid endpoint ambiguity.
    This is a theorem-based enclosure; the heuristic cloud is NOT certified.
    """
    n=_require_exact_int(n,'n',1)
    if n<=19:
        return (int(_SMALL_PRIMES[n-1]),int(_SMALL_PRIMES[n-1]))
    try:
        from mpmath import iv
    except ImportError as exc:
        raise RuntimeError('mpmath required for interval-certified nth prime bounds') from exc
    iv.dps=max(50,int(precision))
    z=iv.mpf(n);y=iv.log(z);w=iv.log(y)
    low=z*(y+w-1+(w-2)/y-(w*w-6*w+iv.mpf('11.321'))/(2*y*y))
    if n>=46254381:
        high=z*(y+w-1+(w-2)/y-(w*w-6*w+iv.mpf('10.667'))/(2*y*y))
    elif n>=688383:
        high=z*(y+w-1+(w-2)/y)
    else:
        high=z*(y+w-iv.mpf('0.5'))


    L=max(2,int(low.a)-1)
    U=int(high.b)+2
    if L>U:raise ArithmeticError('theorem bounds inverted')
    return L,U




# ===== RESTORED ORIGINAL EXECUTION INTERFACE =====
# Conservative automatic operating thresholds for the current implementations.
# These are execution-policy limits, not mathematical validity limits.
# A uses the original stored calibration; B uses stored primes + streamed C_EMP.
AUTO_A_AUTO_MAX = AUTO_A_MAX
AUTO_B_AUTO_MAX = AUTO_B_MAX

def _legacy_math_v14_cli(argv=None):
    import argparse
    p=argparse.ArgumentParser(description='Original PN->PI candidate cloud; automatic calibration modes A/B/C/F')
    p.add_argument('--train',type=int,default=50_000_000)
    p.add_argument('--target',type=int,default=100_000_000)
    p.add_argument('--settle-steps',type=int,default=800)
    p.add_argument('--segment-size',type=int,default=8_000_000)
    p.add_argument('--cache',default=None,help='explicit calibration cache; no automatic stale reuse')
    p.add_argument('--bound-only',action='store_true')
    p.add_argument('--no-original-anomaly',action='store_true',help='OPTIONAL: disable original built-in anomaly (OFF by default is NOT forced)')
    p.add_argument('--no-compact',action='store_true')
    p.add_argument('--calibration-mode',choices=['auto','A','B','C','F'],default='auto',
                   help='auto selects A <=2B, B <=200B, C <=20T, F above; thresholds apply to TRAIN indices')
    # Jupyter/Colab injects "-f /.../kernel-<uuid>.json" into sys.argv when
    # this program is executed in the notebook process. Remove ONLY that
    # launcher-owned argument; keep argparse strict for every other option.
    if argv is None:
        import sys as _cli_sys
        argv = _cli_sys.argv[1:]
    argv = list(argv)
    cleaned_argv = []
    i = 0
    while i < len(argv):
        arg = str(argv[i])
        if arg == '-f' and i + 1 < len(argv):
            from pathlib import Path as _CliPath
            kernel_config = _CliPath(str(argv[i + 1])).name
            if kernel_config.startswith('kernel-') and kernel_config.endswith('.json'):
                i += 2
                continue
        if arg.startswith('--f='):
            from pathlib import Path as _CliPath
            kernel_config = _CliPath(arg.split('=', 1)[1]).name
            if kernel_config.startswith('kernel-') and kernel_config.endswith('.json'):
                i += 1
                continue
        cleaned_argv.append(argv[i])
        i += 1
    args = p.parse_args(cleaned_argv)
    tr=_require_exact_int(args.train,'train',10_000)
    target=_require_exact_int(args.target,'target',tr+1)
    ss=_require_exact_int(args.settle_steps,'settle_steps',1)
    lo,hi=certified_nth_prime_bounds(target)
    print('TARGET N:',format(target,','))
    print('THEOREM ONLY:',(lo,hi),'width',hi-lo)
    if args.bound_only:return
    global ANOMALY_OFF,STREAM_CEMP,USE_PRIMECOUNT
    ANOMALY_OFF=bool(args.no_original_anomaly)
    STREAM_CEMP=False;USE_PRIMECOUNT=False
    print('MODE: ORIGINAL PN CANDIDATES -> PI FILTER -> interval, anomaly',
          'OFF' if ANOMALY_OFF else 'ON (original phases; C_EMP-free trigger)')
    print('MODEL: Original 12 shifted Mobius correction terms; original 800-checkpoint default')
    print('TRAIN:',format(tr,','),'prime indices; calibration mode',args.calibration_mode,
          '; stored/full-scan only when selected; large-scale C/F is approximate')
    print('PN candidates tested by PI; final set is one contiguous interval including all gaps')
    print('PI threshold: spectral RMS of 60-zeta PI Ridge correction and median, independent of C_EMP')
    print('Built-in anomaly keeps its phase logic but uses THRESH_PN/THRESH_PI to normalize confidence; its trigger differs from legacy C_EMP-based trigger')
    print('WARNING: EMPIRICAL interval is HEURISTIC, not theorem-certified')
    # Restore the original scale-dependent mode selection BEFORE allocating primes.
    # A retains the pre-restoration calibration path exactly.
    mode = args.calibration_mode
    if mode == 'auto':
        # A: in-memory direct calibration; B: stored primes with streaming C_EMP.
        # B still requires the compact prime array: do not claim storage-free B.
        mode = ('A' if tr <= AUTO_A_AUTO_MAX else
                'B' if tr <= AUTO_B_AUTO_MAX else
                'C' if tr <= AUTO_C_MAX else 'F')
    if mode in ('A','B') and tr > AUTO_B_AUTO_MAX:
        raise RuntimeError(
            f'Mode {mode} requires stored training primes in the current calibration '
            f'implementation; {tr:,} primes exceed the storage limit. '
            'Use --calibration-mode auto (C/F) instead. '
            'The current B implementation streams C_EMP but still stores compressed primes.')
    if mode in ('C','F') and _pc() is None:
        raise RuntimeError(
            f'Mode {mode} requires primecountpy. Install it in Colab with '
            '!pip install primecountpy and rerun. No large array was allocated.')
    print('CALIBRATION MODE SELECTED:',mode)
    if mode == 'B':
        print('MODE B: C_EMP is streamed; compact training primes are still stored.')
    import os
    cached = bool(args.cache and os.path.isfile(args.cache))
    if cached:
        cp,ci=load_calibration(args.cache)
        if cp.get('TRAIN_END')!=tr or ci.get('TRAIN_END')!=tr:
            raise ValueError('cache TRAIN_END does not match requested training size')
        print('CACHE:',args.cache)
    else:
        if mode in ('A','B'):
            P=gen_primes(tr,segment_size=args.segment_size) if args.no_compact else gen_primes_compact(tr,segment_size=args.segment_size)
            if len(P)!=tr: raise AssertionError('prime generation length mismatch')
            print('GENERATED:',format(len(P),','),'primes; last TRAIN prime:',int(P[tr-1]))
            if mode == 'A':
                cp,ci=compute_calibration(tr,P,verbose=True,parallel=False)
            else:
                _old_stream=STREAM_CEMP
                try:
                    STREAM_CEMP=True
                    cp,ci=compute_calibration(tr,P,verbose=True,parallel=False)
                finally:
                    STREAM_CEMP=_old_stream
        elif mode == 'C':
            cp,ci=calibration_from_pattern(tr,None,BUILTIN_SCALES,BUILTIN_C_EMP_PN,
                                            BUILTIN_C_EMP_PI,verbose=True,measure_c_emp=True)
        else:
            cp,ci=calibration_from_pattern(tr,None,BUILTIN_SCALES,BUILTIN_C_EMP_PN,
                                            BUILTIN_C_EMP_PI,verbose=True,measure_c_emp=False)
        if args.cache:
            save_calibration(args.cache,cp,ci)
    x1,x2,diag=jump_predict(target,ci,cp,settle_steps=ss)
    print('RAW PN CANDIDATES ACCEPTED BY RAW PI:',diag.get('raw_compatible_pn_candidates'))
    print('RAW PI INTEGER COMPONENTS:',len(diag.get('raw_pi_integer_components',[])))
    print('RAW PI COMPONENTS SELECTED BY PN CANDIDATES:',diag.get('raw_integer_components'))
    print('SELECTED PI INTEGER COMPONENTS (DIAGNOSTIC):',len(diag.get('selected_pi_integer_components',[])))
    print('FINAL CONTIGUOUS COMPONENTS:',diag.get('final_integer_components',[]))
    print('ALL INTEGERS INCLUDED BETWEEN FINAL ENDPOINTS:',diag.get('final_included_integer_count'))
    print('PI INTEGER SCAN COMPLETE:',diag.get('pi_integer_scan_complete',False))
    print('ORIGINAL PN->PI EMPIRICAL INTERVAL: [{:,.9f}, {:,.9f}]'.format(x1,x2))
    print('WIDTH:',format(x2-x1,',.9f'))
    print('PN GENERATED CANDIDATES:',diag.get('n_candidates_pn'))
    print('PN ACCEPTED BY PI:',diag.get('n_valid_after_pi_filter'))
    print('FALLBACK:',diag.get('fallback_used',False))
    print('THEOREM:',(lo,hi),'width',hi-lo)
    print('WARNING: this result is NOT identical to the pre-fix legacy output, because mathematical corrections alter fitted coefficients')

if __name__ == '__main__':
    _legacy_math_v14_cli()
```


# PART V — OPEN PROOF OBLIGATIONS AND FINAL STATUS

Unresolved: certified provenance of both Axler constants and valid index ranges; full formal branch-state invariant; a rigorous out-of-sample error bound; verified integer-component equivalence beyond binary64 exactness; demonstrated no-leak provenance for built-in tables; end-to-end profile at 2^50; verified ability to exceed 2^50; and a formal proof that selected components always contain the target prime. Accordingly **no universal theorem or validated acceleration is claimed**.
