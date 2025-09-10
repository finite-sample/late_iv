### LATE IV: Distributional Implications of LATE 

Many randomized encouragement designs have imperfect compliance, where only a fraction of people comply with their assignment. Examples include phone-bank get-out-the-vote (GOTV) campaigns and draft lotteries like the Vietnam Draft Lottery. In these settings, it is common to use instrumental variable (IV) regression for analysis.

Instrumental variables identify the **Local Average Treatment Effect (LATE)** under three standard conditions: random assignment, exclusion, and monotonicity. The exclusion restriction embodies a key intuition: we don't expect GOTV efforts over the phone to affect people we aren't able to reach, nor do we expect people who were merely drawn up in the Vietnam Draft Lottery to have different attitudes towards minorities solely from being drafted—whatever effects we see, we expect them to result from actual service.

This leads to a sharp testable implication. When the instrument $Z$ flips, only **compliers** can change treatment status. Therefore, only compliers can contribute to the treatment effect, which implies that the distribution of the intention-to-treat (ITT) effect follows a specific pattern—one that looks closer to Figure 1d [here (pdf)](http://www.stat.columbia.edu/~gelman/research/unpublished/causal_quartets.pdf), with a lumpy distribution reflecting the complier-only effects.

We turn this insight into testable implications about distributional equalities. One way to check if the data are consistent with these implications is to simulate the lumpy treatment effect pattern and check how closely the empirical distribution matches the theorized one. We provide estimators for the complier outcome distributions, uniform tests of these distributional implications, a complementary GMM test, and simulation evidence to validate the approach.

### Setup and consequence

Assume random assignment. Assume exclusion so that $Y=Y(D)$. Assume monotonicity so that $D(1)\ge D(0)$. Assume SUTVA.

Define the complier share

$$
p_C \;=\; \Pr\{D(1)>D(0)\} \;=\; \mathbb{E}[D\mid Z{=}1]-\mathbb{E}[D\mid Z{=}0].
$$

Let $F_{Y\mid Z=z}$ be the cumulative distribution function of $Y$ under $Z=z$. Let $F_{1C}(y)=\Pr\{Y(1)\le y\mid C\}$ and $F_{0C}(y)=\Pr\{Y(0)\le y\mid C\}$.

**Only‑compliers‑move identity**

$$
F_{Y\mid Z=1}(y)-F_{Y\mid Z=0}(y)\;=\;p_C\Big(F_{1C}(y)-F_{0C}(y)\Big)\quad\text{for all }y\in\mathbb{R}.
$$

This equality is the formal statement that noncompliers produce no distributional movement when $Z$ changes.

### Identification and estimation

With observed covariates $X$, define

$$
e(X)=\Pr(Z{=}1\mid X),\qquad p_z(X)=\mathbb{E}[D\mid Z{=}z,X],\qquad p_C=\mathbb{E}\big[p_1(X)-p_0(X)\big].
$$

Use Abadie‑style weights to identify complier marginals. For any measurable function $g$,

$$
\mathbb{E}[g(Y(1))\mid C] \;=\; \frac{\mathbb{E}\!\left[g(Y)\,\frac{Z}{e(X)}\big(D-p_0(X)\big)\right]}{p_C},
\qquad
\mathbb{E}[g(Y(0))\mid C] \;=\; \frac{\mathbb{E}\!\left[g(Y)\,\frac{1-Z}{1-e(X)}\big(p_1(X)-D\big)\right]}{p_C}.
$$

Set $g_y(u)=\mathbf{1}\{u\le y\}$ to estimate $F_{1C}$ and $F_{0C}$. Enforce valid CDF shape by applying monotone rearrangement on a fine grid and clipping values to $[0,1]$. Estimate all first‑stage nuisance functions with cross‑fitting. Fit on training folds and evaluate on held‑out folds.

#### Tests

**Uniform CDF test**

Compute

$$
\hat\Delta(y)
=
\hat F_{Y\mid Z=1}(y)
-
\hat F_{Y\mid Z=0}(y)
-
\hat p_C\Big(\hat F_{1C}(y)-\hat F_{0C}(y)\Big).
$$

Use either the Kolmogorov–Smirnov statistic

$$
T_{\infty}=\sup_y \big|\hat\Delta(y)\big|
$$

or the Cramér–von Mises statistic

$$
T_{2}=\int \hat\Delta(y)^2\, d\hat H(y),
$$

where $\hat H$ is a pooled empirical measure on a grid. With covariates, obtain critical values using a multiplier bootstrap that holds fold‑specific nuisance estimates fixed.

**GMM moment test**

For indicator basis $g_j(y)=\mathbf{1}\{y\le t_j\}$, define moments

$$m_j = \Big(\mathbb{E}[g_j(Y)\mid Z{=}1]-\mathbb{E}[g_j(Y)\mid Z{=}0]\Big) - p_C\Big(\mathbb{E}[g_j(Y(1))\mid C]-\mathbb{E}[g_j(Y(0))\mid C]\Big)$$

Stack $m=(m_1,\dots,m_J)$. Use a heteroskedastic‑robust GMM $J$‑test

$$T_J = n\, \hat m^{\top} \hat W^{-1}\hat m$$

with a sandwich covariance $\hat W$. This test targets chosen regions of the distribution through the cutpoints $t_j$.

**Placebo falsification using predicted noncompliance**

Define the compliance propensity $c(X)=p_1(X)-p_0(X)$. In regions where $c(X)$ is near zero, LATE implies a null ITT on $Y$. Estimate $c(X)$ with cross‑fitting. Define smoothed weights

$$
w_{\tau}(X)=K\!\left(\frac{\hat c(X)-\tau}{h}\right)
$$

for a kernel $K$ and bandwidth $h$. Test that the weighted conditional mean difference of $Y$ across $Z$ is zero. Construct a max statistic over a grid of thresholds $\tau$. Use a multiplier bootstrap to obtain critical values.

#### Simulation study

The core design is a randomized trial without covariates. Principal strata are drawn with shares that sum to one. Assignment is Bernoulli with rate one half. Potential treatment equals type rules. The baseline outcome $Y(0)$ is standard normal. Complier treatment effects are normal with mean $0.5$ and standard deviation $0.4$. Observed outcomes satisfy exclusion in the null. Four classes of scenarios are studied.

* **Null** that satisfies the LATE conditions with complier shares $p_C=0.30$ and $p_C=0.10$.
* **Exclusion violations** that add a direct effect $\gamma Z$ to $Y$ for noncompliers and a larger variant with $\gamma=0.5$.
* **Defiers** with shares five percent and ten percent.
* **Sensitivity** to sample size and complier share.

The sample size is two thousand unless stated. For each null, the 95th percentiles of test statistics are estimated by Monte Carlo and used as fixed critical values to evaluate empirical size and power.

**Main scenarios, $n=2000$.** Rejection frequencies are at the nominal five percent under the null and rise with violation severity.

| Scenario               | $p_C$ |  Alt | $T_2$ mean | $T_2$ sd | Reject $T_2$ | Reject $T_J$ |
| ---------------------- | :---: | :--: | :--------: | :------: | :----------: | :----------: |
| Null                   |  0.30 | none |  0.000571  | 0.000543 |     0.056    |     0.056    |
| Null                   |  0.10 | none |  0.000795  | 0.000665 |     0.056    |     0.052    |
| Exclusion $\gamma=0.2$ |  0.30 | excl |  0.000779  | 0.000593 |     0.100    |     0.092    |
| Exclusion $\gamma=0.5$ |  0.30 | excl |  0.007091  | 0.002006 |     1.000    |     1.000    |
| Defiers five percent   |  0.30 |  def |  0.000968  | 0.000690 |     0.140    |     0.132    |
| Defiers ten percent    |  0.30 |  def |  0.001470  | 0.000943 |     0.312    |     0.336    |

**Sensitivity.** Power improves with $n$ and varies with $p_C$ in violation‑specific ways.

| Scenario               |  $n$ | $p_C$ | $T_2$ mean | $T_2$ sd | Reject $T_2$ | Reject $T_J$ |
| ---------------------- | :--: | :---: | :--------: | :------: | :----------: | :----------: |
| Exclusion $\gamma=0.2$ | 1000 |  0.30 |  0.000992  | 0.000913 |     0.060    |     0.064    |
| Exclusion $\gamma=0.2$ | 4000 |  0.30 |  0.000580  | 0.000363 |     0.084    |     0.088    |
| Defiers five percent   | 2000 |  0.10 |  0.001518  | 0.001014 |     0.180    |     0.176    |
| Defiers five percent   | 2000 |  0.50 |  0.000562  | 0.000552 |     0.080    |     0.076    |

The uniform CDF test holds size near nominal under the LATE null. It has strong power against large exclusion effects and rising power against defiers as their share grows. The GMM test behaves similarly and allows focus on tails through the choice of cutpoints.

#### How to run

Open JupyterLab and run the replication cell provided in the repository. 

