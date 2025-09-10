## LATE IV: Distributional Implications of LATE

Many randomized encouragement designs have imperfect compliance, where only a fraction of people comply with their assignment. Examples include phone-bank get-out-the-vote (GOTV) campaigns and draft lotteries like the Vietnam Draft Lottery. In these settings, it is common to use instrumental variable (IV) regression for analysis.

Instrumental variables identify the Local Average Treatment Effect (LATE) under three standard conditions: random assignment, exclusion, and monotonicity. The exclusion restriction embodies a key intuition: we don't expect GOTV efforts over the phone to affect people we aren't able to reach, nor do we expect people who were merely drawn up in the Vietnam Draft Lottery to have different attitudes towards minorities solely from being drafted—whatever effects we see, we expect them to result from actual service.

This leads to a sharp testable implication. When the instrument Z flips, only compliers can change treatment status. Therefore, only compliers can contribute to the treatment effect, which implies that the distribution of the intention-to-treat (ITT) effect follows a specific pattern—one with a lumpy distribution reflecting the complier-only effects.

We turn this insight into testable implications about distributional equalities. Under LATE assumptions, the difference in outcome CDFs between instrument values must equal the complier share times the difference in complier potential outcome CDFs. Formally, for all y: $F_{Y|Z=1}(y) - F_{Y|Z=0}(y) = p_C[F_{1C}(y) - F_{0C}(y)]$. This extends LATE from a statement about means to a family of restrictions across the entire distribution.

These distributional restrictions provide leverage for detecting violations of the underlying assumptions. Exclusion violations shift the entire distribution, including regions where compliers are absent. Defiers create opposing movements that disturb the expected monotonicity of the CDF difference. The ITT effect distribution should exhibit concentrated movement where compliers lie and zero movement elsewhere—deviations from this pattern signal assumption failures.

We operationalize these insights through two main tests. The first is a uniform test of the distributional equality using either Kolmogorov-Smirnov or Cramér-von Mises statistics. We estimate complier CDFs using Abadie-style weighting, enforce shape restrictions through monotone rearrangement, and handle covariates via cross-fitting. The second is a GMM test that focuses on specific quantiles, which proves useful when violations are expected to concentrate in the tails of the distribution.

A complementary falsification test exploits heterogeneity in compliance propensity. In covariate regions where compliance approaches zero, LATE predicts null ITT effects. We test this implication by examining weighted conditional mean differences across the instrument, with weights concentrated on low-compliance regions. This provides a direct test of exclusion using observable variation.

Simulations validate the approach across realistic scenarios. Under the null with complier shares ranging from 10% to 30%, tests maintain nominal size. Exclusion violations generate detectable distributional distortions, with power approaching one for moderate direct effects (γ = 0.5). When defiers are present, detection power increases with the defier share—five percent defiers yield 14% power, while ten percent defiers yield 31% power at n = 2000.

The method's power depends predictably on sample size and complier share. Doubling the sample size from 2000 to 4000 approximately doubles power for small exclusion violations. The relationship with complier share proves more complex—very low or very high complier shares reduce power against defier alternatives, as the distributional signature becomes harder to distinguish from sampling variation.

These tests complement existing LATE diagnostics. While covariate balance tests check randomization and first-stage F-statistics assess instrument strength, our approach directly examines the exclusion and monotonicity assumptions that are typically untestable. The distributional perspective reveals violations that mean-based tests might miss—for instance, when positive and negative exclusion violations cancel in expectation but distort distributional shape.

In applications where exclusion holds approximately but not exactly, the magnitude of distributional distortions provides a measure of violation severity. This information guides sensitivity analyses and helps researchers assess whether IV estimates are sufficiently reliable for policy conclusions.

Implementation remains straightforward with standard econometric software. The main computational burden comes from bootstrap inference when covariates are present, but this remains manageable even for moderate sample sizes. Cross-fitting prevents overfitting in first-stage estimation while maintaining valid inference.

The broader methodological point concerns the testable implications of identification assumptions. These assumptions often imply restrictions beyond their immediate targets. By developing appropriate tests for these implications, we strengthen our ability to assess when causal identification strategies succeed or fail. The distributional approach presented here represents one avenue for improving the credibility of instrumental variable analyses.

### Simulation study

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
