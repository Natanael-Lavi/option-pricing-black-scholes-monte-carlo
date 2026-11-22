# European Option Pricing: Black-Scholes vs Monte Carlo

A small Python project that prices European call and put options using two methods:

1. The analytical Black-Scholes formula.
2. Monte Carlo simulation under the risk-neutral measure.

The project compares the numerical estimates with the analytical prices and studies how the Monte Carlo error changes as the number of simulations increases.

## Motivation

I chose this project to compare two approaches to the same problem: an analytical solution through Black-Scholes and a numerical approximation through Monte Carlo simulation. Because the exact price is known, I could go beyond comparing prices and measure how reliable the numerical estimates are and how their accuracy evolves as the number of simulations increases.

## Model and assumptions

The underlying asset follows a geometric Brownian motion under the risk-neutral measure:

$$
dS_t = r \cdot S_t \cdot dt + \sigma \cdot S_t \cdot dW_t.
$$

At maturity, its exact distribution is:

$$
S_T = S_0 \cdot \exp\left(\left(r-\frac{\sigma^2}{2}\right) \cdot T
+\sigma\sqrt{T}\cdot Z\right) 
\qquad Z\sim\mathcal{N}(0,1).
$$

The option value is the discounted expected payoff:

$$
V_0 = e^{-rT} \cdot \mathbb{E}[\mathrm{payoff}(S_T)].
$$

The implementation assumes:

- European exercise (the option can only be exercised at maturity).
- No dividends on the underlying asset.
- Constant risk-free interest rate and volatility.
- A lognormal underlying price under the risk-neutral measure.

The Monte Carlo method simulates terminal prices directly. This is sufficient for standard European calls and puts because their payoffs depend only on the terminal stock price, not on the path followed before maturity.

## Example parameters and results

The default example uses:

| Parameter | Value | Meaning |
|---|---:|---|
| $S_0$ | 100 | Initial underlying price |
| $K$ | 105 | Strike price |
| $T$ | 1 year | Time to maturity |
| $r$ | 3% | Continuously compounded risk-free rate |
| $\sigma$ | 20% | Annualised volatility |
| $N$ | 1 million | Number of Monte Carlo simulations |

Results with one million paths and a fixed seed (the program prints them, and the results are reproducible):

| Option | Black-Scholes | Monte Carlo | 95% CI half-width | Absolute difference |
|---|---:|---:|---:|---:|
| Call | 7.1281 | 7.1320 | ±0.0246 | 0.0039 |
| Put | 9.0248 | 9.0250 | ±0.0217 | 0.0001 |

These differences come from one run with one seed, so they say little about accuracy: the convergence experiment below uses 20 independent runs for each sample size.

The program also checks put-call parity:

$$
C-P=S_0-K \cdot e^{-rT}.
$$

## Convergence experiment

Running the program generates `convergence.png` with two panels:

- **Left:** one reproducible Monte Carlo estimate for each sample size, its approximate 95% pointwise confidence interval, and the Black-Scholes reference price.
- **Right:** the root mean squared error (RMSE) calculated over 20 independent Monte Carlo runs for each sample size. A reference line with slope $-1/2$ illustrates the expected $N^{-1/2}$ error scale.

The RMSE is used instead of the absolute error from a single run because an individual simulation error can move up or down when the number of paths increases. The theoretical convergence rate describes the statistical scale of the error, not a guarantee that every individual estimate improves monotonically.

The confidence intervals use the normal approximation:

$$
\widehat{V}_N \pm 1.96 \cdot \widehat{SE}_N.
$$

![Monte Carlo convergence experiment](convergence.png)

## Key observations

For the chosen parameters and one million simulated paths, the Monte Carlo estimates are close to the analytical Black-Scholes prices, and both differences fall within the corresponding approximate 95% confidence intervals (see the table above).

One point I found particularly interesting is that increasing the number of simulations does not necessarily make the error smaller at every individual step. This is why I used the RMSE over 20 independent runs to study convergence, rather than relying on the absolute error of a single simulation.

The experiment illustrates how Monte Carlo methods can approximate an analytical price, while also showing the role of sampling uncertainty in numerical estimation.

## Why many simulations are needed: most paths pay nothing

The call is out of the money ($K = 105$), so under the risk-neutral measure only $\mathbb{P}(S_T > K) = N(d_2) \approx 42.3\%$ of paths end with a positive payoff. The other 57.7% pay exactly 0. The price of 7.13 is therefore carried by a minority of paths with a long right tail, which is why a precise estimate requires many simulations.

![Payoff distribution](payoff_distribution (1).png)

*Discounted call payoff over one million simulated paths (log scale on the vertical axis).*

## Installation

Python 3.10 or later is recommended.

```bash
python -m venv .venv
source .venv/bin/activate        # macOS / Linux
# .venv\Scripts\activate       # Windows
python -m pip install -r requirements.txt
```

## Run the project

```bash
python pricing.py
```

This prints the call and put prices, checks put-call parity and saves `convergence.png` in the current directory.

## Run the tests

```bash
python -m unittest discover -s tests -v
```

The tests cover reference prices, put-call parity, zero-maturity and zero-volatility cases, invalid inputs, reproducibility and consistency between Monte Carlo estimates and Black-Scholes prices.

## Limitations

This is an educational implementation of the standard Black-Scholes model. It does not include dividends, transaction costs, early exercise, stochastic volatility or jumps. The confidence intervals rely on a normal approximation and are most informative when the number of simulations is sufficiently large.

The numerical results are specific to one set of parameters (an out-of-the-money call with $K = 105$); other strikes or volatilities would change the size of the Monte Carlo error. The simulation uses the plain estimator, without variance reduction.

## Acknowledgements

AI assistants were used during development for code structure, tests and drafting of the documentation. I ran the experiments, checked the numerical results against the closed-form price, and I am responsible for the content of this repository.
