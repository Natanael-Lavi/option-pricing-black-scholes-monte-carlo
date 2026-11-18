# European Option Pricing: Black-Scholes vs Monte Carlo

A small Python project that prices European call and put options using two methods:

1. The analytical Black-Scholes formula.
2. Monte Carlo simulation under the risk-neutral measure.

The project compares the numerical estimates with the analytical prices and studies how the Monte Carlo error changes as the number of simulations increases.

## Model and assumptions

The underlying asset follows a geometric Brownian motion under the risk-neutral measure:

\[
 dS_t = rS_t\,dt + \sigma S_t\,dW_t.
\]

At maturity, its exact distribution is:

\[
 S_T = S_0\exp\left((r-\tfrac12\sigma^2)T + \sigma\sqrt{T}Z\right),
 \qquad Z\sim\mathcal N(0,1).
\]

The option value is the discounted expected payoff:

\[
 V_0 = e^{-rT}\,\mathbb E[\text{payoff}(S_T)].
\]

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
| \(S_0\) | 100 | Initial underlying price |
| \(K\) | 105 | Strike price |
| \(T\) | 1 year | Time to maturity |
| \(r\) | 3% | Continuously compounded risk-free rate |
| \(\sigma\) | 20% | Annualised volatility |
| \(N\) | 1,000,000 | Number of Monte Carlo simulations |

The program prints the Black-Scholes price, the Monte Carlo estimate and an approximate 95% confidence interval for both calls and puts. With the fixed random seed, the results are reproducible.

It also checks put-call parity:

\[
 C-P=S_0-Ke^{-rT}.
\]

## Convergence experiment

Running the program generates `convergence.png` with two panels:

- **Left:** one reproducible Monte Carlo estimate for each sample size, its approximate 95% pointwise confidence interval, and the Black-Scholes reference price.
- **Right:** the root mean squared error (RMSE) calculated over 20 independent Monte Carlo runs for each sample size. A reference line with slope \(-1/2\) illustrates the expected \(N^{-1/2}\) error scale.

The RMSE is used instead of the absolute error from a single run because an individual simulation error can move up or down when the number of paths increases. The theoretical convergence rate describes the statistical scale of the error, not a guarantee that every individual estimate improves monotonically.

The confidence intervals use the normal approximation:

\[
 \widehat V_N \pm 1.96\,\widehat{SE}_N.
\]

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

## Acknowledgement

An AI assistant was used as a coding aid during the initial development. The formulas, assumptions, implementation and tests should be understood and checked by the project author.
