"""
European option pricing: Black-Scholes closed form vs Monte Carlo simulation.

Model: the stock follows a geometric Brownian motion under the risk-neutral measure
    dS_t = r * S_t dt + sigma * S_t dW_t
whose exact solution at maturity T is
    S_T = S_0 * exp((r - sigma^2 / 2) * T + sigma * sqrt(T) * Z),   Z ~ N(0, 1)

Run:  python pricing.py
"""

import math

import matplotlib.pyplot as plt
import numpy as np


# --------------------------------------------------------------------------
# 1. Black-Scholes closed-form price
# --------------------------------------------------------------------------
def norm_cdf(x: float) -> float:
    """Standard normal cumulative distribution function (via the error function)."""
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def black_scholes(S0, K, T, r, sigma, option="call"):
    """Exact Black-Scholes price of a European call or put."""
    d1 = (math.log(S0 / K) + (r + 0.5 * sigma**2) * T) / (sigma * math.sqrt(T))
    d2 = d1 - sigma * math.sqrt(T)
    if option == "call":
        return S0 * norm_cdf(d1) - K * math.exp(-r * T) * norm_cdf(d2)
    return K * math.exp(-r * T) * norm_cdf(-d2) - S0 * norm_cdf(-d1)


# --------------------------------------------------------------------------
# 2. Monte Carlo price
# --------------------------------------------------------------------------
def monte_carlo(S0, K, T, r, sigma, n_paths, option="call", seed=42):
    """
    Monte Carlo price of a European option.

    Idea: price = exp(-rT) * E[payoff(S_T)], and we estimate the expectation
    by the empirical mean over n_paths simulated terminal prices.

    Returns (price, standard_error).
    """
    rng = np.random.default_rng(seed)
    Z = rng.standard_normal(n_paths)
    S_T = S0 * np.exp((r - 0.5 * sigma**2) * T + sigma * math.sqrt(T) * Z)

    if option == "call":
        payoff = np.maximum(S_T - K, 0.0)
    else:
        payoff = np.maximum(K - S_T, 0.0)

    discounted = math.exp(-r * T) * payoff
    price = discounted.mean()
    std_error = discounted.std(ddof=1) / math.sqrt(n_paths)  # shrinks like 1/sqrt(N)
    return price, std_error


# --------------------------------------------------------------------------
# 3. Convergence study + plots
# --------------------------------------------------------------------------
def convergence_plot(S0, K, T, r, sigma, option="call", filename="convergence.png"):
    exact = black_scholes(S0, K, T, r, sigma, option)

    sizes = np.unique(np.logspace(2, 6, 30).astype(int))  # 100 ... 1,000,000 paths
    prices, errors = [], []
    for n in sizes:
        p, se = monte_carlo(S0, K, T, r, sigma, int(n), option)
        prices.append(p)
        errors.append(se)
    prices, errors = np.array(prices), np.array(errors)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

    # Left: MC estimate with 95% confidence band vs the exact price
    ax1.semilogx(sizes, prices, label="Monte Carlo estimate")
    ax1.fill_between(sizes, prices - 1.96 * errors, prices + 1.96 * errors,
                     alpha=0.25, label="95% confidence interval")
    ax1.axhline(exact, color="red", linestyle="--", label=f"Black-Scholes = {exact:.4f}")
    ax1.set_xlabel("Number of simulated paths")
    ax1.set_ylabel("Option price")
    ax1.set_title(f"Monte Carlo convergence ({option})")
    ax1.legend()
    ax1.grid(alpha=0.3)

    # Right: absolute error vs N, compared with the theoretical 1/sqrt(N) rate
    abs_err = np.abs(prices - exact)
    ax2.loglog(sizes, abs_err, "o-", label="|MC - Black-Scholes|")
    ax2.loglog(sizes, errors, "--", label="Standard error (theory)")
    ax2.set_xlabel("Number of simulated paths")
    ax2.set_ylabel("Error")
    ax2.set_title("Error decreases like 1/sqrt(N)")
    ax2.legend()
    ax2.grid(alpha=0.3, which="both")

    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    print(f"Plot saved to {filename}")


# --------------------------------------------------------------------------
# 4. Main
# --------------------------------------------------------------------------
if __name__ == "__main__":
    # Parameters: change these and re-run to see what happens!
    S0 = 100.0     # current stock price
    K = 105.0      # strike
    T = 1.0        # maturity in years
    r = 0.03       # risk-free rate
    sigma = 0.20   # volatility
    N = 1_000_000  # number of Monte Carlo paths

    for opt in ("call", "put"):
        bs = black_scholes(S0, K, T, r, sigma, opt)
        mc, se = monte_carlo(S0, K, T, r, sigma, N, opt)
        print(f"{opt.upper():4s} | Black-Scholes: {bs:.4f} | "
              f"Monte Carlo: {mc:.4f} (+/- {1.96 * se:.4f}) | "
              f"abs. error: {abs(mc - bs):.4f}")

    # Sanity check: put-call parity  C - P = S0 - K * exp(-rT)
    c = black_scholes(S0, K, T, r, sigma, "call")
    p = black_scholes(S0, K, T, r, sigma, "put")
    print(f"Put-call parity check: C - P = {c - p:.6f} vs "
          f"S0 - K*exp(-rT) = {S0 - K * math.exp(-r * T):.6f}")

    convergence_plot(S0, K, T, r, sigma, "call", filename="convergence.png")
