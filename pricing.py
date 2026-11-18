"""European option pricing with Black-Scholes and Monte Carlo simulation.

The stock follows a geometric Brownian motion under the risk-neutral measure:

    dS_t = r S_t dt + sigma S_t dW_t

Its terminal value has the exact distribution

    S_T = S_0 * exp((r - sigma**2 / 2) * T + sigma * sqrt(T) * Z),
    Z ~ N(0, 1).

Run ``python pricing.py`` to print prices and save ``convergence.png``.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def _validate_option(option: str) -> str:
    """Return a normalized option type, or raise a clear error."""
    if not isinstance(option, str):
        raise ValueError("option must be either 'call' or 'put'.")

    option = option.strip().lower()
    if option not in {"call", "put"}:
        raise ValueError("option must be either 'call' or 'put'.")
    return option


def _validate_parameters(S0: float, K: float, T: float, r: float, sigma: float) -> None:
    """Check that model inputs are finite and economically admissible."""
    values = {"S0": S0, "K": K, "T": T, "r": r, "sigma": sigma}
    for name, value in values.items():
        try:
            is_finite = math.isfinite(value)
        except TypeError as exc:
            raise ValueError(f"{name} must be a finite number.") from exc
        if not is_finite:
            raise ValueError(f"{name} must be a finite number.")

    if S0 <= 0:
        raise ValueError("S0 must be strictly positive.")
    if K <= 0:
        raise ValueError("K must be strictly positive.")
    if T < 0:
        raise ValueError("T cannot be negative.")
    if sigma < 0:
        raise ValueError("sigma cannot be negative.")


def norm_cdf(x: float) -> float:
    """Cumulative distribution function of a standard normal variable.

    erfc is used instead of 1 + erf(x) to reduce cancellation in the
    negative tail of the distribution.
    """
    return 0.5 * math.erfc(-x / math.sqrt(2.0))


def black_scholes(
    S0: float,
    K: float,
    T: float,
    r: float,
    sigma: float,
    option: str = "call",
) -> float:
    """Return the Black-Scholes price of a European call or put.

    Assumptions: no dividends, constant risk-free rate and volatility,
    and a lognormal underlying price under the risk-neutral measure.
    """
    _validate_parameters(S0, K, T, r, sigma)
    option = _validate_option(option)

    # At maturity, the option is worth its intrinsic value.
    if T == 0:
        if option == "call":
            return max(S0 - K, 0.0)
        return max(K - S0, 0.0)

    # With zero volatility, the terminal stock price is deterministic.
    if sigma == 0:
        discounted_strike = K * math.exp(-r * T)
        if option == "call":
            return max(S0 - discounted_strike, 0.0)
        return max(discounted_strike - S0, 0.0)

    volatility_time = sigma * math.sqrt(T)
    d1 = (math.log(S0 / K) + (r + 0.5 * sigma**2) * T) / volatility_time
    d2 = d1 - volatility_time

    if option == "call":
        return S0 * norm_cdf(d1) - K * math.exp(-r * T) * norm_cdf(d2)

    return K * math.exp(-r * T) * norm_cdf(-d2) - S0 * norm_cdf(-d1)


def monte_carlo(
    S0: float,
    K: float,
    T: float,
    r: float,
    sigma: float,
    n_paths: int,
    option: str = "call",
    seed: int | None = 42,
) -> tuple[float, float]:
    """Estimate a European option price and its Monte Carlo standard error.

    Only the terminal stock price is simulated. This is sufficient for a
    standard European option because its payoff depends on S_T only.

    Returns:
        (estimated_price, standard_error)
    """
    _validate_parameters(S0, K, T, r, sigma)
    option = _validate_option(option)

    if isinstance(n_paths, bool) or not isinstance(n_paths, (int, np.integer)):
        raise ValueError("n_paths must be an integer greater than or equal to 2.")
    if n_paths < 2:
        raise ValueError("n_paths must be greater than or equal to 2.")

    rng = np.random.default_rng(seed)
    z = rng.standard_normal(n_paths)
    terminal_stock = S0 * np.exp(
        (r - 0.5 * sigma**2) * T + sigma * math.sqrt(T) * z
    )

    if option == "call":
        payoff = np.maximum(terminal_stock - K, 0.0)
    else:
        payoff = np.maximum(K - terminal_stock, 0.0)

    discounted_payoff = math.exp(-r * T) * payoff
    price = float(discounted_payoff.mean())
    standard_error = float(discounted_payoff.std(ddof=1) / math.sqrt(n_paths))
    return price, standard_error


def convergence_plot(
    S0: float,
    K: float,
    T: float,
    r: float,
    sigma: float,
    option: str = "call",
    filename: str | Path = "convergence.png",
    repetitions: int = 20,
    seed: int = 42,
) -> None:
    """Create a convergence plot using independent replications for RMSE.

    The left panel shows one reproducible sequence of estimates and
    approximate pointwise 95% confidence intervals. The right panel shows
    the RMSE over independent runs for each sample size.
    """
    _validate_parameters(S0, K, T, r, sigma)
    option = _validate_option(option)
    if isinstance(repetitions, bool) or not isinstance(repetitions, int) or repetitions < 2:
        raise ValueError("repetitions must be an integer greater than or equal to 2.")

    exact_price = black_scholes(S0, K, T, r, sigma, option)
    sample_sizes = np.unique(np.logspace(2, 6, 20).astype(int))

    # One reproducible run for the price estimate and confidence intervals.
    estimates = []
    standard_errors = []
    for n_paths in sample_sizes:
        estimate, standard_error = monte_carlo(
            S0, K, T, r, sigma, int(n_paths), option, seed=seed
        )
        estimates.append(estimate)
        standard_errors.append(standard_error)

    estimates = np.asarray(estimates)
    standard_errors = np.asarray(standard_errors)

    # Independent replications provide a more meaningful convergence
    # diagnostic than the absolute error from a single random run.
    rmse_values = []
    for size_index, n_paths in enumerate(sample_sizes):
        replication_errors = []
        for replication in range(repetitions):
            replication_seed = seed + 1_000 + size_index * repetitions + replication
            estimate, _ = monte_carlo(
                S0, K, T, r, sigma, int(n_paths), option, seed=replication_seed
            )
            replication_errors.append(estimate - exact_price)
        rmse_values.append(float(np.sqrt(np.mean(np.square(replication_errors)))))

    rmse_values = np.asarray(rmse_values)

    # A reference curve with slope -1/2, anchored near the middle of the RMSE
    # observations. It illustrates the expected rate, not an exact prediction.
    anchor = len(sample_sizes) // 2
    reference_rate = rmse_values[anchor] * (sample_sizes / sample_sizes[anchor]) ** -0.5

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.8))

    ax1.semilogx(sample_sizes, estimates, marker="o", label="Monte Carlo estimate")
    ax1.fill_between(
        sample_sizes,
        estimates - 1.96 * standard_errors,
        estimates + 1.96 * standard_errors,
        alpha=0.22,
        label="Approx. 95% pointwise CI",
    )
    ax1.axhline(
        exact_price,
        color="red",
        linestyle="--",
        label=f"Black-Scholes = {exact_price:.4f}",
    )
    ax1.set_xlabel("Number of simulated paths")
    ax1.set_ylabel("Option price")
    ax1.set_title(f"Monte Carlo price convergence ({option})")
    ax1.legend()
    ax1.grid(alpha=0.3)

    ax2.loglog(
        sample_sizes,
        rmse_values,
        marker="o",
        label=f"RMSE ({repetitions} independent runs)",
    )
    ax2.loglog(
        sample_sizes,
        reference_rate,
        linestyle="--",
        label=r"Reference slope $N^{-1/2}$",
    )
    ax2.set_xlabel("Number of simulated paths per run")
    ax2.set_ylabel("Root mean squared error")
    ax2.set_title("Convergence across independent runs")
    ax2.legend()
    ax2.grid(alpha=0.3, which="both")

    fig.tight_layout()
    output_path = Path(filename)
    fig.savefig(output_path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"Convergence plot saved to: {output_path}")


def main() -> None:
    """Run the example pricing exercise."""
    S0 = 100.0
    K = 105.0
    T = 1.0
    r = 0.03
    sigma = 0.20
    n_paths = 1_000_000

    print("European option pricing: Black-Scholes vs Monte Carlo")
    print(f"S0={S0}, K={K}, T={T}, r={r:.1%}, sigma={sigma:.1%}, N={n_paths:,}\n")

    for option in ("call", "put"):
        exact = black_scholes(S0, K, T, r, sigma, option)
        estimate, standard_error = monte_carlo(
            S0, K, T, r, sigma, n_paths, option, seed=42
        )
        half_width = 1.96 * standard_error
        print(
            f"{option.upper():4s} | Black-Scholes: {exact:.4f} | "
            f"Monte Carlo: {estimate:.4f} +/- {half_width:.4f} "
            f"(approx. 95% CI) | abs. error: {abs(estimate - exact):.4f}"
        )

    call_price = black_scholes(S0, K, T, r, sigma, "call")
    put_price = black_scholes(S0, K, T, r, sigma, "put")
    parity_left = call_price - put_price
    parity_right = S0 - K * math.exp(-r * T)
    print(
        "\nPut-call parity: "
        f"C - P = {parity_left:.10f}; "
        f"S0 - K*exp(-rT) = {parity_right:.10f}; "
        f"difference = {abs(parity_left - parity_right):.2e}"
    )

    convergence_plot(S0, K, T, r, sigma, option="call")


if __name__ == "__main__":
    main()
