"""Basic tests for the Black-Scholes and Monte Carlo functions."""

import math
import unittest

from pricing import black_scholes, monte_carlo


class TestBlackScholes(unittest.TestCase):
    def setUp(self):
        self.parameters = (100.0, 105.0, 1.0, 0.03, 0.20)

    def test_reference_call_and_put_prices(self):
        self.assertAlmostEqual(black_scholes(*self.parameters, "call"), 7.128065, places=5)
        self.assertAlmostEqual(black_scholes(*self.parameters, "put"), 9.024846, places=5)

    def test_put_call_parity(self):
        S0, K, T, r, sigma = self.parameters
        call = black_scholes(S0, K, T, r, sigma, "call")
        put = black_scholes(S0, K, T, r, sigma, "put")
        expected_difference = S0 - K * math.exp(-r * T)
        self.assertAlmostEqual(call - put, expected_difference, places=10)

    def test_zero_maturity_returns_intrinsic_value(self):
        self.assertEqual(black_scholes(100, 95, 0, 0.03, 0.2, "call"), 5.0)
        self.assertEqual(black_scholes(100, 105, 0, 0.03, 0.2, "put"), 5.0)

    def test_zero_volatility_returns_deterministic_price(self):
        S0, K, T, r = 100.0, 105.0, 1.0, 0.03
        discounted_strike = K * math.exp(-r * T)
        expected_call = max(S0 - discounted_strike, 0.0)
        expected_put = max(discounted_strike - S0, 0.0)
        self.assertAlmostEqual(black_scholes(S0, K, T, r, 0.0, "call"), expected_call)
        self.assertAlmostEqual(black_scholes(S0, K, T, r, 0.0, "put"), expected_put)

    def test_option_name_is_validated(self):
        with self.assertRaises(ValueError):
            black_scholes(*self.parameters, "straddle")

    def test_invalid_market_parameters_are_rejected(self):
        with self.assertRaises(ValueError):
            black_scholes(-100, 105, 1, 0.03, 0.2, "call")
        with self.assertRaises(ValueError):
            black_scholes(100, 105, -1, 0.03, 0.2, "call")
        with self.assertRaises(ValueError):
            black_scholes(100, 105, 1, 0.03, -0.2, "call")


class TestMonteCarlo(unittest.TestCase):
    def setUp(self):
        self.parameters = (100.0, 105.0, 1.0, 0.03, 0.20)

    def test_estimate_is_consistent_with_black_scholes(self):
        for option in ("call", "put"):
            exact = black_scholes(*self.parameters, option)
            estimate, standard_error = monte_carlo(
                *self.parameters, 200_000, option, seed=123
            )
            self.assertLessEqual(abs(estimate - exact), 4.0 * standard_error)

    def test_seed_makes_results_reproducible(self):
        first = monte_carlo(*self.parameters, 10_000, "call", seed=7)
        second = monte_carlo(*self.parameters, 10_000, "call", seed=7)
        self.assertEqual(first, second)

    def test_at_least_two_paths_are_required(self):
        with self.assertRaises(ValueError):
            monte_carlo(*self.parameters, 1, "call")

    def test_non_integer_path_count_is_rejected(self):
        with self.assertRaises(ValueError):
            monte_carlo(*self.parameters, 10.5, "call")


if __name__ == "__main__":
    unittest.main()
