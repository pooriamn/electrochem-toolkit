import math
import unittest

import numpy as np

from electrochem_toolkit import kinetics as k


class KineticsTests(unittest.TestCase):
    def test_tafel_slope_value(self):
        # alpha = 0.5 at 25 °C: 2.303 RT/(0.5 F) = 118.3 mV/dec.
        self.assertAlmostEqual(k.tafel_slope_mV(0.5), 118.32, places=2)

    def test_butler_volmer_limits(self):
        self.assertEqual(float(k.butler_volmer(0.0, 1e-5, 0.5)), 0.0)
        # Far cathodic, the back reaction vanishes and the local slope is the Tafel slope.
        eta = np.array([-0.40, -0.50])
        j = -k.butler_volmer(eta, 1e-5, 0.5)
        slope = 1e3 * (0.10) / (math.log10(j[1]) - math.log10(j[0]))
        self.assertAlmostEqual(slope, k.tafel_slope_mV(0.5), places=3)
        # Near equilibrium it is linear: j = j0 f eta (cathodic, so negative).
        f = k.F / (k.R * 298.15)
        self.assertAlmostEqual(float(k.butler_volmer(-1e-5, 1e-5, 0.5)), -1e-5 * f * 1e-5, delta=1e-14)

    def test_full_ir_correction_recovers_kinetics(self):
        el = k.HerElectrode(c_dl=0.0)
        E, j, _ = k.sweep(el, noise=0.0, floor=0.0)
        eta = k.ir_correct(E, j, el.r_u)
        np.testing.assert_allclose(el.kinetic(eta), j, rtol=1e-9, atol=1e-15)

    def test_clean_tafel_line_after_correction(self):
        el = k.HerElectrode(c_dl=0.0)
        E, j, _ = k.sweep(el, noise=0.0, floor=0.0)
        fit = k.fit_tafel_line(k.ir_correct(E, j, el.r_u), j, j_range=(1e-3, 1e-2))
        self.assertAlmostEqual(fit.slope_mV, k.tafel_slope_mV(0.5), delta=0.05)
        self.assertAlmostEqual(fit.j0 / el.j0, 1.0, delta=0.002)

    def test_uncorrected_resistance_inflates_slope(self):
        el = k.HerElectrode(c_dl=0.0)
        E, j, _ = k.sweep(el, noise=0.0, floor=0.0)
        raw = k.fit_tafel_line(E, j, j_range=(1e-2, 5e-2))
        self.assertGreater(raw.slope_mV, 1.5 * k.tafel_slope_mV(0.5))

    def test_partial_compensation_sits_between(self):
        el = k.HerElectrode(c_dl=0.0)
        E0, j0, _ = k.sweep(el, noise=0.0, floor=0.0)
        E85, j85, _ = k.sweep(el, compensated=0.85, noise=0.0, floor=0.0)
        rng = (1e-3, 1e-2)
        raw = k.fit_tafel_line(E0, j0, j_range=rng).slope_mV
        comp = k.fit_tafel_line(E85, j85, j_range=rng).slope_mV
        self.assertLess(k.tafel_slope_mV(0.5), comp)
        self.assertLess(comp, raw)

    def test_tafel_slope_plot_plateau(self):
        el = k.HerElectrode(c_dl=0.0)
        E, j, _ = k.sweep(el, noise=0.0, floor=0.0)
        x, s = k.local_tafel_slope(k.ir_correct(E, j, el.r_u), j)
        slope, lo, hi = k.plateau_slope(x, s, tolerance=0.01)
        self.assertAlmostEqual(slope, k.tafel_slope_mV(0.5), delta=0.5)
        self.assertGreater(hi - lo, 2.0)

    def test_butler_volmer_fit_recovers_parameters(self):
        el = k.HerElectrode()
        E, j, s = k.sweep(el, seed=5)
        fit = k.fit_butler_volmer(E, j, s)
        p = fit.parameters
        self.assertAlmostEqual(10 ** p["log10_j0"] / el.j0, 1.0, delta=0.02)
        self.assertAlmostEqual(p["alpha_c"], 0.5, delta=0.002)
        self.assertAlmostEqual(p["Ru_ohm_cm2"], el.r_u, delta=0.02)
        self.assertAlmostEqual(p["j_capacitive_A_cm2"], -el.c_dl * 0.005, delta=2e-7)

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            k.butler_volmer(0.1, -1, 0.5)
        with self.assertRaises(ValueError):
            k.sweep(k.HerElectrode(), compensated=1.0)
        with self.assertRaises(ValueError):
            k.fit_tafel_line([0, -0.1], [-1e-3, -2e-3])


if __name__ == "__main__":
    unittest.main()
