import math
import unittest

import numpy as np

from electrochem_toolkit import corrosion as co

T = co.Tafel
ACID = co.Electrode((T(5e-5, -0.50, 0.060, True), T(5e-5, -0.50, 0.120, False)))
O2 = T(1e-6, -0.30, 0.120, False, i_lim=3e-5)
NEUTRAL = co.Electrode((T(1e-6, -0.75, 0.070, True), O2, T(1e-6, -0.95, 0.120, False)))


class CorrosionTests(unittest.TestCase):
    def test_astm_g102_carbon_steel(self):
        # 1 uA/cm2 on iron (EW 27.92, rho 7.87) is 0.0116 mm/y.
        cr = co.corrosion_rate_mm_per_year(1e-6, **co.CARBON_STEEL)
        self.assertAlmostEqual(cr, 3.27e-3 * 27.92 / 7.87, places=12)
        self.assertAlmostEqual(cr, 0.0116, places=4)

    def test_corrosion_point(self):
        ecorr, icorr = ACID.corrosion_point()
        self.assertAlmostEqual(ecorr, -0.50, places=9)
        self.assertAlmostEqual(icorr, 5e-5, places=12)
        e, _ = NEUTRAL.corrosion_point()
        self.assertAlmostEqual(float(NEUTRAL.faradaic(e)), 0.0, places=14)

    def test_diffusion_plateau(self):
        self.assertAlmostEqual(float(-O2.current(-2.0)), 3e-5, places=12)

    def test_stern_geary_coefficient(self):
        self.assertAlmostEqual(co.stern_geary_B(0.12, 0.12), 0.12 / (2 * math.log(10)), places=15)
        self.assertAlmostEqual(co.stern_geary_B(0.07, math.inf), 0.07 / math.log(10), places=15)

    def test_lpr_on_clean_activation_data(self):
        E = np.arange(-0.52, -0.48 + 1e-9, 0.0005)
        i = co.measured_current(ACID, E)
        rp, _ = co.polarisation_resistance(E, i, half_width=0.005)
        icorr = co.stern_geary_B(0.060, 0.120) / rp
        # Straight-line fit to an exponential over ±5 mV: under 1 % error.
        self.assertAlmostEqual(icorr / 5e-5, 1.0, delta=0.01)

    def test_tafel_extrapolation_on_clean_data(self):
        E = np.arange(-0.85, -0.15, 0.001)
        i = co.measured_current(ACID, E)
        # Far from Ecorr the opposite branch is negligible: exact recovery.
        far = co.tafel_extrapolation(E, i, window=(0.15, 0.30))
        self.assertAlmostEqual(far.icorr / 5e-5, 1.0, delta=0.002)
        self.assertAlmostEqual(far.beta_a, 0.060, delta=0.0005)
        self.assertAlmostEqual(far.beta_c, 0.120, delta=0.0005)
        # The usual 50-150 mV window still carries some of the anodic current on
        # the cathodic branch, which biases icorr low by about 5 % here.
        near = co.tafel_extrapolation(E, i)
        self.assertAlmostEqual(near.icorr / 5e-5, 0.952, delta=0.005)

    def test_ohmic_drop_is_self_consistent(self):
        el = co.Electrode(ACID.reactions, r_u=200.0)
        E_app = np.linspace(-0.7, -0.3, 41)
        i = co.measured_current(el, E_app)
        E_int = E_app - i * el.r_u
        np.testing.assert_allclose(ACID.faradaic(E_int), i, rtol=1e-9, atol=1e-15)

    def test_scan_adds_capacitive_offset(self):
        E = np.array([-0.5])
        el = co.Electrode(ACID.reactions, c_dl=2e-4)
        self.assertAlmostEqual(float(co.measured_current(el, E, 1e-3)[0]), 2e-7, places=15)

    def test_fit_recovers_icorr_and_ru(self):
        el = co.Electrode(NEUTRAL.reactions, r_u=800.0)
        E, i, s = co.potentiodynamic_scan(el, seed=3)
        _, true_icorr = el.corrosion_point()
        fit = co.fit_mixed_potential(E, i, s, "diffusion")
        self.assertAlmostEqual(10 ** fit.parameters["log10_icorr"] / true_icorr, 1.0, delta=0.06)
        self.assertAlmostEqual(fit.parameters["Ru_ohm_cm2"] / 800.0, 1.0, delta=0.05)

    def test_model_selection_prefers_diffusion_limit(self):
        E, i, s = co.potentiodynamic_scan(NEUTRAL, seed=4)
        a = co.fit_mixed_potential(E, i, s, "activation")
        d = co.fit_mixed_potential(E, i, s, "diffusion")
        self.assertLess(d.aicc + 100, a.aicc)

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            T(0, 0, 0.1, True)
        with self.assertRaises(ValueError):
            co.corrosion_rate_mm_per_year(-1, 27.92, 7.87)
        with self.assertRaises(ValueError):
            co.tafel_extrapolation(np.array([0.0, 1.0]), np.array([1.0, 2.0]))


if __name__ == "__main__":
    unittest.main()
