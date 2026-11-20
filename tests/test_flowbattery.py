import math
import unittest
from dataclasses import replace

import numpy as np

from electrochem_toolkit import flowbattery as fb


class FlowBatteryTests(unittest.TestCase):
    def setUp(self):
        self.cell = fb.FlowCell()
        self.el = fb.Electrolyte()

    def test_ocv_is_formal_potential_at_half_charge(self):
        n = self.el.volume_m3 * self.el.vanadium_mol_m3 / 2
        v, _ = fb.cell_voltage(self.cell, [n, n, n, n], self.el.volume_m3, 0.0, 5e-3)
        self.assertAlmostEqual(v, self.cell.formal_potential_V, places=12)

    def test_nernst_slope(self):
        # Both half-cells at the same SOC: E = E0' + 2RT/F ln(s/(1-s)).
        n, s = self.el.volume_m3 * self.el.vanadium_mol_m3, 0.8
        state = [s * n, (1 - s) * n, (1 - s) * n, s * n]
        v, _ = fb.cell_voltage(self.cell, state, self.el.volume_m3, 0.0, 5e-3)
        expected = self.cell.formal_potential_V + 2 * fb.R * 298.15 / fb.F * math.log(s / (1 - s))
        self.assertAlmostEqual(v, expected, places=12)

    def test_crossover_conserves_vanadium_and_electrons(self):
        state = np.array([0.03, 0.05, 0.06, 0.02])
        rates = fb.crossover_rates(self.cell, state, self.el.volume_m3)
        self.assertAlmostEqual(rates.sum(), 0.0, places=18)
        self.assertAlmostEqual(np.dot([2, 3, 4, 5], rates), 0.0, places=18)

    def test_charge_counting_without_crossover(self):
        cell = replace(self.cell, diffusivity_m2_s=(0, 0, 0, 0))
        start = self.el.initial_state()
        ch = fb.half_cycle(cell, self.el, start, 10.0, 5e-3, 1.65)
        moved = ch["state"][0] - start[0]
        self.assertAlmostEqual(moved * fb.F / 3600, ch["charge_Ah"], places=6)
        self.assertAlmostEqual(ch["V"][-1], 1.65, places=3)

    def test_no_crossover_gives_unit_coulombic_efficiency(self):
        cell = replace(self.cell, diffusivity_m2_s=(0, 0, 0, 0))
        # The first cycle starts below the discharge cut-off; the second starts and
        # ends on the same cut-off, so every coulomb in must come back out.
        _, dis, _ = fb.cycle(cell, self.el, self.el.initial_state(), 10.0, 5e-3)
        _, _, m = fb.cycle(cell, self.el, dis["state"], 10.0, 5e-3)
        self.assertAlmostEqual(m["coulombic_efficiency"], 1.0, delta=1e-3)
        self.assertGreater(m["voltage_efficiency"], 0.5)
        self.assertLess(m["voltage_efficiency"], 1.0)

    def test_crossover_lowers_coulombic_efficiency(self):
        cell0 = replace(self.cell, diffusivity_m2_s=(0, 0, 0, 0))
        st = self.el.initial_state()
        ce0 = fb.cycle(cell0, self.el, st, 10.0, 5e-3)[2]["coulombic_efficiency"]
        ce1 = fb.cycle(self.cell, self.el, st, 10.0, 5e-3)[2]["coulombic_efficiency"]
        self.assertLess(ce1, ce0 - 0.02)

    def test_single_pass_supply_limit(self):
        n = self.el.volume_m3 * self.el.vanadium_mol_m3
        state = [0.5 * n, 0.5 * n, 0.5 * n, 0.5 * n]
        # Flow so low that one pass cannot carry the reactant the current needs.
        u = 10.0 / (fb.F * 800 * self.cell.flow_section_m2) * 0.5
        v, _ = fb.cell_voltage(self.cell, state, self.el.volume_m3, 10.0, u)
        self.assertTrue(math.isnan(v))

    def test_pump_power_scales_with_velocity_squared(self):
        p1, p2 = self.cell.pump_power_W(1e-3), self.cell.pump_power_W(2e-3)
        self.assertAlmostEqual(p2 / p1, 4.0, places=12)

    def test_remix_recovers_pure_crossover_loss(self):
        state = self.el.initial_state()
        first = None
        for _ in range(30):
            _, dis, m = fb.cycle(self.cell, self.el, state, 10.0, 5e-3)
            state = dis["state"]
            first = first or m["discharge_capacity_Ah"]
        faded = m["discharge_capacity_Ah"]
        recovered = fb.cycle(self.cell, self.el, fb.remix(state), 10.0, 5e-3)[2]["discharge_capacity_Ah"]
        self.assertLess(faded, first * 0.98)
        self.assertAlmostEqual(recovered, first, delta=0.01 * first)

    def test_remix_keeps_oxidation_state_shift(self):
        n = 0.16
        state = np.array([0.0, n / 2, n / 2 - 0.01, 0.01]) + np.array([1e-6, 0, 0, 1e-6])
        mixed = fb.remix(state)
        self.assertAlmostEqual(mixed.sum(), state.sum(), places=12)
        self.assertAlmostEqual(np.dot([2, 3, 4, 5], mixed), np.dot([2, 3, 4, 5], state), places=12)
        self.assertAlmostEqual(mixed[:2].sum(), mixed[2:].sum(), places=12)

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            fb.FlowCell(area_m2=0)
        with self.assertRaises(ValueError):
            fb.Electrolyte(initial_soc=1.0)
        with self.assertRaises(ValueError):
            fb.cycle(self.cell, self.el, self.el.initial_state(), -1.0, 5e-3)


if __name__ == "__main__":
    unittest.main()
