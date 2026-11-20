# Changelog

## 0.4.0

- Three new studies, each with a model module, run script, figure and tests:
  - `flowbattery`: vanadium flow cell with four-species crossover (Tang, Bao and Skyllas-Kazacos 2011), single-pass mass transfer, Darcy pumping, voltage-event cycling and tank remixing.
  - `corrosion`: mixed-potential polarisation curves with implicit ohmic drop; Tafel extrapolation, LPR with Stern–Geary, a full-curve fit with AICc model selection, and ASTM G102 corrosion rates.
  - `kinetics`: Butler–Volmer sweeps with ohmic drop and charging current, Tafel slope plots, window-error maps and a full-curve fit.
- Moved the toolkit out of my profile repository into its own home, `pooriamn/electrochem-toolkit`. History is preserved.
- Renamed `projects/` to `examples/` and `journal/` to `notes/`; the package README is now the repository README.
- Redrew both study figures with a shared style (`electrochem_toolkit.plotstyle`): panel letters, ±2σ residual band, ΔAICc and α on the EIS figure, thermoneutral shading and energy-closure note on the electrothermal figure. Numerical results are unchanged.
- One-page website with a live impedance explorer, stack calculator and charge-throughput tool, in light and dark themes.
- Fixed the version reported by `__version__` and `CITATION.cff` (both still said 0.2.0).
- CI now tests Python 3.10, 3.12 and 3.13. Release notes are taken from the matching changelog section, and PyPI publishing is ready to switch on.
- Removed the public-activity snapshot workflow and its weekly bot commits.

## 0.3.0

- Added independently integrated charge (Ah) and time-averaged current (A) to experimental CSV reports, while retaining existing report fields and validation rules.
- Applied the same calculation to the installable toolkit and standalone experimental example; added synthetic regression tests for uneven sampling, changing voltage and zero current.
- Clarified endpoint-trapezoid assumptions and the limits of missing-data interpretation. This release does not estimate Faradaic efficiency or gas production.

## 0.2.0

- Advanced EIS inverse modelling with analytic Jacobian, constrained multistart optimisation, RC/RQ comparison, AICc and local identifiability/uncertainty diagnostics.
- Dynamic electrolyser electrothermal model with Nernst pressure correction, apparent Butler–Volmer kinetics, Arrhenius ASR, energy-state integration and temperature-domain events.
- Reproducible synthetic computational studies and numerical verification; optional advanced dependencies.
- Introductory examples retained as quickstarts.

## 0.1.0

First versioned public toolkit.

- Shared installable Python package with impedance, stack, cycling throughput, polarisation and measurement-report modules.
- CLI for stack, impedance, signed-current integration and experimental reports.
- Browser workbench, documented synthetic examples and a public activity summary.
- Numerical regression tests, input validation tests and continuous integration on Python 3.10 and 3.12.

Known limits: idealised DC balances and simple equivalent circuits; no hardware control, general electrochemical simulation or production validation.
