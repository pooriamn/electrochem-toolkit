# Studies

Each study asks one practical question and answers it on synthetic data with a known answer, so the method can be judged on whether it gets that answer back. Every one has a README with the equations and its limits, a script that regenerates the figure and `summary.json`, and tests that check the model against independent results.

| Study | Question | Methods |
| --- | --- | --- |
| [Impedance inverse modelling](eis-inverse/README.md) | Ideal capacitor or constant-phase element? | Weighted complex least squares, analytic Jacobian, multistart, AICc, local uncertainty |
| [Electrolyser electrothermal model](electrothermal-stack/README.md) | Where does the heat go in a stack under a changing load? | Thermodynamics, apparent Butler–Volmer kinetics, Arrhenius resistance, ODE integration, energy closure |
| [Vanadium flow battery](flow-battery/README.md) | How much electrolyte to pump, and what does remixing recover? | Four-species crossover model, single-pass mass transfer, Darcy pumping, voltage-event cycling |
| [Corrosion rate from a polarisation scan](corrosion/README.md) | When can Tafel extrapolation and LPR be trusted? | Mixed-potential theory, implicit ohmic drop, Stern–Geary, model selection by AICc, ASTM G102 |
| [Tafel analysis](tafel-analysis/README.md) | How much do iR correction and window choice move a Tafel slope? | Butler–Volmer with R_u and charging current, Tafel slope plot, window-error map, full-curve fit |

Run from the repository root with Python 3.10 or newer:

```sh
python -m pip install ".[advanced]"
python studies/eis-inverse/run.py
python studies/electrothermal-stack/run.py
python studies/flow-battery/run.py
python studies/corrosion/run.py
python studies/tafel-analysis/run.py
python -m unittest discover -s tests -v
```

Outputs go to `study-output/`, which git ignores. The copies of `results.svg` and `summary.json` in each folder are what the READMEs show.
