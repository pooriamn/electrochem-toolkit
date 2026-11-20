# Website

Plain HTML, CSS and JavaScript with no build step, served by GitHub Pages from this folder.

Preview locally from the repository root:

```sh
python -m http.server 8000 --directory docs
```

The study figures in `assets/` are copies of `studies/*/results.svg`; CI checks they stay identical. `python scripts/brand.py` regenerates the hero, banner and social preview images.
