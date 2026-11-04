# Supernova cosmology: summer-school project

**Ansh Shinde** | Python, observational cosmology, numerical integration and model fitting

Reconstructs a summer-school analysis of Pantheon+SH0ES Type Ia supernovae. The original Jupyter file was unavailable; this version was rebuilt from its PDF export with ChatGPT assistance and rerun on 4 November 2026. The original exercise included instructor-provided scaffolding. This package is a corrected reconstruction, not the original submission.

## Scientific scope

This is an **educational diagonal-error approximation**, conditional on the existing SH0ES calibration. The data authors require full covariance for cosmological inference. The errors below are conditional fit errors and are **not full cosmological uncertainties**. This project does not independently measure the Hubble tension or establish evolving dark energy.

| Fit | H0 (km/s/Mpc) | Omega_m |
|---|---:|---:|
| Both parameters free | 72.974 +/- 0.262 | 0.3508 +/- 0.0187 |
| Fixed matter density, all entries | 73.530 +/- 0.172 | 0.3 fixed |
| Fixed matter density, z < 0.1 | 73.007 +/- 0.277 | 0.3 fixed |
| Fixed matter density, z >= 0.1 | 73.848 +/- 0.218 | 0.3 fixed |

The derived age at the fitted parameters is 12.358 Gyr. At H0=72.97 and Omega_m=0.3 it is 12.919 Gyr. See the notebook for assumptions and diagnostics.

![Hubble diagram](figures/hubble_diagram.png)

## Run

Use Python 3.12 or newer with the tested package versions in `requirements.txt`.

```bash
python -m pip install -r requirements.txt
python supernova_cosmology.py
```

Alternatively open `supernova_cosmology.ipynb` in Jupyter or Colab and run all cells. JupyterLab or the Notebook interface can be installed separately; ipykernel is included here. On a headless machine use `MPLBACKEND=Agg python supernova_cosmology.py`.

Internet is needed for the initial download of the official data. The checksum is verified. The raw dataset is not bundled; subsequent runs use `data/Pantheon+SH0ES.dat`. Execute from this project's directory. Figures and `results.json` are regenerated.

## Corrections and checks

- Removed the unused fitted parameter that caused infinite covariance in the fixed-density fit.
- Corrected the age sensitivity calculations and age-integral explanation.
- Added the model overlay, error bars, and residual bin diagnostics.
- Replaced a machine-specific Windows path and deprecated table-reader argument.
- Clarified calibration dependence, covariance limitations, redshift approximation, and interpretation.
- Checked distances against Astropy, ages against an analytic solution, and fixed-density fits against an exact weighted-intercept solution.

All seven code cells were executed sequentially in a shared Python namespace, with outputs captured in the notebook. A native Jupyter kernel could not be launched because the verification host blocks socket binding. The notebook schema was validated. These checks validate the chosen approximation, not a precision cosmological likelihood. The original PDFs are not included because they contain superseded results and third-party teaching material. The original instructor/summer-school attribution should be completed when available.

## Sources

- [Official data and covariance documentation](https://github.com/PantheonPlusSH0ES/DataRelease/tree/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR)
- [Brout et al. 2022: Cosmological Constraints](https://arxiv.org/abs/2202.04077)
- [Riess et al. 2022: SH0ES](https://arxiv.org/abs/2112.04510)
- [Planck 2018 cosmological parameters](https://arxiv.org/abs/1807.06209)

Data and original teaching material remain subject to their respective authors' terms. No blanket license is asserted over them.
