# Quantum sensing with interference: data and plotting code

Supporting numerical material for Han Li (`lihan1105`), Group 7, experiment performed on 5 November 2026.

This experiment is kept in `05.11 quantum sensing/` within the shared [Group 7 laboratory repository](https://github.com/lihan1105/Group7_HanLi_Jena_quantum_lab). Later experiments can have separate folders alongside it.

This repository contains the **11 original text files**, the Python code used to generate **all 13 analysis plots** in the laboratory report, and the fitted results. Raw exports are preserved byte-for-byte, including malformed trailing rows. The laboratory photographs and course handout are supplied with the report rather than this numerical repository.

## Reproduce the plots

Use Python 3.12 or newer and a virtual environment. From the repository root, enter this experiment folder first:

```bash
cd "05.11 quantum sensing"
python -m venv .venv
# Windows PowerShell:
.venv\Scripts\Activate.ps1
# macOS/Linux: source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/analyze.py
python scripts/verify.py
```

The dependency versions in `requirements.txt` record the environment used to reproduce the report. No local Windows fonts, TeX installation, or absolute file paths are required for these plots. The script can also run from another working directory. Optional arguments are `--data-dir PATH` and `--output-dir PATH`.

### Files and report correspondence

| Location | Contents |
| --- | --- |
| `raw/data.txt` | Original scan-to-sample log |
| `raw/linear-scan-data_01.txt` to `_10.txt` | Original instrument exports, all 16 columns |
| `raw/SHA256SUMS.txt` | SHA-256 checksums of the 11 original files |
| `scripts/analyze.py` | Data selection, weighted nonlinear fits, uncertainty propagation, plots |
| `scripts/verify.py` | Raw-data integrity and numerical reproduction checks |
| `outputs/glass.png` | Report Figure 1: normalized HOM glass scans |
| `outputs/waveplate.png` | Report Figure 2: normalized HWP scans |
| `outputs/visibility.png` | Report Figure 3: HWP visibility comparison |
| `outputs/fit_01.png` to `fit_10.png` | Report Figures A01-A10: each scan, fitted curve and residuals |
| `outputs/clean_*.csv` | Retained rows, with all original channels |
| `outputs/excluded_rows.csv` | Every excluded row and its original text line number |
| `outputs/results.json` | Fit parameters, covariances, errors and fit statistics |
| `reference/results.json` | Results used in the report, for comparison after rerunning |

## Scan identities

Identities follow `raw/data.txt`, not guesses based on photograph filenames.

| Scan | Measurement |
| --- | --- |
| 01 | Michelson reference, no sample |
| 02 | Michelson with sample |
| 03 | HOM reference, no sample |
| 04 | HOM glass sample 2 |
| 05 | HOM glass sample 1, thinnest |
| 06 | HOM glass sample 3 |
| 07 | HWP 0 degrees, broad scan |
| 08 | HWP 10 degrees |
| 09 | HWP 0 degrees, adjusted fine scan |
| 10 | HWP 23 degrees |

## Data selection and statistical conventions

- The first five lines are instrument headers. Column 1 is stage position in mm; columns 2 and 3 are singles rates; column 6 is coincidence channel 01. Integration time is 1 s. Exported rates are used without the display's `x10` multiplier.
- A retained row has 16 finite entries, absolute stage position below 20 mm, and positive values in both active singles channels and coincidence 01. These checks identify corrupted trailing blocks in these specific exports; they are not a general rule for discarding legitimate zero-count observations.
- Exclusions: scan 03 has 39 rows, scan 04 has 2, and scan 07 has 122. There are 1,081 retained rows and 163 excluded rows. No point is removed because of its fit residual.
- Michelson scans use a Gaussian-envelope cosine. HOM scans use a Gaussian dip. Fits use Poisson uncertainties `sqrt(N)` and `scipy.optimize.curve_fit` with absolute uncertainties.
- Michelson visibility is `A/B`. HOM visibility uses the course convention `A/(2B-A)`, which differs from the fractional dip depth `A/B`.
- Reported covariance is the Poisson fit covariance multiplied by `max(1, reduced chi-square)`. Both enlarged and Poisson-only visibility uncertainties are retained. This is not a complete systematic uncertainty budget.
- The HWP comparison uses scan 09 for the zero-degree fine-scan point. Scan 07 remains available separately. The dashed curve is an ideal comparison anchored at zero degrees, not a fitted curve.

The report discusses the limits of the attenuation model, unresolved calibration errors and nonmonotonic HWP data. Do not treat statistical error bars as absolute measurement accuracy.

## Provenance and AI assistance

The measurements originate from the supplied laboratory exports. AI assisted with the analysis code, numerical fitting, plot generation, documentation and report preparation. No experimental measurement was generated by AI. The student author is responsible for reviewing the analysis and its interpretation.

This repository is shared for inspection and reproduction of the report. No additional reuse license is asserted for the supplied experimental materials.
