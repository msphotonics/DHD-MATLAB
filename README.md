# DHD-MATLAB

MATLAB code for the illustrative examples and the performance comparison in

> M. F. Spencer, "Digital-Holographic Detection," in *Encyclopedia of Modern
> Optics III*, Elsevier (2027).

The scripts simulate digital-holographic detection with wave optics for four
recording geometries and compare the measured signal-to-noise ratios (SNRs)
and field-estimated Strehl ratios with the closed-form expressions in the
chapter. Appendices C and D of the chapter describe the numerical model and
performance metrics.

## Getting started

1. Download this repository (green **Code** button → **Download ZIP**) or
   clone it with Git.
2. In MATLAB, change to the `scripts` folder.
3. Run one of the scripts, for example

   ```matlab
   cd scripts
   script_IPRG_illustrative
   ```

Each script changes to its own folder, adds `../functions` to the path,
prints its results to the command window, and saves its figures (`.pdf` and
`.png`) and results (`.mat`) to `../outputs/`.

## Scripts

| Script | Chapter | Recording geometry | FPA location | Figures saved |
|---|---|---|---|---|
| `script_IPRG_illustrative.m` | Section 4.6 | Off-axis IPRG | Image plane | `IPRG_*` |
| `script_PPRG_illustrative.m` | Section 5.6 | Off-axis PPRG | Pupil plane | `PPRG_*` |
| `script_PPPSRG_illustrative.m` | Section 6.6 | On-axis PSRG | Pupil plane | `PPPSRG_*` |
| `script_IPPSRG_illustrative.m` | Section 6.6 | On-axis PSRG | Image plane | `IPPSRG_*` |
| `script_comparison.m` | Section 7 | All four | Both | `comparison_2x2` |
| `script_truth_fields.m` | Appendix C | All four | — | `truth_fields` |
| `script_check_epsilon_W.m` | Section 7 | All four | Both | none (prints the encircled-energy check) |

Each illustrative script runs a point-source object and a bulls-eye
(extended) object. `script_comparison.m` sweeps 10 SNR values from 0.1 to 100
with 100 realizations each and takes about 25 minutes on a desktop
computer. Reduce `nReals` near the top of the script for a quick check.

## Functions

The `functions` folder holds the shared helpers:

- **Propagation:** `fresnel_prop_object_to_pupil`, `fresnel_prop_pupil_to_image`,
  `fft2c`, `ifft2c`
- **Objects and aberrations:** `make_isoplanatic_phase`
- **Detection:** `record_holograms`, `detect_dn`, `quantize_dn`, `pixel_bin`,
  `pixel_tf`, `register_pixel_centres`
- **Reconstruction:** `reconstruct_estimate_iprg`, `reconstruct_estimate_pprg`,
  `reconstruct_estimate_psrg`
- **Metrics:** `pixel_eta` (total-system efficiency), `field_strehl`,
  `align_global_phase`
- **Reporting and plotting:** `print_run_header`, `print_object_report`,
  `make_figure`, `show_*`, `add_zoom_inset`, `centered_indices`

## Viewing the PDFs on macOS

MATLAB's vector PDF export draws each image as a tiling pattern, which
macOS Preview and other Apple PDF viewers can render with an offset (the
images and colorbars appear wrapped). Adobe Acrobat and most other viewers
render them correctly. To make a PDF that renders correctly everywhere, run

```
python tools/fix_pdf_patterns.py outputs/IPRG_pt.pdf outputs/IPRG_pt_fixed.pdf
```

The script requires Python and `pypdf` (`pip install pypdf`). It keeps the
same image pixels and only changes how the PDF draws them.

## Requirements

- MATLAB R2020b or later (`tiledlayout`, `exportgraphics`).
- Statistics and Machine Learning Toolbox (optional). With it, the scripts
  draw photon noise from a Poisson distribution with `poissrnd`; without it,
  they detect its absence and use a Gaussian approximation instead.

## Citation

If you use this code, please cite the chapter above.

## License

See [`LICENSE`](LICENSE).
