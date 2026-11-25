# Gitis (2010): audited experimental-symbol digitization

The CSVs contain **27 Figure 6 symbols and 60 Figure 7 symbols** at 0.82, 1.02 and 1.62 m. They are digitized experimental data, not original/raw laboratory files. No model curve was sampled, and no kinetic model was fitted here.

The complete article text was read in the cached publisher PDF and coauthor-posted copy. The mirror PDF could not be downloaded locally (HTTP 503). Actual figure images were recovered from coauthor Isaak Rubinstein's Academia copy. Original Figure 6 is `figure_006.jpg`; original Figure 7 is `figure_008.jpg` because the site's automated figure numbering differs from the printed captions. The original-image hashes and sizes are recorded in `axis_calibration.json`.

## Method

1. Inspect native JPGs and identify filled circles (0.82 m), open circles (1.02 m), and downward triangles (1.62 m). The printed caption explicitly distinguishes experimental symbols from calculated curves.
2. Manually select visible symbol centers. Filled triangles use the approximate center of their bounding box. Each selected point retains its native pixel coordinates, marker type and series identifier in the CSV.
3. Calibrate linear axes from every labeled major tick. The horizontal graph axis extends left of the time-zero tick: **its left-hand vertical axis is not time zero**. Figure 6's zero tick is approximately x=303.5 px, and Figure 7's is x=264 px.
4. Convert pixels to hours and C/C0 using the affine transformations recorded in `axis_calibration.json`. This is axis calibration only, not physical-model fitting. Figure 6's first symbols align with the displayed zero tick; a calibration residual below 0.01 h is recorded as zero for those three symbols. Other times remain digitized, rather than snapped to an assumed sampling schedule.
5. Visually inspect all selections against the source images using the QA overlays. Both overlays were inspected and every selected point corresponds to an experimental marker. No additional unmarked points were inserted.

`digitize_gitis.py` reproduces the CSVs, calibration JSON, pixel-center JSON and QA images from the two native images. `figure6_experimental_points.csv` is the initial calibration candidate. `figure7_experimental_points.csv` is reserved as a separate-run comparison candidate; its incomplete run metadata limits interpretation as strict validation.

## Digitization uncertainty

A conservative **±3 native pixels** envelope combines manual center selection, tick location, JPEG blur and nearby curve interference. It is an extraction envelope, **not** an experimental standard deviation or confidence interval. Approximate resulting envelopes are:

| Image | Time | C/C0 |
|---|---:|---:|
| Figure 6 | ±0.03576 h | ±0.001538 |
| Figure 7 | ±0.04276 h | ±0.003110 |

Figure 7's dense clusters and some late Figure 6 filled/open circles touch or nearly touch. Pixel centers remain manually selected; the envelope intentionally covers small shifts within those visible marker footprints. Figure 7's central points carry a quality flag. Laboratory measurement scatter exceeds extraction uncertainty in some regions; it must not be replaced by the pixel envelope when assessing statistical fit.

## Conditions and limits

`metadata.json` separates reported experimental values, source model parameters, calculated geometry values and unreported quantities. The measured porosity is **0.44**. The source's generic illustrative numerical porosity 0.4, feed 10 mg/L and particle density 1500 kg/m³ are **not** reported run-specific measurements for Figures 6–7.

Both captions specify 10 m/h superficial velocity and L=1.62 m; the experimental description rounds packed-bed depth to 1.6 m. The experiments used kaolin in Cincinnati tap water and downward contact filtration. Feed mass concentration, particle-size distribution, particle/deposit density and run-specific coagulant conditions are not given in the relevant experiment section. The different fitted coefficients in the two printed captions do not document a controlled change in any experimental input.

The paper explicitly states that these two runs ended before the breakthrough stage. Consequently, these points can assess ripening/operable behavior, but cannot identify terminal retention capacity or a separate aging regime. They provide one particulate-suspension benchmark, not a validation of two species or all seven saturation states.

Figure 6's first symbols occur at displayed time zero with nonzero concentration at all depths. The source does not explain whether the clock starts after feed arrival. The initial advective transit times calculated from porosity and velocity are 0.036, 0.045 and 0.071 h. Any clock-offset or initial-profile treatment in calibration must be declared rather than silently shifting the data.

The source paper's model includes dispersion. A calibration of the authorized dispersion-free model must report any resulting mismatch and its restricted scope. Figure 7 may test transfer to a separate run, but unknown feed/chemistry differences prevent claiming controlled same-condition replication.

## Provenance

- Gitis, V.; Rubinstein, I.; Livshits, M.; Ziskind, G. Chemical Engineering Journal 2010, 163, 78–85. DOI: https://doi.org/10.1016/j.cej.2010.07.044 .
- Coauthor copy: https://www.academia.edu/108262058/Deep_bed_filtration_model_with_multistage_deposition_kinetics .
- Figure 6: https://www.academia.edu/figures/9480255/figure-6-experimental-data-dots-and-model-prediction-curves .
- Figure 7: https://www.academia.edu/figures/9480280/figure-7-experimental-data-dots-and-model-prediction-curves .
- Experimental description: Section 5, pp. 83–84; run coefficients: printed Figure 6/7 captions.

Source images are retained for extraction audit. The full copyrighted paper has not been placed in the deliverable directory.
