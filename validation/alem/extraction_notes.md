# Alem 2015 experimental extraction

The extracted series come from the lowest tested feed concentration, 0.25 g/L, under constant volumetric flow and constant head. No parameter fitting was performed here.

## Source provenance

Alem, A.; Ahfir, N.-D.; Elkawafi, A.; Wang, H.Q. Hydraulic Operating Conditions and Particle Concentration Effects on Physical Clogging of a Porous Medium. Transport in Porous Media 2015, 106(2), 303–321. DOI: 10.1007/s11242-014-0402-8.

Publisher: https://link.springer.com/article/10.1007/s11242-014-0402-8

Readable primary author manuscript: https://www.academia.edu/52371161/Alem_et_al_2014_Transp_Porous_Med

Native Figure 2: https://media.springernature.com/full/springer-static/image/art%3A10.1007%2Fs11242-014-0402-8/MediaObjects/11242_2014_402_Fig2_HTML.gif

Native Figure 6: https://media.springernature.com/full/springer-static/image/art%3A10.1007%2Fs11242-014-0402-8/MediaObjects/11242_2014_402_Fig6_HTML.gif

The figure endpoints return PNG bytes although their URLs end in GIF. File magic, image dimensions and SHA256 hashes are recorded in metadata.json. Full-width and lw1200 requests returned identical 427-pixel native width; enlarged audit images add no source resolution. No complete paper PDF or raw experimental spreadsheet was obtained.

## What was extracted

- Figure 2a: 26 open-circle experimental concentration markers, constant flow.
- Figure 2b: 29 filled-circle experimental concentration markers, constant head.
- Figure 6a: 7 open-circle hydraulic-ratio locations, constant flow; the final four overlap the 0.50 g/L series and are flagged.
- Figure 6b: 7 filled-circle hydraulic-ratio markers, constant head.
- Table 1: two numerical end-of-run mass balances at 83 injected pore volumes. The constant-flow retained fraction is 0.6795 and should be retained as a held-out target, not used in the first concentration/hydraulic fit.

The figures present experimental symbols with visual connectors. Only identifiable symbol centers were sampled; no points were generated along connecting lines. Congested startup symbols and ambiguous parts of the constant-head concentration series were omitted. The CSV intentionally has gaps rather than artificial observations. In particular, no constant-head concentration points between approximately 6 and 26 pore volumes were invented.

C/C0 is the calibrated-turbidity-based total effluent concentration relative to the feed. K/K0 is the authors' experimental hydraulic ratio calculated from measured pressures and flow with Darcy's law, not a directly measured intrinsic permeability and not a current-model output. The two particle populations are not separately identified by these aggregate curves.

## Coordinate mapping and uncertainty

Pixels use the top-left image origin. Each CSV row retains its native x,y coordinates. Axes are linear, x=0 to100 pore volumes and y=0 to1:

x_data = 100*(pixel_x-x0)/(x100-x0)

y_data = (y0-pixel_y)/(y0-y1)

All four panel calibrations are in metadata.json and digitize_alem.py. Centers were manually selected and visually checked; isolated filled-circle positions were additionally refined by centroids of selected dark-pixel components. The figures and overlaid audit images let a reviewer assess every selection.

Ordinary points have a heuristic resolution bound of ±2 native pixels, including marker-center and axis placement, corresponding to about ±0.55 pore volumes and ±0.0075 in C/C0 or K/K0. The four overlapping flow hydraulic points have ±4-pixel vertical bounds, about ±0.0149. These are digitization bounds, not confidence intervals or experimental standard deviations. The source does not provide individual observation error bars. Small apparent nonmonotonicities are preserved rather than smoothed.

For conservative analysis, exclude the four overlapping hydraulic rows or report sensitivity to including them. Do not assign these points greater precision than allowed by the native raster.

## Parameters and conversions

Reported geometry, operating conditions, porosity, density and initial hydraulic conductivity are in metadata.json. The reported 0.25 g/L divided by particle density2600 g/L gives total injected particle volume fraction9.61538e-5. The reported initial hydraulic conductivity has unitsm/s; it must not be inserted directly as intrinsic permeability in m2. Converting to intrinsic permeability requires fluid viscosity, density and gravity, with any assumed fluid temperature declared.

The two component inlet fractions, fixed deposit porosities and kinetic thresholds are not provided by this experiment. They cannot be claimed measured. Vertical hydraulic-head differences must be treated as hydraulic potentials consistently with the model pressure variable.

## Table 1 and profile limits

The source reports constant-flow retained/injected0.6795, effluent/injected0.3162, and total recovery0.9958 at83PV; corresponding constant-head values are0.7394,0.1964 and0.9358. These are transcribed table values. The0.0001 difference between the rounded constant-flow sum and reported total is preserved. Experimental mass discrepancies are distinct from numerical conservation errors.

A spatial retained-mass profile has now been added for the 0.25 g/L constant-flow run in data_retention_profile_table1_flow.csv. Its nine Table 1 percentages were independently checked against the primary author manuscript: 9.86, 16.94, 11.31, 8.34, 6.08, 5.85, 3.73, 3.69 and 2.15% of injected mass. They sum to the reported retained fraction 0.6795. Section 2 specifies consecutive sample widths of 2 cm, seven intervals of 5 cm, and 3 cm, giving boundaries at 0, 2, 7, 12, 17, 22, 27, 32, 37 and 40 cm from the inlet. These are reported sample dimensions, not lengths inferred from schematic pixels.

The nominal mass profile converts, using the reported 83 injected pore volumes, total feed particle volume fraction 0.25/2600 and the section widths, to a length-weighted mean deposited-solid fraction of 0.0054229326923 per original pore volume. Figure 5 was not digitized to create this profile: the values are numerical table observations, with provenance recorded in taps_and_profile_check.md. Only the constant-flow section values were independently checked and used. The constant-head section row was not used because its author-page OCR values have an unresolved discrepancy with the reported aggregate; its global retained fraction remains separately transcribed from Table 1, with that limitation preserved.

## Scope for validation

This is a broad kaolinite size distribution, not a controlled binary feed. A two-bin interpretation needs explicit assumptions or further size-resolved data. Keep the seven original balances, fixed deposit porosities, existing liquid exchange and permeability relation unchanged. An unsuccessful joint fit is evidence about the unchanged framework and must not be repaired by altering observations or adding prohibited physics.

Scripts and CSVs are extraction/reproducibility materials. Source figure rasters and audit images are internal evidence; their inclusion in a published supplement would require a rights check. Redraw the validation comparison from the extracted observations instead of republishing source figures.

Reported hydraulic baselines require care: Q=145mL/min and column diameter4.5cm imply initial Darcy velocity0.0015195m/s. Using reported K0=0.0007m/s and interpreting105cm as the total-head drop over40cm gives0.0018375m/s. These values do not coincide. Resolve pressure-head and total-head references, or use a clearly declared initial hydraulic calibration; do not claim all values agree exactly.
