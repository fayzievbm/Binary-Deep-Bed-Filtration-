from pathlib import Path
from PIL import Image, ImageDraw
import csv, json, hashlib, math
P=Path(__file__).parent
# Marker centers read manually in native raster pixel coordinates. Origin is top-left.
# Only circle series (0.25 g/L) are selected; congested start-up markers are excluded.
sets={
 'fig2_flow': {'file':'source_fig2.png','figure':'2a','mode':'constant_flow','quantity':'C_over_C0','axis':{'x0':48.5,'x100':417.,'y0':273.,'y1':3.5},'points':[(65.5,229.),(78.,229.),(90.,225.5),(102.,226.),(113.5,223.),(125.,215.5),(136.5,211.5),(148.5,212.),(160.,208.5),(171.5,198.),(183.,191.),(195.,187.5),(207.,183.5),(219.,173.),(231.,173.),(243.,173.),(255.,169.5),(267.,162.),(279.,158.5),(291.,159.),(303.,151.5),(315.,148.5),(327.,148.5),(339.,148.5),(348.,145.5),(354.,145.5)]},
 'fig2_head': {'file':'source_fig2.png','figure':'2b','mode':'constant_head','quantity':'C_over_C0','axis':{'x0':50.5,'x100':416.5,'y0':602.5,'y1':333.},'points':[(65.21,557.5),(73.09,557.22),(148.41,556.74),(155.5,554.5),(163.,555.17),(192.14,551.),(199.22,548.91),(206.76,548.8),(213.61,547.48),(220.91,548.91),(228.12,548.83),(235.5,544.21),(242.52,544.39),(249.52,544.39),(256.5,542.),(263.5,542.),(270.56,541.6),(277.71,542.),(284.5,542.),(291.92,542.08),(298.61,542.04),(305.5,542.),(312.39,542.04),(319.5,542.),(326.,540.),(333.,540.),(340.,540.),(347.,540.),(354.,540.)]},
 'fig6_flow': {'file':'source_fig6.png','figure':'6a','mode':'constant_flow','quantity':'K_over_K0','axis':{'x0':48.5,'x100':414.,'y0':271.5,'y1':3.},'points':[(94.,14.),(150.,17.5),(183.5,18.5),(216.5,19.),(250.,19.5),(284.,19.5),(318.,20.5)]},
 'fig6_head': {'file':'source_fig6.png','figure':'6b','mode':'constant_head','quantity':'K_over_K0','axis':{'x0':50.5,'x100':416.5,'y0':602.5,'y1':334.},'points':[(117.95,339.29),(152.19,337.95),(218.,350.29),(250.09,358.09),(281.6,359.04),(313.,362.86),(353.88,367.17)]}
}
for name,s in sets.items():
 im=Image.open(P/s['file']).convert('RGB').resize((1281,Image.open(P/s['file']).height*3))
 d=ImageDraw.Draw(im)
 for j,(x,y) in enumerate(s['points'],1):
  x*=3;y*=3;d.ellipse((x-9,y-9,x+9,y+9),outline=(255,0,0),width=1);d.line((x-3,y,x+3,y),fill=(255,0,0));d.line((x,y-3,x,y+3),fill=(255,0,0));d.text((x+5,y+4),str(j),fill=(220,0,0))
 im.save(P/(name+'_audit.png'))

# Linear pixel-to-axis mapping. Uncertainty bounds are digitization resolution,
# not reported experimental standard deviations.
rows=[]
for name,s in sets.items():
 a=s['axis']
 for j,(x,y) in enumerate(s['points'],1):
  overlap=(name=='fig6_flow' and j>=4)
  px_bound=4.0 if overlap else 2.0
  rows.append({'study':'Alem2015','figure':s['figure'],'series':'0.25 g/L circles','mode':s['mode'],'quantity':s['quantity'],'point_id':j,'pore_volumes':round(100*(x-a['x0'])/(a['x100']-a['x0']),6),'value':round((a['y0']-y)/(a['y0']-a['y1']),6),'pixel_x':x,'pixel_y':y,'digitization_x_bound_PV':round(100*2/(a['x100']-a['x0']),6),'digitization_value_bound':round(px_bound/(a['y0']-a['y1']),6),'overlapping_series':overlap,'source_file':s['file'],'method':'manually selected experimental circle marker center','experimental_uncertainty':'not reported for individual plotted observations'})
fields=list(rows[0])
for file,subset in [('data_observations.csv',rows),('data_breakthrough.csv',[r for r in rows if r['quantity']=='C_over_C0']),('data_hydraulics.csv',[r for r in rows if r['quantity']=='K_over_K0'])]:
 with (P/file).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(subset)
# Numerical endpoints transcribed from original Table1; not digitized from a fit.
with (P/'data_retention_table1.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['study','source','mode','feed_g_L','end_pore_volumes','retained_fraction_of_injected','effluent_fraction_of_injected','reported_total_fraction','method','experimental_uncertainty'])
 w.writerow(['Alem2015','Table1','constant_flow',.25,83,.6795,.3162,.9958,'transcribed printed numerical table','individual errors unavailable; published whole mass closure0.9958'])
 w.writerow(['Alem2015','Table1','constant_head',.25,83,.7394,.1964,.9358,'transcribed printed numerical table','individual errors unavailable; published whole mass closure0.9358'])

# Explicit identifiers rather than constructing implicit image origins.
source_urls={'source_fig1.png':'https://media.springernature.com/full/springer-static/image/art%3A10.1007%2Fs11242-014-0402-8/MediaObjects/11242_2014_402_Fig1_HTML.gif','source_fig2.png':'https://media.springernature.com/full/springer-static/image/art%3A10.1007%2Fs11242-014-0402-8/MediaObjects/11242_2014_402_Fig2_HTML.gif','source_fig6.png':'https://media.springernature.com/full/springer-static/image/art%3A10.1007%2Fs11242-014-0402-8/MediaObjects/11242_2014_402_Fig6_HTML.gif'}
metadata={
 'citation':'Alem, A.; Ahfir, N.-D.; Elkawafi, A.; Wang, H.Q. Hydraulic Operating Conditions and Particle Concentration Effects on Physical Clogging of a Porous Medium. Transport in Porous Media2015,106(2),303–321. DOI10.1007/s11242-014-0402-8.',
 'doi':'10.1007/s11242-014-0402-8','retrieved_date':'2026-11-25','publisher_url':'https://link.springer.com/article/10.1007/s11242-014-0402-8','primary_author_manuscript_url':'https://www.academia.edu/52371161/Alem_et_al_2014_Transp_Porous_Med','full_pdf_downloaded':False,
 'parameters_reported':{'column_length_m':.40,'column_inner_diameter_m':.045,'initial_porosity':.37,'initial_porosity_reported_plusminus':.01,'sand_grain_range_um':[315,630],'sand_median_diameter_um':410,'kaolinite_particle_range_um':[1.7,40],'kaolinite_median_diameter_um':10.5,'particle_material_density_kg_m3':2600,'feed_mass_concentration_g_L':.25,'constant_volumetric_flow_mL_min':145,'constant_head_drop_m':1.05,'initial_hydraulic_conductivity_m_s':7e-4,'initial_hydraulic_conductivity_reported_plusminus_m_s':5e-5,'total_injected_pore_volumes':83,'reported_initial_Darcy_velocity_m_s':.00152,'reported_Reynolds_number':.62},
 'parameter_source_locator':'Primary author manuscript Section2; Table1; Figures2 and6 captions.',
 'derived_not_reported':{'cross_section_m2':math.pi*.045**2/4,'initial_Darcy_velocity_from_Q_m_s':(145e-6/60)/(math.pi*.045**2/4),'initial_pore_volume_m3':.37*.40*math.pi*.045**2/4,'feed_total_particle_volume_fraction':.25/2600},
 'unknowns':['two-bin inlet mass fractions','species-resolved breakthrough','immobile-liquid/deposit porosities','kinetic thresholds','water temperature/viscosity for exact intrinsic-permeability conversion','individual observational error bars'],
 'observation_meaning':{'C_over_C0':'Measured effluent concentration normalized by feed, using calibrated turbidity. Selected plotted experimental circles, not an interpolating/model curve.','K_over_K0':'Experimental hydraulic conductivity/permeability ratio inferred by authors from measured pressures and flow via Darcy law. Not a direct independent measurement of intrinsic permeability and not a current-model prediction.','retained_fraction_of_injected':'Post-experiment retained mass normalized by injected mass; numerical source Table1.'},
 'images':{},'axes':{n:s['axis'] for n,s in sets.items()},'point_counts':{n:len(s['points']) for n,s in sets.items()},
 'digitization':{'native_pixels':True,'method':'Manual experimental marker identification; native pixel centers; affine calibration from labelled linear axes. Filled-circle centroids refined within selected local image components.','ordinary_pixel_center_bound':2,'hydraulics_flow_overlap_pixel_y_bound':4,'bounds_are_not':'experimental error bars, confidence intervals or standard deviations','startup_omissions':'Highly crowded symbols near origin and ambiguous species overlaps excluded; 0.25g/L head concentration markers roughly6–26PV cannot be separated reliably and were not invented.','high_resolution_attempt':'Publisher full and lw1200 routes return same427px-wide raster; zoomed audit images add no information.'},
 'limitations':['Broad polydisperse kaolinite is not a controlled two-size binary feed.','Original pressure model uses pressure potential; vertical experimental heads require hydraulic-potential treatment.','Two species-specific closures are not identifiable from aggregate concentration and conductivity alone.','Table1 mass recovery is0.9958 and0.9358 in selected flow/head cases; digitization error is much smaller than these experimental mass discrepancies.','Late internal cake is reported for high-concentration runs, not used here.'],
 'initial_hydraulic_consistency_note':'Using reported Q145mL/min and diameter4.5cm givesU0=.0015195m/s. Reported conductivityK0=.0007m/s times1.05m/.4m gives.0018375m/s if1.05m is total-head difference. The paper pressure-head versus total-head reference and rounded values must be resolved; do not silently force incompatible absolute baselines.',
 'hydraulic_domain_interpretation':'Figure1 depicts inlet/outlet and intermediate piezometers. Text calls Figure6 overall permeability of the column; full40cm comparison is supported, but exact numerical tap depths/head computation are not specified in schematic. See taps_and_profile_check.md.',
 'no_fitting_performed':True
}
for file,url in source_urls.items():
 im=Image.open(P/file);metadata['images'][file]={'url':url,'width':im.width,'height':im.height,'sha256':hashlib.sha256((P/file).read_bytes()).hexdigest(),'file_format_verified':im.format}
(P/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
print(json.dumps({'counts':metadata['point_counts'],'total':len(rows),'table_rows':2}))
