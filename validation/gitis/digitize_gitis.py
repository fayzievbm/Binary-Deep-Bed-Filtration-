"""Convert manually inspected symbol centers to data coordinates; no model curves fitted."""
from pathlib import Path
import csv, json, hashlib
import numpy as np
from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parent
# Coordinates refer to the unmodified native-resolution coauthor-posted JPGs.
FIGURES = {
 '6': {'file':'figure_006.jpg', 'x_ticks':[[303.5,0],[470.5,2],[638.5,4],[806.5,6],[974.5,8]],
       'y_ticks':[[507.5,0],[410.0,.05],[312.5,.10],[215.0,.15],[117.5,.20],[20.0,.25]],
       'points': {
        '0.82':[[303.5,195],[387,161],[429,151],[513,424],[555,444],[639,460],[723,462],[807,467],[863,468]],
        '1.02':[[303.5,347],[387,298],[429,364],[513,440],[555,471],[639,476],[723,468],[807,477],[863,475]],
        '1.62':[[303.5,438],[387,441],[429,469],[513,474],[555,495],[639,495],[723,491],[807,491],[863,485]]}},
 '7': {'file':'figure_008.jpg', 'x_ticks':[[264,0],[404.5,2],[544.5,4],[685,6],[825.5,8],[965.5,10]],
       'y_ticks':[[590.5,0],[493.5,.1],[397,.2],[300.5,.3],[204,.4],[108,.5],[11.5,.6]],
       'points': {
        '0.82':[[275,155],[287,163],[298,138],[310,150],[334,192],[345,159],[357,231],[369,227],[381,198],[392,298],[404,303],[416,319],[427,318],[439,315],[474,377],[521,422],[580,437],[685,441],[743,435],[790,435]],
        '1.02':[[275,295],[287,295],[298,307],[310,309],[334,328],[345,382],[357,389],[369,385],[381,394],[392,420],[404,437],[416,467],[427,453],[439,441],[474,465],[521,476],[580,478],[685,501],[743,507],[790,501]],
        '1.62':[[275,358],[287,371],[298,371],[310,380],[334,423],[345,447],[357,460],[369,458],[381,453],[392,497],[404,485],[416,510],[427,507],[439,509],[474,513],[521,522],[580,519],[685,536],[743,541],[790,534]]}}
}
# Manual center selection +/-2 px, tick center +/-1 px, plus JPEG/line overlap +/-1 px.
# Use a conservative combined +/-3 px positional envelope, not an experimental SD.
POSITION_ENVELOPE_PX = 3.0
MARKERS = {'0.82':'filled_circle','1.02':'open_circle','1.62':'filled_down_triangle'}
COLORS = {'0.82':(220,40,40),'1.02':(30,110,230),'1.62':(0,150,100)}
calibrations={}
all_rows=[]
for fig, spec in FIGURES.items():
    im=Image.open(ROOT/'sourcefigs'/spec['file']).convert('RGB')
    xt=np.array(spec['x_ticks']);yt=np.array(spec['y_ticks'])
    # Only axis calibration, not physical model fitting.
    ax,bx=np.polyfit(xt[:,0],xt[:,1],1); ay,by=np.polyfit(yt[:,0],yt[:,1],1)
    calibrations[fig]={
      'source_image':spec['file'],'width_px':im.width,'height_px':im.height,
      'source_sha256':hashlib.sha256((ROOT/'sourcefigs'/spec['file']).read_bytes()).hexdigest(),
      'x_ticks_pixel_value':spec['x_ticks'],'y_ticks_pixel_value':spec['y_ticks'],
      'time_h_from_pixel':{'slope':float(ax),'intercept':float(bx)},
      'concentration_ratio_from_pixel':{'slope':float(ay),'intercept':float(by)},
      'tick_fit_max_abs_residual_h':float(np.max(np.abs(ax*xt[:,0]+bx-xt[:,1]))),
      'tick_fit_max_abs_residual_ratio':float(np.max(np.abs(ay*yt[:,0]+by-yt[:,1]))),
      'positional_uncertainty_envelope_px':POSITION_ENVELOPE_PX,
      'time_uncertainty_envelope_h':float(abs(ax)*POSITION_ENVELOPE_PX),
      'ratio_uncertainty_envelope':float(abs(ay)*POSITION_ENVELOPE_PX),
      'uncertainty_type':'conservative digitization envelope; not experimental uncertainty, SD or confidence interval'
    }
    rows=[]
    overlay=im.resize((im.width*2,im.height*2))
    draw=ImageDraw.Draw(overlay)
    for depth,points in spec['points'].items():
      for i,(x,y) in enumerate(points,1):
        t=float(ax*x+bx)
        if abs(t)<0.01:t=0.0 # The three Fig6 first symbols align with the t=0 tick.
        row={'figure':int(fig),'series':f'fig{fig}_depth_{depth}m','point_id':f'F{fig}_D{depth}_{i:02d}',
             'depth_m':float(depth),'time_h':round(t,5),'concentration_ratio':round(float(ay*y+by),6),
             'pixel_x':x,'pixel_y':y,'marker':MARKERS[depth],
             'time_digitization_envelope_h':round(abs(ax)*POSITION_ENVELOPE_PX,5),
             'ratio_digitization_envelope':round(abs(ay)*POSITION_ENVELOPE_PX,6),
             'data_origin':'digitized_experimental_symbol','raw_measurement':False,
             'curve_sampled':False,'quality_flag':'touching_or_near_curve_symbols; manual center envelope covers visual ambiguity' if fig=='7' and 6<=i<=14 else 'manual_symbol_center'}
        rows.append(row)
        xx,yy=2*x,2*y
        draw.ellipse((xx-12,yy-12,xx+12,yy+12),outline=COLORS[depth],width=2)
        draw.line((xx-5,yy,xx+5,yy),fill=COLORS[depth],width=2)
        draw.line((xx,yy-5,xx,yy+5),fill=COLORS[depth],width=2)
    for pixel,value in spec['x_ticks']:
        xx=2*pixel;yy=(595 if fig=='6' else 590.5)*2
        draw.line((xx,yy-12,xx,yy+12),fill=(170,0,200),width=2)
    with (ROOT/f'figure{fig}_experimental_points.csv').open('w',newline='') as f:
      writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    overlay.save(ROOT/'sourcefigs'/f'figure{fig}_symbol_selection_QA.png')
    all_rows+=rows
with (ROOT/'all_experimental_points.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(all_rows[0]));writer.writeheader();writer.writerows(all_rows)
(ROOT/'axis_calibration.json').write_text(json.dumps(calibrations,indent=2)+'\n')
(ROOT/'manual_symbol_centers.json').write_text(json.dumps(FIGURES,indent=2)+'\n')
print(json.dumps({'figure6_points':len([r for r in all_rows if r['figure']==6]),'figure7_points':len([r for r in all_rows if r['figure']==7]),'total':len(all_rows),'calibration':calibrations},indent=2))
