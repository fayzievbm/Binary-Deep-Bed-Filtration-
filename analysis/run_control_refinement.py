"""Refine the matched-capacity controls to distinguish quality windows from smearing."""
from pathlib import Path
from dataclasses import replace
import json,time
import numpy as np
import run_mode_controls as study
from run_mode_controls import Params

def main():
    start=time.perf_counter();cases={}
    reference=json.loads((study.OUT/'results.json').read_text())['cases']
    for mode in ['pressure','rate']:
        for kin in ['capped_linear','charging_linear','multistage']:
            name=f'control_{mode}_{kin}'
            cases[name]=[dict(n=160,summary=reference[name])]
            for n in [80,320]:
                r=study.run(name+f'_n{n}',Params(mode=mode,kinetics=kin),n,14.)
                cases[name].append(dict(n=n,summary=study.summary(r)))
            cases[name].sort(key=lambda v:v['n'])
    report=dict(scope='Fixed-parameter spatial/time-step-CFL refinement of matched-capacity controls; all runs use original t* equations and the same PV14 endpoint.',
        quality_note='Quality thresholds are illustrative, not water standards. Effluent PV excludes the first initial PV, but small residual windows near the advective front can still be dominated by numerical spreading.',
        cases=cases,elapsed_seconds=time.perf_counter()-start)
    (study.OUT/'control_refinement.json').write_text(json.dumps(study.clean_json(report),indent=2))
    print(json.dumps(study.clean_json(report),indent=2))

if __name__=='__main__':main()
