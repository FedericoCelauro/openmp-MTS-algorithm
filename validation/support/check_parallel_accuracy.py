#!/usr/bin/env python3
"""Check whether parallel energy differences start as same-state force errors.

Compares serial and 1-/4-thread MTS trajectories, independently sums physical
energy every step, repeats 4-thread runs and perturbs serial inputs by one ULP.
Both trajectories also supply checkpoints for identical-state one-step tests.
Raw jobs are temporary. Summaries, a separation plot and metadata are retained.
Requires the optional analysis dependencies and standard C++ build dependencies.
"""
from pathlib import Path
import csv
import subprocess
import os
import json
import math
import shutil
import numpy as np

import argparse
import tempfile
import hashlib

def persist(result, output):
    output.mkdir(parents=True,exist_ok=True)
    (output/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for key,filename in [('summary','run_summary.csv'),('same_state_tests','same_state_checks.csv'),
                         ('trajectory_difference','trajectory_separation.csv'),('perturbations','input_perturbations.csv')]:
        rows=result[key]
        with (output/filename).open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4.5))
    separation=[x for x in result['trajectory_difference'] if x['rms_position_difference']>0]
    axes[0].semilogy([x['time'] for x in separation],[x['rms_position_difference'] for x in separation],'o-')
    axes[0].set_xlabel('Physical time')
    axes[0].set_ylabel('RMS position-component difference')
    axes[0].set_title('Serial vs one 4-thread trajectory')
    groups=[('serial','Serial'),('omp1','OpenMP 1T'),('omp4_','OpenMP 4T'),('serial_ulp_','Serial + 1 ULP')]
    for index,(prefix,label) in enumerate(groups):
        rows=[r for r in result['summary'] if r['name']==prefix or r['name'].startswith(prefix) and prefix.endswith('_')]
        values=[r['max_abs_error'] for r in rows]
        offsets=np.linspace(-.08,.08,len(values)) if len(values)>1 else [0]
        axes[1].scatter(np.array(offsets)+index,values,label=label)
    axes[1].set_xticks(range(4))
    axes[1].set_xticklabels([label for _,label in groups],rotation=15)
    axes[1].set_ylabel('Maximum |E(t) - E(0)|')
    axes[1].set_title('Finite-run maxima, dt=0.04, duration=50')
    for ax in axes:ax.grid(alpha=.25)
    fig.tight_layout();fig.savefig(output/'parallel_sensitivity.png',dpi=170);plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--eigen-include',default=os.environ.get('EIGEN_INC','/usr/include/eigen3'))
    parser.add_argument('--repeats',type=int,default=8)
    parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parents[2]/'validation/results/parallel_accuracy')
    args=parser.parse_args()
    if args.repeats<2:parser.error('Use at least two repetitions')
    if not (Path(args.eigen_include)/'Eigen/Core').is_file():parser.error('Eigen/Core not found')
    with tempfile.TemporaryDirectory(prefix='md-parallel-') as temporary:
        root=Path(__file__).resolve().parents[2]; base=Path(temporary)
        eigen=args.eigen_include
        variants=('MDMultipleTimeSteps','MDMultipleTimeStepsOpenMP')
        initial=root/'MDNormal/jobs/sim1/initial_conf.csv'
        cfg={}
        for line in (root/'MDNormal/jobs/sim1/config.txt').read_text().splitlines():
            parts=line.replace('=',' ').split()
            if len(parts)==2 and parts[0]!='delimiter':cfg[parts[0]]=parts[1]
        cfg['dt']='0.04';cfg['totsteps']='1250'

        def job(name, source=initial, steps=1250):
            p=base/name;p.mkdir(exist_ok=True)
            (p/'config.txt').write_text('\n'.join(f'{k} {steps if k=="totsteps" else v}' for k,v in cfg.items())+'\n')
            shutil.copyfile(source,p/'initial_conf.csv')
            return p

        def run(variant, threads, p, prefix):
            subprocess.run([str(base/(variant+'_trace')),str(p),str(prefix)],check=True,
                           env=dict(os.environ,OMP_NUM_THREADS=str(threads),OMP_DYNAMIC='FALSE',OMP_PROC_BIND='close',OMP_PLACES='cores',OMP_WAIT_POLICY='ACTIVE'))
            data=np.genfromtxt(str(prefix)+'.energy.csv',delimiter=',',names=True)
            error=data['E']-data['E'][0]
            result={'name':prefix.name,'variant':variant,'threads':threads,'max_abs_error':float(np.max(np.abs(error))),
                    'final_error':float(error[-1]),'max_energy_sum_discrepancy':float(np.max(np.abs(data['E']-data['E_fixed']))),
                    'max_U_sum_discrepancy':float(np.max(np.abs(data['U']-data['U_fixed']))),
                    'max_K_sum_discrepancy':float(np.max(np.abs(data['K']-data['K_fixed'])))}
            print(result,flush=True)
            return result

        for v in variants:
            subprocess.run(['make','-C',str(root/v),'-j2','BUILD_MODE=releaseO3','EIGEN_INC='+eigen],check=True,stdout=subprocess.DEVNULL)
            flags=['-std=c++17','-O3','-I'+str(root/v/'src/include'),'-I'+eigen]
            if v.endswith('OpenMP'): flags.append('-fopenmp')
            objects=sorted(p for p in (root/v/'build/releaseO3').rglob('*.o') if p.name!='MolecularDynamics.o')
            subprocess.run(['g++',*flags,str(root/'validation/support/trace_mts_diagnostics.cpp'),*map(str,objects),'-o',str(base/(v+'_trace')),'-lstdc++fs'],check=True)

        summaries=[]
        perturbations=[]
        for name,v,t in [('serial',variants[0],1),('omp1',variants[1],1)]+[(f'omp4_{i}',variants[1],4) for i in range(args.repeats)]:
            summaries.append(run(v,t,job(name),base/name))
        # ULP changes probe sensitivity of serial itself, not a finite physical change.
        with initial.open() as stream: source=list(csv.reader(stream))
        for i in range(args.repeats):
            p=job(f'serial_ulp_{i}')
            data=[r.copy() for r in source]
            particle=1+i%6;column=5+i%3
            original=float(data[particle][column])
            changed=float(np.nextafter(original, math.inf if i%2==0 else -math.inf))
            data[particle][column]=repr(changed)
            perturbations.append({'run':p.name, 'particle_id':int(data[particle][0]),
                                  'field':source[0][column], 'before':original, 'after':changed,
                                  'difference':changed-original})
            with (p/'initial_conf.csv').open('w',newline='') as stream:csv.writer(stream).writerows(data)
            summaries.append(run(variants[0],1,p,base/p.name))

        # Compare fresh identical states sampled from both trajectories, at init and after one step.
        state_tests=[]
        for trajectory in ('serial','omp4_0'):
            with Path(str(base/trajectory)+'.state.csv').open() as stream: checkpoints=list(csv.DictReader(stream))
            for step in sorted({int(x['step']) for x in checkpoints}):
                p=job('fixed_state',steps=1)
                with (p/'initial_conf.csv').open('w',newline='') as stream:
                    writer=csv.writer(stream);writer.writerow(['id','m','x','y','z','vx','vy','vz'])
                    writer.writerows([x[k] for k in ('id','m','x','y','z','vx','vy','vz')] for x in checkpoints if int(x['step'])==step)
                outputs=[]
                for v,t in ((variants[0],1),(variants[1],4)):
                    prefix=base/('same_'+v)
                    subprocess.run([str(base/(v+'_trace')),str(p),str(prefix)],check=True,env=dict(os.environ,OMP_NUM_THREADS=str(t),OMP_DYNAMIC='FALSE'))
                    state=np.genfromtxt(str(prefix)+'.state.csv',delimiter=',',names=True)
                    outputs.append(state)
                a,b=outputs
                row={'trajectory':trajectory,'checkpoint_step':step}
                for group,cols in [('initial_force',('fx','fy','fz')),('one_step_state',('x','y','z','vx','vy','vz')),('one_step_force',('fx','fy','fz','fxAux','fyAux','fzAux'))]:
                    mask=a['step']==(0 if group=='initial_force' else 1)
                    actual=np.column_stack([a[c][mask] for c in cols]);other=np.column_stack([b[c][mask] for c in cols])
                    row[group+'_max_absolute_difference']=float(np.max(np.abs(actual-other)))
                    row[group+'_relative_scale_difference']=float(np.max(np.abs(actual-other))/max(1,np.max(np.abs(actual))))
                state_tests.append(row)
        print('Same-state comparison maximums:',{k:max(x[k] for x in state_tests) for k in state_tests[0] if k.endswith('difference')},flush=True)
        # First divergence of checkpoint positions, using minimum-image position differences.
        a=np.genfromtxt(str(base/'serial')+'.state.csv',delimiter=',',names=True)
        b=np.genfromtxt(str(base/'omp4_0')+'.state.csv',delimiter=',',names=True)
        trajectory_difference=[]
        for step in sorted(set(a['step'])):
            mask=a['step']==step
            pos=np.column_stack([a[c][mask]-b[c][mask] for c in ('x','y','z')]);box=np.array([int(cfg[k]) for k in ('nx','ny','nz')])*float(cfg['rho'])**(-1/3);pos-=box*np.round(pos/box)
            row={'step':int(step),'time':float(step*.04),'rms_position_difference':float(np.sqrt(np.mean(pos**2))),'max_position_component_difference':float(np.max(np.abs(pos)))}
            trajectory_difference.append(row)
        print('Trajectory separation:',trajectory_difference,flush=True)
        serial=np.genfromtxt(str(base/'serial')+'.energy.csv',delimiter=',',names=True)
        one=np.genfromtxt(str(base/'omp1')+'.energy.csv',delimiter=',',names=True)
        result={'config':cfg,'summary':summaries,'same_state_tests':state_tests,'trajectory_difference':trajectory_difference,
                'serial_vs_omp1_energy_exact':bool(np.array_equal(serial['E'],one['E'])),
                'perturbations':perturbations,
                'metadata':{'initial_csv_sha256':hashlib.sha256(initial.read_bytes()).hexdigest(),
                            'particles':len(source)-1,'dt':0.04,'duration':50,'repeats':args.repeats,
                            'independent_energy_sum':'Serial fixed pair order with long-double accumulation',
                            'note':'No integration formulas changed. ULP changes probe sensitivity; this is not an ensemble-equivalence proof.'}}
        persist(result,args.output_dir)


if __name__=='__main__':
    main()
