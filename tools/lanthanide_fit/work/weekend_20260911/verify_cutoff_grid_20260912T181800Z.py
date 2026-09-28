"""Integrate actual NoSpherA2 exports across their documented screening boundaries."""
import os
os.environ.update(OMP_NUM_THREADS='14',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
import ast,datetime as dt,hashlib,json,subprocess,sys,time
from pathlib import Path
import numpy as np
from scipy.io import FortranFile
b=Path(__file__).resolve().parent;p=b/'fleet_candidate_20260912T152903Z'
sys.path.insert(0,str(b.parent.parent))
from ptb_lnf.denmat import read_denmat,evaluate_density
out=p/'cutoff_grid_20260912T181800Z';out.mkdir(exist_ok=False)
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
ptb=Path('D:/git/ptb/build/ptb_windows_ifx_localized_final_20260907.exe');nos=Path('D:/git/NoSpherA2/build/release-windows/bin/NoSpherA2.exe')
hashes={str(f):sha(f) for f in (ptb,nos,p/'weekend_candidate.atompara',p/'basis_vDZP',Path(__file__))}
original=json.loads((p/'physical_validation_manifest.json').read_text())
assert all(hashes[str(f)]==original['hashes'][str(f)] for f in (ptb,nos,p/'weekend_candidate.atompara',p/'basis_vDZP'))
def save(name,data):
 q=out/name;t=q.with_suffix(q.suffix+'.new');t.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');t.replace(q)
def boundaries(path):
 # read_ptb at wfn_class.cpp:9246; occupied primitive coefficient maximum at 1662.
 with FortranFile(path,'r') as f:
  assert f.read_ints(np.int32)[0]==1
  ncent,nbf,nmo,nprim=map(int,f.read_ints(np.int32));assert ncent==1
  f.read_record(np.uint8)
  for k in range(3):assert f.read_reals(np.float64)[0]==0
  f.read_ints(np.int32)
  types=np.array([f.read_ints(np.int32)[0] for _ in range(nprim)])
  centers=np.array([f.read_ints(np.int32)[0] for _ in range(nprim)]);assert np.all(centers==1)
  ipao=np.array([f.read_ints(np.int32)[0] for _ in range(nprim)])
  exps=f.read_reals(np.float64);contr=f.read_reals(np.float64);occ=f.read_reals(np.float64);f.read_reals(np.float64)
  mo=f.read_reals(np.float64).reshape(nmo,nbf)
 coeff=mo[:,ipao-1]*contr
 maximum=float(np.max(np.abs(coeff[occ!=0])))
 cutoff=float(np.log(5e-5/maximum));assert cutoff<0
 radii=np.sqrt(-cutoff/exps)
 edges=np.unique(np.r_[0,.125,.5,2,6,20,80,radii[(radii>0)&(radii<80)]])
 return edges,dict(maximum_occupied_primitive_coefficient=maximum,exponent_cutoff=cutoff,edges=edges.tolist(),export_occupation_sum=float(occ.sum()))
def grid(edges,order):
 u,wu=np.polynomial.legendre.leggauss(8);phi=np.arange(16)*np.pi/8
 dirs=np.array([[np.sqrt(1-t*t)*np.cos(a),np.sqrt(1-t*t)*np.sin(a),t] for t in u for a in phi]);aw=np.repeat(wu,16)*np.pi/8
 x,w=np.polynomial.legendre.leggauss(order);rr=[];ww=[]
 for lo,hi in zip(edges[:-1],edges[1:]):
  r=lo+(x+1)*(hi-lo)/2;rr.extend(r);ww.extend(w*(hi-lo)/2*r*r)
 rr=np.array(rr);coords=(rr[:,None,None]*dirs).reshape(-1,3);weights=(np.array(ww)[:,None]*aw).ravel();radii=np.repeat(rr,len(aw))
 for limit in (20,80):assert abs(weights[radii<limit].sum()/(4*np.pi*limit**3/3)-1)<1e-12
 return coords,weights,radii
def metrics(rho,w,r):
 assert np.isfinite(rho).all() and rho.min()>=-1e-9
 counts={str(limit):float(np.dot(rho[r<limit],w[r<limit])) for limit in (20,80)}
 return dict(counts=counts,rms_radius=float(np.sqrt(np.dot(rho*w,r*r)/counts['80'])))
def call(cmd,cwd,log):
 with log.open('w') as f:subprocess.run(cmd,cwd=cwd,stdout=f,stderr=subprocess.STDOUT,timeout=90,check=True,creationflags=subprocess.CREATE_NO_WINDOW)
save('manifest.json',dict(utc=dt.datetime.now(dt.timezone.utc).isoformat(),hashes=hashes,method='NoSpherA2 unchanged CPU executable; 32/64 radial Gauss orders on segments split at every primitive screening radius, angular8x16 exact for atom-centred polynomial products through f. Fixed 1e-5 count convergence and 1e-3 counts. No scoring or physical limits changed.'))
results=[];start=time.monotonic()
for e,z in [('La',57),('Sm',62),('Tb',65),('Eu',63),('Er',68),('Tm',69)]:
 entry=dict(element=e,states=[],status='checking')
 try:
  for charge in (0,2,3):
   expected=z-46-charge;models={}
   for model in ('anchor','candidate'):
    assert time.monotonic()-start<600 and sha(nos)==hashes[str(nos)] and sha(ptb)==hashes[str(ptb)]
    work=p/'atomic'/e/f'charge{charge}'/model;work.mkdir(parents=True,exist_ok=True)
    if not (work/'ptb.denmat').exists():
     xyz=work/'atom.xyz';xyz.write_text(f'1\nAtomic count/size check\n{e} 0 0 0\n')
     params=p/'anchors'/(e+'.atompara') if model=='anchor' else p/'weekend_candidate.atompara'
     call([str(ptb),str(xyz),'-par',str(params),'-bas',str(p/'basis_vDZP'),'-chrg',str(charge),'-uhf',str(expected%2),'-stda','-denmat',str(work/'ptb.denmat')],work,work/'ptb.log')
    dm=read_denmat(work/'ptb.denmat');trace=float(np.trace(dm.density@dm.overlap));assert abs(trace-expected)<1e-6
    edges,screen=boundaries(work/'wfn.xtb');assert abs(screen['export_occupation_sum']-expected)<1e-6
    dest=out/e/f'charge{charge}'/model;dest.mkdir(parents=True)
    runs=[]
    for order in (32,64):
     coords,w,r=grid(edges,order);points=dest/f'grid{order}.points'
     with points.open('w') as f:f.write(str(len(coords))+'\n');np.savetxt(f,coords,fmt='%.16e')
     call([str(nos),'-no_gpu','-cpus','14','-wfn',str(work/'wfn.xtb'),'-rho_at_points',str(points)],dest,dest/f'nosphera{order}.log')
     rho=np.loadtxt(str(points)+'.rho');assert rho.shape==w.shape
     v=metrics(rho,w,r);assert all(abs(c-expected)<1e-3 for c in v['counts'].values())
     runs.append(dict(order=order,**v))
     if order==32 and e in ('Eu','Er','Tm'):
      physical=metrics(evaluate_density(dm,coords),w,r)
    delta=max(abs(runs[0]['counts'][k]-runs[1]['counts'][k]) for k in ('20','80'))
    models[model]=dict(trace=trace,screening=screen,runs=runs,convergence=delta,export_sha256=sha(work/'wfn.xtb'))
    if e in ('Eu','Er','Tm'):models[model]['ptb_denmat']=physical
    assert delta<1e-5,(e,charge,model,'convergence',delta)
   radius=models['candidate']['runs'][-1]['rms_radius']/models['anchor']['runs'][-1]['rms_radius']
   state=dict(charge=charge,radius_ratio=radius,models=models)
   if e in ('Eu','Er','Tm'):state['ptb_radius_ratio']=models['candidate']['ptb_denmat']['rms_radius']/models['anchor']['ptb_denmat']['rms_radius']
   entry['states'].append(state)
   assert .9<=radius<=1.1,(e,charge,'atomic_size',radius)
  entry['status']='passed'
 except Exception as exc:entry.update(status='failed',error=repr(exc))
 results.append(entry);save('result.json',dict(results=results,elapsed_seconds=time.monotonic()-start,release=False));print(e,entry['status'],entry.get('error',''),flush=True)
