"""Fit the ten original scans and reproduce all thirteen analysis figures.

The scientific calculations are unchanged from the report-generation script.
Only dependency loading and input/output paths were made portable.
Raw input files are read only. See README.md for model assumptions.
"""
from pathlib import Path
import argparse
import numpy as np, json, csv
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description="Reproduce all 13 quantum sensing report plots from the original scan exports.")
parser.add_argument("--data-dir", type=Path, default=ROOT / "raw")
parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs")
args = parser.parse_args()
D = args.data_dir.resolve()
W = args.output_dir.resolve()
if W == D:
    parser.error("The output directory must differ from the raw data directory.")
W.mkdir(parents=True, exist_ok=True)
labels=['MI reference','MI sample','HOM reference','HOM sample 2','HOM sample 1 (thinnest)','HOM sample 3','HWP 0 deg broad scan','HWP 10 deg','HWP 0 deg fine scan','HWP 23 deg']
def hom(x,b,a,m,s):return b-a*np.exp(-0.5*((x-m)/s)**2)
def mi(x,b,a,m,s,f,p):return b+a*np.exp(-0.5*((x-m)/s)**2)*np.cos(2*np.pi*f*x+p)
res=[];arrays={};excluded=[]
for i in range(1,11):
 p=D/f'linear-scan-data_{i:02}.txt';raw=np.loadtxt(p,skiprows=5)
 if raw.ndim != 2 or raw.shape[1] != 16:raise ValueError(f'{p.name}: expected 16 columns')
 valid=np.isfinite(raw).all(axis=1)&(np.abs(raw[:,0])<20)&(raw[:,1]>0)&(raw[:,2]>0)&(raw[:,5]>0)
 for j in np.where(~valid)[0]:excluded.append([i,int(j+6),*raw[j].tolist()])
 a=raw[valid];arrays[i]=a;x=a[:,0];y=a[:,5];err=np.sqrt(y)
 if i<3:
  best=None
  for f in [22.0,22.3,22.6]:
   for phase in [0,3.14]:
    try:
     par,cov=curve_fit(mi,x,y,p0=[np.median(y),np.ptp(y)/2,x[np.argmax(abs(y-np.median(y)))],.35,f,phase],sigma=err,absolute_sigma=True,bounds=([0,0,x.min(),.03,21,-100],[1e5,1e5,x.max(),5,24,100]),maxfev=1500)
     chi=np.sum(((y-mi(x,*par))/err)**2)
     if best is None or chi<best[0]:best=(chi,par,cov)
    except RuntimeError:pass
  chi,par,cov=best;fun=mi;v=par[1]/par[0];grad=np.array([-par[1]/par[0]**2,1/par[0],0,0,0,0]);vc=float(grad@cov@grad)
 else:
  fun=hom;par,cov=curve_fit(fun,x,y,p0=[np.percentile(y,80),np.ptp(y),x[y.argmin()],.008],sigma=err,absolute_sigma=True,bounds=([0,0,x.min(),.0002],[1e5,1e5,x.max(),.2]),maxfev=20000)
  chi=np.sum(((y-fun(x,*par))/err)**2);v=par[1]/(2*par[0]-par[1]);den=2*par[0]-par[1];grad=np.array([-2*par[1]/den**2,2*par[0]/den**2,0,0]);vc=float(grad@cov@grad)
 red=chi/(len(y)-len(par));scale=max(1,red);covr=cov*scale
 item=dict(scan=i,label=labels[i-1],n=len(a),excluded=int((~valid).sum()),xmin=float(x.min()),xmax=float(x.max()),step=float(np.median(abs(np.diff(x)))),params=par.tolist(),errors=np.sqrt(np.diag(covr)).tolist(),cov=covr.tolist(),chi2red=red,V=v,uV=np.sqrt(vc*scale),uV_poisson=np.sqrt(vc),fwhm=2.354820045*par[3],ufwhm=2.354820045*np.sqrt(covr[3,3]))
 res.append(item)
 fig,axs=plt.subplots(2,1,figsize=(8,5.4),gridspec_kw={'height_ratios':[3,1]},sharex=True)
 axs[0].errorbar(x,y,yerr=err,fmt='.',ms=3,lw=.5,color='#194b73',label='Measured coincidences (1 s)')
 xx=np.linspace(x.min(),x.max(),4000);axs[0].plot(xx,fun(xx,*par),color='#bd512c',lw=1.2,label='Gaussian-envelope cosine' if i<3 else 'Gaussian dip fit')
 axs[0].set_ylabel('Coincidence rate (s$^{-1}$)');axs[0].set_title(f'{i:02}  {labels[i-1]}');axs[0].legend(fontsize=8)
 axs[1].plot(x,(y-fun(x,*par))/err,'.',ms=3);axs[1].axhline(0,color='black',lw=.6);axs[1].set_ylabel('Residual / $\\sqrt{N}$');axs[1].set_xlabel('Stage position (mm)')
 fig.tight_layout();fig.savefig(W/f'fit_{i:02}.png',dpi=160);plt.close(fig)
 np.savetxt(W/f'clean_{i:02}.csv',a,delimiter=',',header='position_mm,single0_s-1,single1_s-1,ch2,ch3,coinc01_s-1,02,03,12,13,23,012,013,023,123,0123',comments='')
(W/'results.json').write_text(json.dumps(res,indent=2),encoding='utf8')
with (W/'excluded_rows.csv').open('w',newline='') as f:
 cw=csv.writer(f);cw.writerow(['scan','original_text_line','position','ch0','ch1','ch2','ch3','01','02','03','12','13','23','012','013','023','123','0123']);cw.writerows(excluded)
for r in res:print(r['scan'],r['n'],'excluded',r['excluded'],'params',np.round(r['params'],6),'V',round(r['V'],5),'+-',round(r['uV'],5),'chi2',round(r['chi2red'],2))
for name,ids in [('glass',[3,5,4,6]),('waveplate',[9,8,10])]:
 fig,ax=plt.subplots(figsize=(8,4.7))
 for i in ids:
  r=res[i-1];a=arrays[i];ax.errorbar(a[:,0],a[:,5]/r['params'][0],yerr=np.sqrt(a[:,5])/r['params'][0],fmt='.',ms=3,lw=.5,label=labels[i-1]);xx=np.linspace(a[:,0].min(),a[:,0].max(),800);ax.plot(xx,hom(xx,*r['params'])/r['params'][0],lw=1)
 ax.set(xlabel='Stage position (mm)',ylabel='Coincidence rate / fitted baseline');ax.legend(fontsize=8);fig.tight_layout();fig.savefig(W/f'{name}.png',dpi=160);plt.close(fig)
fig,ax=plt.subplots(figsize=(7,4))
rr=[res[i-1] for i in [9,8,10]];ax.errorbar([0,10,23],[r['V'] for r in rr],yerr=[r['uV'] for r in rr],fmt='o',capsize=3,label='Manual visibility')
angles=np.linspace(0,45,200);v0=rr[0]['V'];q=2*v0/(1+v0)*np.cos(np.deg2rad(2*angles))**2;ax.plot(angles,q/(2-q),'--',label='Ideal HWP trend anchored at 0 deg')
ax.set(xlabel='Recorded half-wave plate angle (deg)',ylabel='Visibility (manual convention)',ylim=(0,.9));ax.legend(fontsize=8);fig.tight_layout();fig.savefig(W/'visibility.png',dpi=160)

plt.close(fig)
