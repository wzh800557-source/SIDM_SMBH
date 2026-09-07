from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import argparse
p=argparse.ArgumentParser(description="Regenerate the manuscript domain and scan diagrams from measured tables.")
p.add_argument('--results-csv',type=Path,required=True)
p.add_argument('--outdir',type=Path,required=True)
args=p.parse_args()
ROOT=args.outdir;ROOT.mkdir(parents=True,exist_ok=True)
D=list(csv.DictReader(args.results_csv.open()))
if len(D)!=1200: raise ValueError('The manuscript heat maps require 1200 measured rows')
H=sorted({float(r['represented_halo_mass_msun']) for r in D});B=sorted({float(r['black_hole_mass_msun']) for r in D});S=[10,100,1000]
H=np.array(H);B=np.array(B);hi={v:i for i,v in enumerate(H)};bi={v:i for i,v in enumerate(B)}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'axes.linewidth':.65,'pdf.fonttype':42,'savefig.dpi':300})
def edges(a):
 l=np.log10(a);return 10**np.r_[l[0]-(l[1]-l[0])/2,(l[:-1]+l[1:])/2,l[-1]+(l[-1]-l[-2])/2]
HE,BE=edges(H),edges(B)
def mat(s,key):
 z=np.full((20,20),np.nan)
 for r in D:
  if float(r['sigma0_over_m_cm2_g'])==s:z[bi[float(r['black_hole_mass_msun'])],hi[float(r['represented_halo_mass_msun'])]]=float(r[key])
 return z
def failed(s,both=False):
 z=np.zeros((20,20))
 for r in D:
  if float(r['sigma0_over_m_cm2_g'])==s:
   good=r['fp_geometric_admissible']=='True' and (r['bridge_status']=='OK' if both else True)
   z[bi[float(r['black_hole_mass_msun'])],hi[float(r['represented_halo_mass_msun'])]]=not good
 return z
def setup(ax,j,bottom=True):
 ax.set_xscale('log');ax.set_yscale('log');ax.set_xlim(HE[0],HE[-1]);ax.set_ylim(BE[0],BE[-1]);ax.tick_params(which='both',direction='out',length=2.5,labelleft=j==0,labelbottom=bottom)
 if j==0:ax.set_ylabel(r'$M_\bullet\ [\mathrm{M}_\odot]$')
def hatch(ax,z):
 if np.any(z):
  # Pad to outer cell edges so the masking covers the complete sampled plane.
  xx=np.r_[HE[0],H,HE[-1]];yy=np.r_[BE[0],B,BE[-1]];zz=np.pad(z,1,mode='edge')
  cs=ax.contourf(xx,yy,zz,levels=[.5,1.5],colors='none',hatches=['////'])
  cs.set_edgecolor('#bbbbbb');cs.set_linewidth(.3)
def save(fig,name):
 fig.savefig(ROOT/(name+'.pdf'),bbox_inches='tight',pad_inches=.04);fig.savefig(ROOT/(name+'.png'),bbox_inches='tight',pad_inches=.04);plt.close(fig)
fig=plt.figure(figsize=(7.2,2.75))
norm=LogNorm(min(float(r['r_in_over_r_h']) for r in D),max(float(r['r_in_over_r_h']) for r in D))
for j,s in enumerate(S):
 ax=fig.add_axes([.075+j*.269,.24,.228,.63]);z=mat(s,'r_in_over_r_h')
 im=ax.pcolormesh(HE,BE,z,norm=norm,cmap='viridis',rasterized=True);hatch(ax,failed(s))
 if z.min()<.5<z.max():ax.contour(H,B,z,levels=[.5],colors='white',linewidths=1.1)
 ax.set_title(rf'$\sigma_0/m={s:g}\ \mathrm{{cm^2\,g^{{-1}}}}$',fontsize=8.5,pad=5);setup(ax,j)
cax=fig.add_axes([.891,.24,.014,.63]);cb=fig.colorbar(im,cax=cax);cb.set_label(r'$r_{\rm in}/r_h$',labelpad=4)
fig.text(.45,.04,r'Represented halo mass, $M_h\ [\mathrm{M}_\odot]$',ha='center');save(fig,'fig_scan20_interface_revised')
conv=1.98847e33*1e10/(1e6*365.25*86400)
fields=[('first_shell_gravothermal_luminosity_msun_kms2_per_myr',r'$|L(r_0)|\ [\mathrm{erg\,s^{-1}}]$'),('maximum_absolute_gravothermal_luminosity_msun_kms2_per_myr',r'$\max_r|L(r)|\ [\mathrm{erg\,s^{-1}}]$')]
fig=plt.figure(figsize=(7.2,4.6))
for i,(field,label) in enumerate(fields):
 vals=np.array([abs(float(r[field]))*conv for r in D]);norm=LogNorm(vals.min(),vals.max());y=.57 if i==0 else .15
 for j,s in enumerate(S):
  ax=fig.add_axes([.075+j*.269,y,.228,.325]);im=ax.pcolormesh(HE,BE,abs(mat(s,field))*conv,norm=norm,cmap='magma',rasterized=True);hatch(ax,failed(s,True));setup(ax,j,bottom=i==1)
  if i==0:ax.set_title(rf'$\sigma_0/m={s:g}\ \mathrm{{cm^2\,g^{{-1}}}}$',fontsize=8.5,pad=5)
 cax=fig.add_axes([.891,y,.014,.325]);cb=fig.colorbar(im,cax=cax);cb.set_label(label,labelpad=4)
fig.text(.45,.015,r'Represented halo mass, $M_h\ [\mathrm{M}_\odot]$',ha='center');save(fig,'fig_scan20_luminosity_revised')
fig,ax=plt.subplots(figsize=(7.0,2.3));fig.subplots_adjust(left=.08,right=.98,bottom=.23,top=.8)
x=np.linspace(-2.8,2.8,500);N=4*np.exp(-abs(x));cross=np.log(4)
ax.axvspan(x.min(),-cross,color='#DDECF2');ax.axvspan(-cross,cross,color='#F7E8CB');ax.axvspan(cross,x.max(),color='#E6EADF')
ax.plot(x,N,color='#333333',lw=1.8);ax.axhline(1,color='#777777',lw=.8,ls='--')
for xx in (-cross,cross):ax.axvline(xx,ls=':',color='#555555',lw=1)
ax.text(-2.1,3.9,'Inner orbit domain\nFP solver',ha='center',va='bottom',fontsize=9)
ax.text(0,4.25,'Collision-dominated region\nFluid solver',ha='center',va='bottom',fontsize=9)
ax.text(2.1,3.9,'Outer orbit domain\nLMFP fluid closure',ha='center',va='bottom',fontsize=9)
ax.set_xticks([-cross,cross],[r'$r_{\rm in}$',r'$r_{\rm out}$']);ax.set_yticks([1],[r'$1$']);ax.set_ylim(0,4.3);ax.set_xlim(x.min(),x.max());ax.set_xlabel('Radius (schematic logarithmic coordinate)');ax.set_ylabel('Orbital transport depth')
ax.spines[['top','right']].set_visible(False);save(fig,'fig_domains_revised')
print('Wrote 3 vector figures and PNG previews')
