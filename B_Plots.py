import numpy as n
L = n.genfromtxt
from numpy import nanmax, nanmin, abs
import matplotlib.pyplot as p
from matplotlib.collections import LineCollection
from matplotlib import cm
import figfun as f
import os
from sys import argv
p.close('all')

try:    
    worthless, savefigs, expt = argv
    savefigs = bool(save)
    expt = int(expt)
except:
    savefigs = False
    expt = 1
    os.chdir('../BT-{}_FS40SS10'.format(expt))

to = .0474    
    
#########
# Results.dat
# '[0]BulgeHeight, [1]Sphere Rad [2]Radius-Rolling [3]Radius-Transv. [4]MajorLogStn [5]MinorLogStn \n[6]ex [7]ey [8]Gamma  [9]LEeq [10]Thickness (incomp.) [11]Thickness (iterated) [12]Sdev MajStn [13]TotalPts [14]Filtered Pts'  
STLP = L('STLP.dat',delimiter=',')
d = L('Results.dat',delimiter=',')
profStgs = L('profStgs.dat',delimiter=',').astype(int)
profs = L('Profiles_Filtered.dat',delimiter=',')

bbox_props = dict(boxstyle="circle", fc="w", ec="r", alpha=1)

###################################################    
###### FIGURE 1 - Quad Plot ###########
###################################################    
p.style.use('mysty-quad')

fig1, ax111, ax112, ax121, ax122 = f.makequad()
###### ax111 - Top Left - Pressure vs h/r ###########
p.sca(ax111)
p.plot(d[:,0]/3,STLP[:,-1],'b')
p.plot(d[profStgs,0]/3,STLP[profStgs,-1],'o',mec='r',mfc='r',ms=4)

for k,i in enumerate(profStgs):
    xloc = d[i,0]/3+.02
    yloc = -80+STLP[i,-1]
    if i == profStgs[-2]:
        xloc = d[i,0]/3-.01

    p.text(xloc,yloc,'{}'.format(k+1),color='r',
            ha='center',va='top',size=10,
            bbox=bbox_props)
    p.annotate('',xy=[d[i,0]/3,STLP[i,-1]],xytext=[xloc,yloc],
                       color='r',ha='center',va='top',
                       arrowprops=dict(arrowstyle='-', color='r'),transform=ax111.transData)
f.eztext(ax111,'Al-6022-T43','ul')
f.eztext(ax111,'BT-{}'.format(expt),'lr')
p.xlabel('h/R')
p.ylabel('P\n(psi)')
p.xticks(n.arange(0,2,.2))
ax111.set_xlim([0,1.1*n.max(d[:,0]/3)])

###### ax112 - Top Right - Pressure vs e ###########
p.sca(ax112)
ex, = p.plot(d[:,6],STLP[:,-1],'b',label='$\\mathsf{e}_{\\mathsf{x}\\prime}$')
ey, = p.plot(d[:,7],STLP[:,-1],'r',label='$\\mathsf{e}_{\\mathsf{y}\\prime}$')
f.eztext(ax112,'Al-6022-T43','ul')
f.eztext(ax112,'BT-{}'.format(expt),'lr')
p.xlabel('e')
p.ylabel('P\n(psi)')
p.legend(loc='right',frameon=False,fontsize=14)
ax112.set_xlim(left=0)

###### ax121 - Bottom Left - rho/R vs h/R ###########
p.sca(ax121)
rx, = p.plot(d[:,0]/3, d[:,2],'b',label='$\\rho_{\\mathsf{x}\\prime}$')
ry, = p.plot(d[:,0]/3, d[:,3],'r',label='$\\rho_{\\mathsf{y}\\prime}$')
f.eztext(ax121,'Al-6022-T43','ll')
f.eztext(ax121,'BT-{}'.format(expt),'ur')
p.xlabel('h/R')
p.ylabel('$\\frac{\\rho}{\\mathsf{R}}$')
p.legend(loc='right',frameon=False,fontsize=14)
p.xticks(n.arange(0,2,.2))
ax121.set_xlim([0,1.1*n.max(d[:,0]/3)])

###### ax122 - Bottom Right - P vs dt/t ###########
p.sca(ax122)
dtto = (to - d[:,[10,11]])/to
v1, = p.plot(dtto[:,0],STLP[:,-1],'r',label='$\\nu = \\mathsf{0}.\\mathsf{5}$')
v2, = p.plot(dtto[:,1],STLP[:,-1],'b',label='$\\nu = \\mathsf{0}.\\mathsf{326}$')
f.eztext(ax122,'Al-6022-T43','ul')
f.eztext(ax122,'BT-{}'.format(expt),'lr')
p.xlabel('$\\delta\\mathsf{t}/\\mathsf{t}_{\\mathsf{o}}$',fontsize=18)
p.ylabel('P\n(psi)')
p.legend(loc='right',frameon=False,fontsize=14)
ax122.set_xlim(left=0)
f.myax(ax111)
f.myax(ax112,lambda x: f.ksi2Mpa(x/1000)*10,'P\n(bar)')
f.myax(ax121)
f.myax(ax122,lambda x: f.ksi2Mpa(x/1000)*10,'P\n(bar)')
p.xticks( n.arange(0,1,.08) )
if savefigs:
    p.savefig('Fig1 - Responses.png',bbox_inches='tight')

   
###################################################    
###### FIGURE 2 - COLORED DOME PROFILES ###########
###################################################    

p.style.use('mysty')
fig2 = p.figure(facecolor='w',figsize=(14,6) )    
ax2 = fig2.add_axes([.1,.12,.9,.68])
cmap = 'viridis'
for i in range(len(profStgs)):
    data = profs[:,3*i:3*i+3]
    data = data[~n.any(n.isnan(data),axis=1),:]
    x,z,e = (data * [1/3,1/3,1]).T
    if i == 0:
        lobound = n.min(z)
    elif i == len(profStgs) - 1:
        upbound = n.min(z)
    else:
        pass
    xneg = -x
    points = n.array([x, z]).T.reshape(-1, 1, 2)
    segments = n.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap=cmap,lw=3,array=e,norm=p.Normalize(.1,n.max(d[:,9])))
    ax2.add_collection(lc)
    points = n.array([xneg, z]).T.reshape(-1, 1, 2)
    segments = n.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap=cmap,lw=3,array=e,norm=p.Normalize(.1,n.max(d[:,9])))
    ax2.add_collection(lc)
ax2.axis(ymin=0,xmin=-1*n.max(x)*1.1,xmax=n.max(x)*1.1)
ax2.axis('equal')
ypos = n.linspace(lobound,upbound,len(profStgs))
for i in range(len(profStgs)):    
    p.text(.98 , ypos[i]/ax2.get_ylim()[1], 
            '{}'.format(i+1), color='r',bbox=bbox_props,size=10,va='center',ha='right',transform=ax2.transAxes)
ax2.set_xticks(n.arange(-.8,1,.2))            
cbar = p.colorbar(lc,ticks=n.arange(.1,.7,.1))
cbar.set_label('e$_\\mathsf{e}$', rotation=0)
p.xlabel('r/R')
p.ylabel('$\\frac{\\mathsf{z}\\prime}{\\mathsf{R}}$')
f.eztext(ax2,'Al-6022-T43','ul')
f.eztext(ax2,'BT-{}'.format(expt),'tr')
f.myax(fig2,TW=.0025,HW=.25,HL=.045,OH=.2)
#f.colorbar(ax2,cbar,TW=.0025,HW=.25,HL=.045,OH=.2)
if savefigs:
    p.savefig('Fig2 - Colored Profiles.png',bbox_inches='tight')

###################################################    
###### FIGURE 3 - Equiv. Stn vs r/R ###########
################################################### 
fig3 = p.figure(3)
ax3 = fig3.add_subplot(111)
for i in range(len(profStgs)):
    p.plot(profs[:,3*i]/3,profs[:,3*i+2],'b')
    if i == 0:
        lobound = profs[n.nanargmax(profs[:,3*i]),3*i+2]
    elif i == len(profStgs) - 1:
        upbound = profs[n.nanargmax(profs[:,3*i]),3*i+2]
    else:
        pass

ax3.set_xlim([ 0,1.1*n.nanmax(profs[:,3*i])/3 ])
ax3.set_ylim([0, 1.1*n.nanmax(profs[:,-1]) ])
ypos = n.linspace(lobound,upbound,len(profStgs))
for i in range(len(profStgs)):    
    p.text(1.05*n.nanmax(profs[:,-3])/3 , ypos[i], 
            '{}'.format(i+1), color='r',bbox=bbox_props,size=8)
p.xlabel('r/R')
p.ylabel('$\\mathsf{e}_\\mathsf{e}$')
f.eztext(ax3,'Al-6022-T43','ul')
f.eztext(ax3,'BT-{}'.format(expt),'ur')
f.myax(fig3)
if savefigs:
    fig3.savefig('Fig3 - Stn Profile.png',bbox_inches='tight')

fig4 = p.figure(4)
ax4 = fig4.add_subplot(111)
ax4.plot(STLP[1:,1],STLP[1:,-1]-60,'b')
f.eztext(ax4,'BT-{}'.format(expt),'ul')
f.eztext(ax4,'$\\mathsf{{P}}_{{\\mathsf{{max}}}}$ = {:.0f}'.format(n.max(STLP[:,-1])),'br')
ax4.set_xlabel('Time (s)')
ax4.set_ylabel('P\n(psi)')
f.myax(ax4)
if savefigs:
        fig4.savefig('Fig4 - Pressure-Time.png',bbox_inches='tight')

fig5 = p.figure(5)
ax5 = fig5.add_subplot(111)
ax5.plot(d[:,0]/3,STLP[:,-1],'b')
f.eztext(ax5,'BT-{}'.format(expt),'ul')
f.eztext(ax5,'$\\mathsf{{P}}_{{\\mathsf{{max}}}}$ = {:.0f}'.format(n.max(STLP[:,-1])),'br')
ax5.set_xlabel('h/R')
ax5.set_ylabel('P\n(psi)')
f.myax(ax5)
if savefigs:
        fig5.savefig('Fig5 - Pressure-height.png',bbox_inches='tight')
        
os.chdir('../PyFiles')



