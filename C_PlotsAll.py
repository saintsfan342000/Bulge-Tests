import numpy as n
L = n.genfromtxt
from numpy import nanmax, nanmin, abs
import matplotlib.pyplot as p
from matplotlib.collections import LineCollection
from matplotlib import cm
import figfun as f
import os


expts = n.array([1,2])
to = .0474    
p.style.use('mysty-quad')
from cycler import cycler
p.rcParams['axes.prop_cycle'] = cycler('color',['b','r','g'])
p.rcParams['lines.linewidth'] = 1.5
alpha = 0.5

for K,ex in enumerate(expts):

    #########
    # Results.dat
    # '[0]BulgeHeight, [1]Sphere Rad [2]Radius-Rolling [3]Radius-Transv. [4]MajorLogStn [5]MinorLogStn \n[6]ex [7]ey [8]Gamma  [9]LEeq [10]Thickness (incomp.) [11]Thickness (iterated) [12]Sdev MajStn [13]TotalPts [14]Filtered Pts'  
    STLP = L('../BT-{}_FS40SS10/STLP.dat'.format(ex), delimiter=',')
    d = L('../BT-{}_FS40SS10/Results.dat'.format(ex),delimiter=',')
    stgs = L('../BT-{}_FS40SS10/profStgs.dat'.format(ex),delimiter=',',dtype=int)
    
    ###################################################    
    ###### FIGURE 1 - Quad Plot ###########
    ###################################################    
    if K == 0:
        fig1, ax111, ax112, ax121, ax122 = f.makequad()

    ###### ax111 - Top Left - Pressure vs h/r ###########
    pline, = ax111.plot(d[:,0]/3,STLP[:,-1],label='BT-{}'.format(ex),alpha=alpha)
    if K == len(expts) - 1:
        f.eztext(ax111,'Al-6022-T43','ul',fontsize=12)
        ax111.set_xlabel('h/R')
        ax111.set_ylabel('P\n(psi)')
        ax111.set_xlim([0,1.2*n.max(d[:,0]/3)])
        leg111 = ax111.legend(loc='right')
        [tx.set_color(leg111.get_lines()[k].get_color()) for (k,tx) in enumerate(leg111.get_texts())]
        p.setp(leg111.get_lines(),lw=0)
        f.myax(ax111,nudge=('left',.2,0))

    ###### ax112 - Top Right - eeq vs h/R###########
    ax112.plot(d[:,0]/3,d[:,9],label='BT-{}'.format(ex),color=pline.get_color(),alpha=alpha)
    if ex == expts[-1]:
        f.eztext(ax112,'Al-6022-T43','ul',fontsize=12)
        ax112.set_ylabel('e$_{\\mathsf{e}}$')
        ax112.set_xlabel('h/R')
        leg112 = ax112.legend(loc='right')
        [tx.set_color(leg112.get_lines()[k].get_color()) for (k,tx) in enumerate(leg112.get_texts())]
        p.setp(leg112.get_lines(),lw=0)
        f.myax(ax112)

    ###### ax121 - Bottom Left - rho/R vs h/R ###########
    ax121.plot(d[:,0]/3, d[:,2]/3,label='BT-{}'.format(ex),color=pline.get_color(),alpha=alpha)
    if ex == expts[-1]:
        ax121.set_xlim([0,1.1*n.max(d[:,0]/3)])
        #ax121.set_xticks(n.arange(0,2,.2))
        f.eztext(ax121,'Al-6022-T43','ll',fontsize=12)
        ax121.set_xlabel('h/R')
        ax121.set_ylabel('$\\frac{\\rho}{\\mathsf{R}}$')
        leg121 = ax121.legend(loc='right')
        [tx.set_color(leg121.get_lines()[k].get_color()) for (k,tx) in enumerate(leg121.get_texts())]
        p.setp(leg121.get_lines(),lw=0)
        f.myax(ax121)
    
    ###### ax121 - Bottom Left - ey' vs ex' ###########
    m,b = n.polyfit(d[stgs[0]:stgs[2]+1,6],d[stgs[0]:stgs[2]+1,7],1)
    ax122.plot(d[:,6],d[:,7 ],
                label='BT-{}\ne$_{{\\mathsf{{y}}\\prime}}$/e$_{{\\mathsf{{x}}\\prime}} = ${:.3f}'.format(ex,m),
                color=pline.get_color(),alpha=alpha)
    x = n.array([d[stgs[0],6], d[stgs[2],6]])
    y = m*x+b
    ax122.plot(x,y,color='k',lw=1,zorder=50+K)
    if ex == expts[-1]:
        f.eztext(ax122,'Al-6022-T43','ul',fontsize=12)
        ax122.set_xlabel('e$_{\\mathsf{x}\\prime}$')
        ax122.set_ylabel('e$_{\\mathsf{y}\\prime}$')
        leg122 = ax122.legend(loc='lower right')
        #ax122.set_xticks(n.arange(0,1,
        [tx.set_color(leg122.get_lines()[k].get_color()) for (k,tx) in enumerate(leg122.get_texts())]
        p.setp(leg122.get_lines(),lw=0)
        f.myax(ax122)
        xlim = ax122.get_xlim()
        ax122.set_xticks(n.arange(0,.7,.07))
        ax122.set_yticks(n.arange(0,.7,.07))
        ax122.set_xlim(xlim)
        ax122.set_ylim(xlim)
        
fig1.savefig('../All-1.png',boox_inches='tight',dpi=150)
p.show('all')