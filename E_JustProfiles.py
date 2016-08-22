import numpy as n
from numpy import (any, abs, nanmean, nanstd, isnan, nanmax, nanmin,
                    array, vstack, hstack,polyfit,exp,linspace)
from numpy.linalg import eigvalsh, eigh
import matplotlib.pyplot as p
from pandas import read_csv
from scipy.interpolate import interp1d, griddata, LinearNDInterpolator
from scipy.spatial.distance import cdist
from scipy.signal import savgol_filter as SG
import numexpr as ne
import os, glob
from mysqrtm import mysqrtm 
from CircleFitByPratt import CircleFitByPratt
from SphereFit import SphereFit
from TrueSts import TrueSts

'''
Just the profile portion of the A_DICAnalysis file.
'''

expt = 2
FS = 40
SS = 10
savedata = True
plotdata = True
makecontours = False
BIN = True
saveAram  = False

savepath = '../BT-{}_FS{}SS{}'.format(expt,FS,SS)
if BIN:
    arampath = '{}/AramisBinary'.format(savepath)
    aramprefix = 'BT-{}_'.format(expt)
else:
    arampath = '{}/AllPts'.format(savepath)
    aramprefix = 'BT-{}-Stage-0-'.format(expt)
    
os.chdir(savepath)
last = len(glob.glob('{}/{}*'.format(arampath,aramprefix))) - 1

STLP = n.genfromtxt('STLP.dat',delimiter=',')
  
##############################################
### Generate what we need for meridian strains
##############################################

# Want 4 stages before limit load, one at, and two after
LL = n.argmax(STLP[:,-1])
h = n.genfromtxt('Results.dat',delimiter=',')[:,0]
hstgs = n.hstack(( n.linspace(h[0],h[LL],7)[1:-1], n.linspace(h[LL],h[-1],2) ))
profStgs = n.zeros(len(hstgs))
for k,i in enumerate(hstgs):
    profStgs[k] = n.nonzero(h>=i)[0][0]
profStgs = profStgs.astype(int)

interpnumber = 200
profs = n.empty((interpnumber,0))
profs_filt = n.empty_like(profs)
k=0
# Iterate backwards b/c we need it to find out meridian line
for i in n.flipud(profStgs):
    if BIN:
        A = n.load('./AramisBinary/{}{}.npy'.format(aramprefix,i))
    else:
        A = read_csv('./AllPts/{}{}.txt'.format(aramprefix,i),comment='#',index_col=None,header=None,sep=',',na_values=' ',skiprows=3).values
        A = A[ ~any(isnan(A),axis=1), :]

    F=A[:,-4:].reshape(-1,2,2) 
    FtF = n.einsum('...ji,...jk',F,F) #Same as that commented out above
    U = mysqrtm( FtF )     #Explicit calculation of sqrtm
    eigU = eigvalsh(U)  
    LEmin, LEmaj = n.sort( n.log(eigU), axis=1 ).T    #Major log and minor log strain
    LEeq = ne.evaluate('( 2/3 * ( LEmin**2 + LEmaj**2 + (-LEmin-LEmaj)**2 ) )**0.5') #epeq
   
    # Shift x,y so max point is at zero in deformed coords
    loc_zmax = n.argmax(A[:,7])
    A[:,[2,3,5,6]]-=A[loc_zmax,[2,3,5,6]]
   
    # Find the point furthest from the max
    if i == last:
        loc_dmax = n.argmax(cdist( A[:,[5,6]], A[loc_zmax,[5,6]][None,:] ).flatten())
        xend,yend = A[loc_dmax,[2,3]]
        
    # Need to interp deformed x,y,z,and epeq based on the undeformed x,y line
    ### The line will being at (0, 0) and pass thru the apex (xend, yend)
    x = linspace(0,xend,interpnumber)
    y = (yend/xend)*x
    xy = vstack((x,y)).T
    
    xint,yint,zint,LEint = (LinearNDInterpolator( A[:,[2,3]], hstack(( A[:,[5,6,7]],LEeq[:,None] )) ).__call__(xy)).T
    dist = ne.evaluate('sqrt(xint**2+yint**2)')
    data = vstack(( dist,zint,LEint )).T
    #data = data[ ~any(isnan(data),axis=1), :]
    # Stack new stages first because we started iterating from last
    profs = hstack(( data, profs ))
    
    # Filter the profiles
    winlen=13
    rgn = ~n.isnan(data[:,0])
    rold = data[rgn,0]
    zold = data[rgn,1]
    eold = data[rgn,2]
    rmin, rmax = n.min(rold),n.max(rold)
    r = n.linspace(rmin,rmax,interpnumber)
    zint = interp1d(rold,zold).__call__(r)
    eint = interp1d(rold,eold).__call__(r)
    e = SG(eint,winlen,1)
    profs_filt = hstack((r[:,None],zint[:,None],e[:,None],profs_filt))
    
    if makecontours:
        if not os.path.exists('./Contours'):
            os.mkdir('./Contours')
        p.figure(facecolor='w',figsize=(12,9))
        p.tricontourf(A[:,2],A[:,3],LEeq,256,cmap='viridis')
        p.axis('equal')
        cbar = p.colorbar()
        cbar.set_label('$\\mathsf{e}_{\\mathsf{e}}$',rotation=0,fontsize=20)
        p.savefig('./Contours/Stg{}.png'.format(i),bbox_inches='tight')
        p.close()
    
####################################
############ Save data #############    
if savedata:
    # Profiles.dat
    headerline = 'Profiles at profStgs\n[0]Distance along profile [1]Z-coord along profile [2]Equiv. Stn'
    n.savetxt('Profiles.dat', X=profs, fmt='%.8f', delimiter=',', header=headerline)
    # Profiles_Filtered.dat
    headerline = 'Profiles at profStgs\n[0]Distance along profile [1]Z-coord along profile [2]Equiv. Stn'
    n.savetxt('Profiles_Filtered.dat', X=profs_filt, fmt='%.8f', delimiter=',', header=headerline)
    # ProfStages.dat
    headerline = 'Stages at which profiles were generated'
    n.savetxt('profStgs.dat',X=profStgs,fmt='%.0f',delimiter=',')
    
os.chdir('../PyFiles')