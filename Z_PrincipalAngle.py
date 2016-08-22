import numpy as n
from numpy import (any, abs, nanmean, nanstd, isnan, nanmax, nanmin,
                    array, vstack, hstack,polyfit,exp,linspace)
from numpy.linalg import eigh
import matplotlib.pyplot as p
from pandas import read_csv
from scipy.interpolate import interp1d, griddata, LinearNDInterpolator
from scipy.spatial.distance import cdist
import numexpr as ne
import os, glob
from mysqrtm import mysqrtm 
from CircleFitByPratt import CircleFitByPratt
from SphereFit import SphereFit
from TrueSts import TrueSts

'''
Calculate the principal strain angle
And make histograms
'''

expt = 1
FS = 40
SS = 10
thickness = 0.0474
E, nu, yldsts = 10136.41*1000, 0.3261, 17*1000
savedata = False
plotdata = False
BIN = True

savepath = '../BT-{}_FS{}SS{}'.format(expt,FS,SS)
if BIN:
    arampath = '{}/AramisBinary'.format(savepath)
    aramprefix = 'BT-{}_'.format(expt)
else:
    arampath = '{}/AllPts'.format(savepath)
    aramprefix = 'BT-{}-Stage-0-'.format(expt)
    
os.chdir(savepath)

last = len(glob.glob('{}/{}*'.format(arampath,aramprefix))) - 1

##################################
######## BEGIN STAGE LOOP ########
##################################
tnew = thickness
stages = n.sort( n.unique( hstack( (n.arange(20,last+20,20),[last]) ) ) )
stages = stages[ stages <= last+1]
#for i in range(len(last+1)):
for i in stages:
    print(i,end=',  ')
    if i == 0:
        pass
    else:
        if BIN:
            A = n.load('./AramisBinary/{}{}.npy'.format(aramprefix,i))
        else:
            A = read_csv('./AllPts/{}{}.txt'.format(aramprefix,i),comment='#',index_col=None,header=None,sep=',',na_values=' ',skiprows=3).values
            A = A[ ~any(isnan(A),axis=1), :]
            if saveAram:
                n.save('./AramisBinary/BT-{}_{}'.format(expt,i),A)

        #[0]Index_x [1]Index_y [2,3,4]Undef_X,Y,Z inches [5,6,7]Def_X,Y,Z inches [8,9,10,11]DefGrad (11 12 21 22) *)

        # Find max point
        h = n.max(A[:,7])
        loc = n.argmax(A[:,7])
        # Shift x,y so max point is at zero in deformed coords
        A[:,[2,3,5,6]]-=A[loc,[2,3,5,6]]
        # Points within 0.3 inches of max        
        dists = cdist( A[:,[5,6]], A[loc,[5,6]][None,:] ).ravel()
        A3 = A[ dists <= 0.3, :]

        # Calculate strains...do the whole field since we'll need epeq along a diagonal meridian
        F=A3[:,-4:].reshape(-1,2,2) 
        FtF = n.einsum('...ji,...jk',F,F) #Same as that commented out above
        U = mysqrtm( FtF )     #Explicit calculation of sqrtm
        LE = n.log(U)           # Log strains from stretch
        LEx, LEy = LE[:,0,0], LE[:,1,1] #ex, ey
        rat = LEx/LEy
        ratmean = nanmean( rat )
        ratsdev = nanstd( rat )
        keepers = (rat >= (ratmean-.75*ratsdev)) & (rat <= (ratmean+.75*ratsdev) )
        U = U[keepers,:,:]
        valU, vecU = eigh(U)
        sorted_vals =  n.argsort( valU, axis=1 )
        if any( sorted_vals[:,0] == 1 ):
            for i in range(len(valU[:,0])):
                if (sorted_vals[i,0]==1):
                    vecU[i,:,:] = n.fliplr( vecU[i,:,:] )
        # Now the stack of eigenvectors are sorted such that second column in each pizza box
        # is the eigenvector corresponding the the major direction
        # So I need to take the arctan of [1,1]/[0,1] in each pizza box
        princAng = n.arctan2(vecU[:,1,1],vecU[:,0,1])*180/n.pi
        princAng[princAng<0]+=180
        print(len(princAng))
        p.figure(facecolor='w')
        p.hist(princAng,bins=180/15,normed=True,alpha=0.6,range=(0,180))
        p.title('Principal Angle Distribution: Stage {}'.format(i))
        p.savefig('Histograms/Stage {}'.format(i))
        p.close()
            