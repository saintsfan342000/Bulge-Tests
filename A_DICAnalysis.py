import numpy as n
from numpy import (any, abs, nanmean, nanstd, isnan, nanmax, nanmin,
                    array, vstack, hstack,polyfit,exp,linspace)
from numpy.linalg import eigvalsh, eigh
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
Bulge Test Analysis
-Homogeneous best-fit sphere radius
-Radii in rolling and transverse direction
-Major, minor log stn
-Rolling, transv. stn
-Gamma
-StdDev of major and minor
-dt

See http://stackoverflow.com/questions/19868548/set-line-colors-according-to-colormap
and http://stackoverflow.com/questions/8945699/gnuplot-linecolor-variable-in-matplotlib/18516488#18516488
and http://matplotlib.org/examples/pylab_examples/multicolored_line.html
for instructions on plotting line with color based on third value
'''

expt = 1
FS = 40
SS = 10
thickness = 0.0474
E, nu, yldsts = 10136.41*1000, 0.3261, 17*1000
savedata = True
plotdata = True
makecontours = True
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

if (saveAram) and (not os.path.exists('AramisBinary')):
    os.mkdir('AramisBinary')

if not os.path.exists('STLP.dat'):
    ###########################################
    ########### OFFICIAL POLICY ###############
    ### is to not correct the labview file. ###
    ### Leave the pressure -60 in the LV file #
    ### Correct it in STLP.dat #################
    ST = read_csv('./Misc/ArST.dat',sep=',',header=None,comment='#').values
    #[0] Stage [1]Time
    LV = read_csv('./Misc/LV_raw.dat',comment='#',header=None,delim_whitespace=True).values
    #[0]Time [1]LVDT (in) [2]Pressure (psi)
    LV[:,2] += 60   #Correct the pressure
    LVint = interp1d(LV[:,0],LV[:,1:],axis = 0).__call__(ST[:,1])
    STLP = n.hstack( (ST,LVint) )
    header = '[0]DIC Stage [1]Time [2]LVDT (in) [3]Corrected pressure (psi)'
    n.savetxt('STLP.dat',X=STLP,fmt='%.0f, %.4f, %.8f, %.8f',header=header)
else:    
    #'[0]DIC Stage [1]Time [2]LVDT (in) [3]Corrected pressure (psi)'
    STLP = read_csv('STLP.dat',sep=',',header=None,comment='#').values
    
headerline = '[0]BulgeHeight, [1]Sphere Rad [2]Radius-Rolling [3]Radius-Transv. [4]MajorLogStn [5]MinorLogStn \n[6]ex [7]ey [8]Gamma  [9]LEeq [10]Thickness (incomp.) [11]Thickness (iterated) [12]Sdev MajStn [13]TotalPts [14]Filtered Pts'

export = n.empty( (last+1, 13) )

numpts = n.empty( (last+1,2) )
    
##################################
######## BEGIN STAGE LOOP ########
##################################
tnew = thickness
for i in range(0,last+1):
    print(i,end=',  ')
    if BIN:
        A = n.load('./AramisBinary/{}{}.npy'.format(aramprefix,i))
    else:
        A = read_csv('./AllPts/{}{}.txt'.format(aramprefix,i),comment='#',index_col=None,header=None,sep=',',na_values=' ',skiprows=3).values
        A = A[ ~any(isnan(A),axis=1), :]
        if saveAram:
            n.save('./AramisBinary/BT-{}_{}'.format(expt,i),A)
        #[0]Index_x [1]Index_y [2,3,4]Undef_X,Y,Z inches [5,6,7]Def_X,Y,Z inches [8,9,10,11]DefGrad (11 12 21 22) *)
            
    if i == 0:
        export[i] = [0]*13
        export[i,[1,2,3]] = n.nan
        export[i,[10,11]] = [thickness, thickness]
    else:
        # Find max point
        h = n.max(A[:,7])
        loc = n.argmax(A[:,7])
        # Shift x,y so max point is at zero in deformed coords
        A[:,[2,3,5,6]]-=A[loc,[2,3,5,6]]
        # Points within 0.6 and 0.3 inches of max        
        dists = cdist( A[:,[5,6]], A[loc,[5,6]][None,:] ).ravel()
        A6 = A[ dists <= 0.6, :]
        A7 = A[ dists <= 0.7, :]
        A3 = A[ dists <= 0.3, :]

        # Best-fit sphere rad to 0.6
        R_sphere = SphereFit( A6[:,[5,6,7]] )
        # Rolling and transverse radii
        ### Example: X-direction will all have y-coord zero
        ### Then the X-Z plane is the fitting plane
        xspace = linspace(nanmin(A6[:,5]),  nanmax(A6[:,5]), 250)
        rolldir_Zint = griddata( A7[:,[5,6]], A7[:,7],( xspace[None,:], array([[0]]) ) ).ravel()
        R_x = CircleFitByPratt( n.vstack( (xspace,rolldir_Zint) ).T )[-1]
        yspace = linspace(nanmin(A6[:,6]),  nanmax(A6[:,6]), 250)
        transvr_Zint = griddata( A7[:,[5,6]], A7[:,7], (array([[0]]), yspace[:,None] ) ).ravel()
        R_y = CircleFitByPratt( n.vstack( (yspace,rolldir_Zint) ).T )[-1]

        # Calculate strains...do the whole field since we'll need epeq along a diagonal meridian
        F=A3[:,-4:].reshape(-1,2,2) 
        FtF = n.einsum('...ji,...jk',F,F) #Same as that commented out above
        U = mysqrtm( FtF )     #Explicit calculation of sqrtm
        eigU, V = eigh(U)
        LEprin = n.log(eigU)
        # Now rotate the principle log strains back to x, y using the eigenvectors
        # It turns out that LEx, and LEy are very close to just taking n.log(U)
        # And the resulting shears are very close to the off-diagonals of U
        LExy = n.einsum('...ij,...j,...jk',V,   LEprin,V)
        LEx, LEy = LExy[:,0,0], LExy[:,1,1] 
        LEmin, LEmaj = n.sort( n.log(eigU), axis=1 ).T    #Major log and minor log strain
        LEeq = ne.evaluate('( 2/3 * ( LEmin**2 + LEmaj**2 + (-LEmin-LEmaj)**2 ) )**0.5') #epeq
        NEx = U[:,0,0] - 1
        NEy = U[:,1,1] - 1
        NExy = U[:,0,1]
        gam = n.arctan(NExy/(1+NEx)) + n.arctan(NExy/(1+NEy))

        strain = n.vstack(( A3[:,5], A3[:,6],LEmaj,LEmin,LEx,LEy,gam,LEeq)).T
        #[0]Major [1]Minor [2]Rolling [3]Transverse [4]Shear
        rat = LEx/LEy
        ratmean = nanmean( rat )
        ratsdev = nanstd( rat )
        keepers = (rat >= (ratmean-.75*ratsdev)) & (rat <= (ratmean+.75*ratsdev) )
        strain = strain[keepers, :]
        X,Y,LEmaj,LEmin,LEx,LEy,gam,LEeq = strain.T
        
        #Interpolate LEx and LEy such that they are along their respective directions
        xspace = n.linspace(n.min(X),n.max(X),len(X))
        LExOnX = griddata( n.vstack((X,Y)).T, LEx,( xspace[None,:], array([[0]]) ) ).ravel()
        yspace = n.linspace(n.min(Y),n.max(Y),len(Y))
        LEyOnY = griddata( n.vstack((X,Y)).T, LEy, (array([[0]]), yspace[:,None] ) ).ravel()
        
        nomsts = STLP[i,-1]*R_sphere/(2*tnew)
        if nomsts <= yldsts:
            e3 = -2*nu*nomsts/E
            tnew = thickness*exp(e3)
        else:
            #TrueSts(P,e1,e2,R1,R2,E,v,to)
            tnew_it = TrueSts(STLP[i,-1], nanmean(LEmaj), nanmean(LEmin), R_sphere, R_sphere, E, nu, thickness)
            #tnew_inc = TrueSts(STLP[i,-1], nanmean(LEmaj), nanmean(LEmin), R_sphere, R_sphere, E, 0.5, thickness)
            tnew_inc = tnew=thickness*exp(-(nanmean(LEmaj)+nanmean(LEmin) ) )
            

        #'[0]BulgeHeight, [1]Sphere Rad [2]Radius-Rolling [3]Radius-Transv. [4]MajorLogStn [5]MinorLogStn \n[6]ex [7]ey [8]Gamma 
        #[9]Equiv.Stn [10]Thickness (incomp.) [11]Thickness (iterated) [12]Sdev MajStn [13]TotalPts [14]Filtered Pts'
        export[i] = [h, R_sphere, R_x, R_y, nanmean(LEmaj), nanmean(LEmin), nanmean(LExOnX), nanmean(LEyOnY), nanmean(abs(gam)), nanmean(LEeq), tnew_inc, tnew_it, nanstd(LEmaj)]

##############################################
### Generate what we need for meridian strains
##############################################

# Want 4 stages before limit load, one at, and two after
LL = n.argmax(STLP[:,-1])
h = export[:,0]
hstgs = n.hstack(( n.linspace(h[0],h[LL],6)[1:-1], n.linspace(h[LL],h[-1],3) ))
profStgs = n.zeros(len(hstgs))
for k,i in enumerate(hstgs):
    profStgs[k] = n.nonzero(h>=i)[0][0]
profStgs = profStgs.astype(int)

interpnumber = 200
profiles = n.empty((interpnumber,0))
k=0
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
    profiles = hstack(( data,profiles ))
    
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
    # Results.dat
    headerline = '[0]BulgeHeight, [1]Sphere Rad [2]Radius-Rolling [3]Radius-Transv. [4]MajorLogStn [5]MinorLogStn \n[6]ex [7]ey [8]Gamma  [9]LEeq [10]Thickness (incomp.) [11]Thickness (iterated) [12]Sdev MajStn'    
    n.savetxt('Results.dat',X=export,fmt='%.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f',header=headerline)
    # Profiles.dat
    headerline = 'Profiles at profStgs\n[0]Distance along profile [1]Z-coord along profile [2]Equiv. Stn'
    n.savetxt('Profiles.dat', X=profiles, fmt='%.8f', delimiter=',', header=headerline)
    # ProfStages.dat
    headerline = 'Stages at which profiles were generated'
    n.savetxt('profStgs.dat',X=profStgs,fmt='%.0f',delimiter=',')

if plotdata:
    a = os.system('python ../PyFiles/B_Plots.py {} {:.0f}'.format(str(plotdata),expt))
    if a == 1:
        os.system('python3 ../PyFiles/B_Plots.py {} {:.0f}'.format(str(plotdata),expt))
