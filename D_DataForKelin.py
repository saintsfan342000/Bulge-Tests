import numpy as n
L = n.genfromtxt
fit = n.polyfit
abs=n.abs
from numpy import nanmax, nanmin
import os
import pandas as pd
import shutil

exporthead = 'DIC Stage, Time, Pressure(psi), Sphere Radius, Radius X, Radius Y, LogStn X, LogStn Y, Gamma, Log Mises Stn (total), Bulge Height'
exporthead = exporthead.split(',')
#So I just need to open LV_DIC and InSGZone

expts = ['1_FS40SS10','2_FS40SS10']

fid = pd.ExcelWriter('../../AlcoaData_Kelin/BulgeTest.xlsx')
misc = n.empty(len(expts)) #Store the misc properties

for G in range(len(expts)):
    print(expts[G])
    # STLP: Stage, Time, Pressure(psi),
    STLP = L('../BT-{}/STLP.dat'.format(expts[G]),delimiter=',')[:,[0,1,3]]
        
    #Results.dat
    # [0]BulgeHeight, [1]Sphere Rad [2]Radius-Rolling [3]Radius-Transv. [4]MajorLogStn [5]MinorLogStn 
    # [6]ex [7]ey [8]Gamma  [9]LEeq [10]Thickness (incomp.) [11]Thickness (iterated) [12]Sdev MajStn [13]TotalPts [14]Filtered Pts
    A = L('../BT-{}/Results.dat'.format(expts[G]),delimiter=',')[:,[1,2,3,6,7,8,9,0]]
    
    d = n.hstack(( STLP,A ))
    
    #print(expdata.shape)
    #print(len(exporthead))
    pd.DataFrame(d).to_excel( fid, sheet_name='BT-{}'.format(expts[G].split('_')[0]), index=False, header=exporthead)

fid.save()