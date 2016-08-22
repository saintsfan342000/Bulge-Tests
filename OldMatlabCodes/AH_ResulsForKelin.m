clear
clc
close all

curdir=pwd;
addpath(sprintf('%s\\MATLAB\\extras',curdir(1:2)));


paths = {'BT-2_Recalc_Results/','BT-3_StgZeroDeleted_Results/','BT-4_Results/','BT-5_Results/R1p8/'};
paths = {paths{4}};
lstring = {'BT-2','BT-3','BT-4'};

% (1)Stage (2)Pressure (3)Bulge Height (4)Rx (5)Ry (6)ex (7)ey (8)eps_xy

path = ['../../Bulge Tests/' 'BT-6_Results/'];

testres = load([path 'TestResults.dat']); 
    R = testres(:,1);
    Z = testres(:,11);

STLP = load( [path 'STLP.dat'] );
    P = STLP(:,4);

xyrd = load( [path 'Roll-Trans-Radius.dat'] );
    Rx = xyrd(:,1);
    Ry = xyrd(:,2);
    ex = xyrd(:,3);
    ey = xyrd(:,6);
    exy = mean(xyrd(:,[7 8]),2);
    
exp = [STLP(:,1)  P  Z  Rx  Ry  ex  ey  exy R];

fid = fopen([path 'BT-6_ForKelin.dat'],'w');
fprintf(fid,'%.0f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f\n',exp');
fclose(fid);
