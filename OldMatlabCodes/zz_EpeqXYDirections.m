% Here's what I used to generate Epeq-XYDirections-DifferentRadii.dat
savefiles=0;


close all

%%%% Path of the experiment folder
    frompath='E:\Martin_Experiments\Sandia_Anis\Bulge Tests\BT-2_FS75SS25\Thresh1p57';
    savepath='E:\Martin_Experiments\Sandia_Anis\Bulge Tests\BT-2_FS75SS25\Thresh1p57';
%%%% Relative path and prefix of the cleaned aramis files
    prefix='AramisExport_MissingRemoved\BT-2_FS75S25';
%%%% Last stage
    last=218;
%%%% Yield Sts
    yld=44083;
%%%% Facet size
    FS=75;
    SS=25;
    BT=2;
    
% Add extras because I'll need spherefit and AutoAlignFigures
curdir=pwd;
addpath(sprintf('%s\\MATLAB\\extras',curdir(1:2)));

% Aramis file columns
% (1)Index_X (2)Index_Y (3-5)DeformedCoord_X,Y,Z (6-7)Major,Minor Stn
% (8-10)MajorStnDirection_X,Y,Z (11-12) Log Stn_X,Y (13)EpsXY

% STLP columns
% (1)Stage  (2)Time (3)LVDT (4)Pressure

% Load Stage-Time-LVDT-Force
    STLP=load(sprintf('%s\\STLP.dat',frompath));
    %Correct pressure!
    STLP(:,4)=STLP(:,4)+60;

%Open up the rolling-dir'n and transverse dir'n radii
%(1)X-Axis-Aligned Section BF-circle radius (2)Y-aligned radius ....
    xyrad=load(sprintf('%s\\Roll-Trans-Radius.dat',frompath));

%Open up the rolling-dir'n and transverse dir'n mean ex and ey
%(1)Roll Dir. eX (2)Roll Dir eY (3)Trans Dir eX (4) Trans Dir eY
    XY=load(sprintf('%s\\XY_epsX_epsY.dat',frompath));    

z=1;
pssx=[];pssy=[];psst=[];
for i=1:length(STLP(:,1));
        % Go to Tru sts-plastic stn algorithm
    if STLP(i,4)*mean(xyrad(i,[1 2]))/(2*.0399) > yld
        %(P,e1,e2,R1,R2,E,v,to);
        % Considering rolling (x) direction
        %[tru1, tru2,e1p,e2p,e3p,epeq,k]=TrueStress_NotIsoTropic(STLP(i,4),XY(i,1),XY(i,2),xyrad(i,1),xyrad(i,2),10191,.316,.0399);
        %pssx=[pssx;i e1p e2p e3p epeq tru1/1000 tru2/1000];
        % COnsidering transverse (y) direction
        %clear tru1 tru2 e1p e2p e3p epeq k
        %[tru1, tru2,e1p,e2p,e3p,epeq,k]=TrueStress_NotIsoTropic(STLP(i,4),XY(i,3),XY(i,4),xyrad(i,1),xyrad(i,2),10191,.316,.0399);
        %pssy=[pssy;i e1p e2p e3p epeq tru1/1000 tru2/1000];
        % COnsidering the point at the middle
        clear tru1 tru2 e1p e2p e3p epeq k tnew
        [tru1, tru2,e1p,e2p,e3p,epeq,k,tnew]=TrueStress_NotIsoTropic(STLP(i,4),XY(i,1),XY(i,4),xyrad(i,1),xyrad(i,2),10191,.316,.0399);
        psst=[psst;i e1p e2p e3p epeq tru1/1000 tru2/1000 tnew];
    end
end;
 
plot(psst(:,3),psst(:,6),'Color',[.5 .5 .5])
hold on
plot(psst(:,3),psst(:,7),'.','Color',[.5 .5 .5])
plot(pssx(:,3),pssx(:,6),'r')
plot(pssx(:,3),pssx(:,7),'r.')
plot(pssy(:,3),pssy(:,6),'b')
plot(pssy(:,3),pssy(:,7),'b.')
