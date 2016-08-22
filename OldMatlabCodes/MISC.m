function MISC();

%Calls the anisotropic tru-sts-epeq function for four different cases
%Uses the radii are calculated using all points above thresh=1.545
% Strains, though, are taken from even higher, thresh=1.57

savefiles=0;


%%%% Path of the experiment folder
    frompath='F:\Martin_Experiments\Sandia_Anis\Bulge Tests\BT-2_Results\Final Computation - Radius0p6-Stn0p3';
    
%%%% Yield Sts
    yld=44083;

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
    xyrad=load(sprintf('%s\\Roll-Trans-Radius.dat',frompath));
    mmrad=load(sprintf('%s\\Maj-Min-Radius.dat',frompath));
    
% Exported are all points above Z=1.545 in the final stage.  To those
% points I will fit the sphere.  But I will furhter reduce to a higher
% Z-coordinate threshhold.  Now I'll open up the last stage, find these IJ
% indices, and keep only these points to compute strain

z=1;
pss_1=[];pss_2=[];pss_3=[];pss_4=[]; dt=[];
for i=1:length(STLP(:,1));
    %load the export file
        clear A mmRatio xyRatio phi tru epeq
    P2t=STLP(i,4)/(2*.0399);
        
    % Go to Tru sts-plastic stn algorithm
    if P2t*mmrad(i,1) > yld && P2t*mmrad(i,2) > yld
        % X and Y stn, X and Y rad
        clear trunew1 trunew2 e1p e2p e3p epeq tnew
        [~, ~,~,~,~,~,~,tnew]=TrueStress_NotIsoTropic(STLP(i,4),mean(mmrad(i,[3])),mean(mmrad(i,[6])),mmrad(i,1),mmrad(i,2),10191,.316,.0399);
        dt = [dt;(.0399 - tnew)/.0399];
    else
        clear tnew
        ez = -.316*(mean(mmrad(i,[3])) + mean(mmrad(i,[6])));
        dt=[dt; abs(ez)];
    end;

end;    



%{
figure(1)
plot(pss_4(:,2),pss_4(:,end-1),'r')
hold on
plot(pss_4(:,3),pss_4(:,end-0),'b')
plot(pss_3(:,2),pss_3(:,end-1),'k')
plot(pss_3(:,3),pss_3(:,end-0),'g')
%plot([pss_4(60,3) pss_4(60,3)],[50 80],'k--')
title(sprintf('\\tau - e:  Rad=%.2f,Stn=%.2f',distfrommaxz,distfrommaxzforstn),'fontsize',14)
set(gca,'Ylim',[45 90])
l=legend('MajorStn Direction','MinorStn Direction','Rolling Direction','Transverse Direction')
set(l,'location','southeast')
% print(gcf,'-dpdf',sprintf('%s\\Tru-Sts-Plastic-Stn',savepath))

figure(2)
plot(pss_4(:,2),'r')
hold on
plot(pss_4(:,3),'b')
plot(pss_3(:,2),'k')
plot(pss_3(:,3),'g')
l=legend('MajorStn Direction','MinorStn Direction','Rolling Direction','Transverse Direction')
set(l,'location','Northwest')
title(sprintf('Strain at Each Stage:  Rad=%.2f,Stn=%.2f',distfrommaxz,distfrommaxzforstn))
% print(gcf,'-dpdf',sprintf('%s\\Plastic Strain Evolution',savepath))
% 
% save(sprintf('%s\\Case1.dat',savepath),'pss_1','-ascii')
% save(sprintf('%s\\Case2.dat',savepath),'pss_2','-ascii')
% save(sprintf('%s\\Case3.dat',savepath),'pss_3','-ascii')
% save(sprintf('%s\\Case4.dat',savepath),'pss_4','-ascii')
%}