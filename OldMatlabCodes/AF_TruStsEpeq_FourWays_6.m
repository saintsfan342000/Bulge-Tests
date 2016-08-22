function zz_TruStsEpeq_FourWays_6();

% Must turn mmrad back on after principal direction is decided!

%Calls the anisotropic tru-sts-epeq function for four different cases
%Uses the radii are calculated using all points above thresh=1.545
% Strains, though, are taken from even higher, thresh=1.57

savefiles=1;


%%%% Path of the experiment folder
    frompath='F:\Martin_Experiments\Sandia_Anis\Bulge Tests\BT-6_Results';
    savepath=frompath;
%%%% Relative path and prefix of the cleaned aramis files
    prefix='AramisExport_MissingRemoved\BT-6-Stage-0-';
%%%% Last stage
    last=412;
    to = .04;
    E = 10191;
    nu = 0.316;
%%%% Yield Sts
    yld=44083;
%%%% Facet size
    FS=75;
    SS=8;
    BT=2;
    
    distfrommaxz = 1.2/2;
    distfrommaxzforstn = 0.3;

    
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
    STLP(:,4)=STLP(:,4);

%Open up the rolling-dir'n and transverse dir'n radii
    xyrad=load(sprintf('%s\\Roll-Trans-Radius.dat',savepath));
   % mmrad=load(sprintf('%s\\Maj-Min-Radius.dat',savepath));
    
z=1;
pss_1=[];pss_2=[];pss_3=[];pss_4=[];
for i=1:length(STLP(:,1));
    STLP(i,1)
    %load the export file
        clear A mmRatio xyRatio phi tru epeq
        A=load(sprintf('%s\\%s%d.txt',frompath,prefix,STLP(i,1)));
        A(:,5) = A(:,5) + 0.91885186/25.4;
   %FILTER:  Get rid of major/minor strains or x/y strains that are less than zero
        A(A(:,6)<0 | A(:,7)<0 | A(:,11)<0 | A(:,12)<0 ,:)=[];    

     % FILTER: Strain ratio tolerance
     % Major minor ratio
         mmRatio=A(:,7)./A(:,6);
     % x-y Ratio
         xyRatio=A(:,12)./A(:,11);
     % Eliminate erroneous ratios from A and ratio
         A(xyRatio>1.5 | xyRatio<1/1.5 | mmRatio<1/1.5,:)=[];
         clear mmRatio xyRatio
    % emaj=A(:,6);emin=A(:,7);ex=A(:,11);ey=A(:,12);exy=A(:,13);

    [~,locz]=max(A(:,5));
    dz=sqrt((A(:,3)-A(locz,3)).^2+(A(:,4)-A(locz,4)).^2);
    A(dz>distfrommaxz,:)=[];
    dz(dz>distfrommaxz,:)=[];
    
    %Compute best fit sphere based on the whole exported point cloud
        XYZ=A(:,[3 4 5]);
        [ctr(i,:),rad(i,1)]=sphereFit(XYZ);
        
        A(dz>distfrommaxzforstn,:)=[];
        
    P2t=STLP(i,4)/(2*.04);
        
    % Go to Tru sts-plastic stn algorithm
    if P2t*rad(i) > yld
    %[trunew1, trunew2,e1p,e2p,e3p,epeq,k,tnew]=TrueStress_NotIsoTropic(P,e1,e2,R1,R2,E,v,to);
        
        % Iso Stn, Iso Radius
        clear trunew1 trunew2 e1p e2p e3p epeq tnew
        [trunew1, trunew2,e1p,e2p,e3p,epeq,k,tnew]=TrueStress_NotIsoTropic(STLP(i,4) , mean([A(i,6),A(i,7)]) , mean([A(i,6),A(i,7)]) , rad(i) , rad(i) , E , nu , to);
        pss_1=[pss_1;i 0 0 0 0 STLP(i,4)*rad(i)/(2*to) STLP(i,4)*rad(i)/(2*to) tnew];
        
        % Maj-Min stn, Iso Radius
        clear trunew1 trunew2 e1p e2p e3p epeq tnew
        [trunew1, trunew2,e1p,e2p,e3p,epeq,k,tnew]=TrueStress_NotIsoTropic(STLP(i,4) , A(i,6) , A(i,7) , rad(i) , rad(i) , E , nu , to);
        pss_2=[pss_2;i e1p e2p e3p epeq trunew1/1000 trunew2/1000 tnew];
    else
        pss_2=[pss_2;i nan nan nan nan STLP(i,4)*rad(i)/(2*to) STLP(i,4)*rad(i)/(2*to) -2*nu*(STLP(i,4)*rad(i)/(2*to))/E];
    end;
    if P2t*xyrad(i,1) > yld && P2t*xyrad(i,2) > yld
        % X and Y stn, X and Y rad
        clear trunew1 trunew2 e1p e2p e3p epeq tnew
        [trunew1, trunew2,e1p,e2p,e3p,epeq,k,tnew]=TrueStress_NotIsoTropic(STLP(i,4),mean(xyrad(i,[3])),mean(xyrad(i,[6])),xyrad(i,1),xyrad(i,2),E,.316,to);
        pss_3=[pss_3;i e1p e2p e3p epeq trunew1/1000 trunew2/1000 tnew];
    else
        pss_3=[pss_3;i nan nan nan nan STLP(i,4)*xyrad(i,1)/(2*to)/1000 STLP(i,4)*xyrad(i,2)/(2*to)/1000 -(nu/E)*STLP(i,4)*(xyrad(i,1)+xyrad(i,2))/(2*to) ];
    end;
%     if  P2t*mmrad(i,1) > yld && P2t*mmrad(i,2) > yld      
%         % Maj and Min stn, Maj and Min rad
%         clear trunew1 trunew2 e1p e2p e3p epeq tnew
%         [trunew1, trunew2,e1p,e2p,e3p,epeq,k,tnew]=TrueStress_NotIsoTropic(STLP(i,4) , mean(mmrad(i,[3])) , mean(mmrad(i,[6])) , mmrad(i,1) , mmrad(i,2) , E , nu , to);
%         pss_4=[pss_4;i e1p e2p e3p epeq trunew1/1000 trunew2/1000 tnew];
%        
%     end
end;    

%  figure(1)
% % plot(pss_4(:,2),pss_4(:,end-1),'r')
% % hold on
% % plot(pss_4(:,3),pss_4(:,end-0),'b')
% plot(pss_3(:,2),pss_3(:,end-1),'k')
% hold on
% plot(pss_3(:,3),pss_3(:,end-0),'g')
% %plot([pss_4(60,3) pss_4(60,3)],[50 80],'k--')
% title(sprintf('\\tau - e:  Rad=%.2f,Stn=%.2f',distfrommaxz,distfrommaxzforstn),'fontsize',14)
% set(gca,'Ylim',[45 90])
% legend('Rolling Direction','Transverse Direction')
% %l=legend('MajorStn Direction','MinorStn Direction','Rolling Direction','Transverse Direction')
% %set(l,'location','southeast')
% print(gcf,'-dpdf',sprintf('%s\\Tru-Sts-Plastic-Stn',savepath))
% 
% figure(2)
% %plot(pss_4(:,2),'r')
% hold on
% %plot(pss_4(:,3),'b')
% plot(pss_3(:,2),'k')
% plot(pss_3(:,3),'g')
% legend('Rolling Direction','Transverse Direction')
% %l=legend('MajorStn Direction','MinorStn Direction','Rolling Direction','Transverse Direction')
% %set(l,'location','Northwest')
% title(sprintf('Strain at Each Stage:  Rad=%.2f,Stn=%.2f',distfrommaxz,distfrommaxzforstn))
% print(gcf,'-dpdf',sprintf('%s\\Plastic Strain Evolution',savepath))

fid1 = fopen(sprintf('%s\\Case1.dat',savepath),'w');
fprintf(fid1,'%.0f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f\n',pss_1')
fclose(fid1);

fid2 = fopen(sprintf('%s\\Case2.dat',savepath),'w');
fprintf(fid2,'%.0f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f\n',pss_2')
fclose(fid2);

fid3 = fopen(sprintf('%s\\Case3.dat',savepath),'w');
fprintf(fid3,'%.0f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f\n',pss_3')
fclose(fid3);

figure()
plot(pss_3(:,2),pss_3(:,6),'b')
hold on
plot(pss_3(:,3),pss_3(:,7),'r')
legend('Roll','Trans','location','southeast')
set(gcf,'color','w')
axis([00 .16 20 90])
export_fig([savepath '/Sts-Stn.png'],'-nocrop','-r250')

% fid4 = fopen(sprintf('%s\\Case4.dat',savepath),'w');
% fprintf(fid4,'%.0f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f\n',pss_4')
% fclose(fid4);
