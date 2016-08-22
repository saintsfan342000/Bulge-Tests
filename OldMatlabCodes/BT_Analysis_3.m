function BT_Analysis_3()

savefiles=0;

close all

%%%% Path of the experiment folder
    frompath='E:\Martin_Experiments\Sandia_Anis\Bulge Tests\BT-2_Recalc_Results';
    savepath=frompath;
%%%% Relative path and prefix of the cleaned aramis files
    prefix='AramisExport_MissingRemoved\BT2-Recalc-Stage-0-';
%%%% Last stage
    to = .04;
    last=355;
%%%% Yield Sts
    yld=44083;
%%%% Facet size
    FS=75;
    SS=8;
    BT=5;
    
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
    
%Find stages nearest to the Pressures at which we want to plot histograms by finding the nearest pressure
    %Pressures at which we want to plot things
pvals=[250 500 750 1000:75:1450];
for i=1:length(pvals);
    [~,ploc(i)]=min(abs(STLP(:,4)-pvals(i)));
    %ploc is an array that contains the row indexes in STLP where the
    %pressure is closest.  So STLP(ploc(i),1) gives the corresponding
    %stage number
end;

z=1;
a=[];b=[];c=[];
for i=1:length(STLP(:,1));
    i
    %load the export file
        clear A mmRatio xyRatio phi tru epeq
        A=load(sprintf('%s\\%s%d.dat',frompath,prefix,STLP(i,1)));
        
       A(:,5) = A(:,5) + 0.91885186/25.4;
        
        [~,locz] = max(A(:,5));
        
        Acopy = A;
        
        dz=sqrt((A(:,3)-A(locz,3)).^2+(A(:,4)-A(locz,4)).^2);
        A(dz>distfrommaxzforstn,:)=[];    
        Acopy(dz>distfrommaxz,:)=[];
        
        a=length(A(:,1));
        
   %FILTER:  Get rid of major/minor strains or x/y strains that are less than zero
        A(A(:,6)<0 | A(:,7)<0 | A(:,11)<0 | A(:,12)<0 ,:)=[];    
   % Go through A, row-by-row to convert major strain direction to unit vector
   for j=1:size(A,1);
        A(j,[8 9 10])=A(j,[8 9 10])/norm(A(j,[8 9 10]));  %Convert the major strain direction to unit vector
        if A(j,8)<0  && acosd(dot(A(j,[8 9 10]),[1 0 0])) > 135;
            phi(j,1)=180 - acosd(dot(A(j,[8 9 10]),[1 0 0]));
        elseif A(j,9)<0  && acosd(dot(A(j,[8 9 10]),[1 0 0])) > 45 && acosd(dot(A(j,[8 9 10]),[1 0 0])) <135;
            phi(j,1)=180 - acosd(dot(A(j,[8 9 10]),[1 0 0]));
        else
            phi(j,1)=acosd(dot(A(j,[8 9 10]),[1 0 0]));
        end;
    end;

    % FILTER: Strain ratio tolerance
    % Major minor ratio
        mmRatio=A(:,7)./A(:,6);
    % x-y Ratio
        xyRatio=A(:,12)./A(:,11);
    % Eliminate erroneous ratios from A and ratio
        A(xyRatio>1.5 | xyRatio<1/1.5 | mmRatio<1/1.5,:)=[];
        
        b = length(A(:,1));
        c = [c;a b];
        
        clear mmRatio xyRatio
    % emaj=A(:,6);emin=A(:,7);ex=A(:,11);ey=A(:,12);exy=A(:,13);

    %Compute best fit sphere based on the whole exported point cloud
        XYZ=Acopy(:,[3 4 5]);
        [ctr(i,:),rad(i,1)]=sphereFit(XYZ);
        
        mmRatio=A(:,7)./A(:,6);
        xyRatio=A(:,12)./A(:,11);
    
        thickness = to*exp( -mean( A(:,6) )- mean( A(:,7) ) ) ;
    
    %(1)Rad (2)MeanMajorStn (3)MeanMinorStn (4)Mean min/maj Ratio (5)Sdev Maj/Min Ratio (6)Mean principal direction
    %(7)Mean ex (8)Mean ey (9)Mean ey/ex (10)std(ey/ex) 
    D(i,:)=[rad(i) mean(A(:,6)) mean(A(:,7)) mean(mmRatio) std(mmRatio) mean(phi) mean(A(:,11)) mean(A(:,12)) mean(xyRatio) std(xyRatio) max(A(:,5)) thickness to-thickness];
    STDS(i,:) = [nanstd(A(:,11)) nanstd(A(:,12))];
            
    if any(i == ploc) && savefiles==1
        figure
        subplot(1,3,1)
        hist(phi,20)
        title(sprintf('Stage %d. Pressure ~%.0f psi. \\sigma %.0f ksi',STLP(i,1),pvals(z),round(D(i,1))))
        xlabel('\theta w.r.t. Rolling Axis')
        ylabel('N');
        set(gca,'Xlim',[0 135])
        subplot(1,3,2)
        hist(mmRatio,20)
        xlabel('e_m_i_n / e_m_a_j')
        ylabel('N');
        subplot(1,3,3)
        hist(xyRatio,20)
        title(sprintf('BT-%d - FS%d - SS%d',BT,FS,SS))
        xlabel('e_t_r / e_r_o_l_l')
        ylabel('N');
        set(gcf,'color','w')
        set(gcf, 'Units','Inches','Position', [3 3 9 6])
        export_fig(sprintf('%s\\Histograms\\Hist_Stg%d.png',savepath,STLP(i,1)));close;
        close
        z=z+1;
    end;
    
end;
 

    %(1)NomSts (2)MeanMajorStn (3)MeanMinorStn (4)Mean emaj/emin (5)std(emaj/emin (6)Mean principal direction
    %(7)Mean ex (8)Mean ey (9)Mean ey/ex (10)std(ey/ex) (11) Dome Heigh[in] (12)Thickness (13)Change-Thickness



    %(1)NomSts (2)MeanMajorStn (3)MeanMinorStn (4)Mean emaj/emin (5)std(emaj/emin (6)Mean principal direction
    %(7)Mean ex (8)Mean ey (9)Mean ey/ex (10)std(ey/ex) 

    %Strain Ratio vs Pressure
    figure
    mm=plot(STLP(:,4),D(:,4),'color','b','linewidth',2);
    hold on
    xy=plot(STLP(:,4),D(:,9),'linewidth',2,'color','k');
    [maxs,locs]=max(D(:,1));
    plot(STLP(:,4),D(:,4)+D(:,5),'Color',[.5 .5 1])
    plot(STLP(:,4),D(:,4)-D(:,5),'Color',[.5 .5 1])
    plot(STLP(:,4),D(:,9)+D(:,10),'Color',[.5 .5 .5])
    plot(STLP(:,4),D(:,9)-D(:,10),'Color',[.5 .5 .5])
    plot(STLP(ploc,4),D(ploc,4),'ko','MarkerFaceColor','k','MarkerEdgeColor','k');
    plot(STLP(locs,4),D(locs,4),'rs','MarkerFaceColor','r','MarkerEdgeColor','r');
    plot(STLP(ploc,4),D(ploc,9),'ko','MarkerFaceColor','k','MarkerEdgeColor','k');
    plot(STLP(locs,4),D(locs,9),'rs','MarkerFaceColor','r','MarkerEdgeColor','r');
    title(sprintf('Strain Ratio - BT-%d - FS%d - SS%d',BT,FS,SS),'Fontsize',14)
    xlabel ('Pressure (psi)','Fontsize',14)
    l=legend([mm xy],{'e_2/e_1','e_t_r/e_r_o_l_l'});
    set(l,'Location','Southeast')
    if savefiles==1;
    print(gcf,'-dpdf',sprintf('%s\\Strain Ratio',savepath));close;
    end;
    

if savefiles==1;
    % stgdata
        fid=fopen(sprintf('%s\\TestResults.dat',savepath),'w');
        fprintf(fid,'%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f\n',D');
        fclose(fid);clear fid;
    %Spherefit
        fid=fopen(sprintf('%s\\SphereFit.dat',savepath),'w');
        fprintf(fid,'%.8f %.8f %.8f %.8f\n',[ctr rad]');
        fclose(fid);clear fid;

    %Column Labels
        fid=fopen(sprintf('%s\\Column Labels.dat',savepath),'w');
        fprintf(fid,'Test Results:\n(1)BF-Sphere Radius (2)MeanMajorStn (3)MeanMinorStn (4)Mean emaj/emin (5)std(emaj/emin (6)Mean principal direction\n(7)Mean ex (8)Mean ey (9)Mean ey/ex (10)std(ey/ex) (11)Bulge Ht (12)thickness (13)delta-thickness');
        fclose(fid);clear fid;
end;

fname = [savepath '/StdDev.dat'];
fid = fopen(fname,'w');
fprintf(fid,'%.8f,%.8f\n',STDS');
fclose(fid)

        
    % % Calculation of facet size and step based on first stage
% A=load(sprintf('%s\\%s_%d.dat',frompath,prefix,STLP(1,1)));
% x=A(:,3); y=A(:,4); z=A(:,5);
% 
%     min_disp = 1000;    %Initialize to erroneously large value
%     for i = 1:length(A);
%         for j=1:length(A);  
%             dist=((x(i)-x(j))^2 + (y(i)-y(j) )^2 + ( z(i)-z(j) )^2 )^0.5;
%             if dist==0;
%                 dist=1000;
%             end;
%             min_disp = min ( min_disp , dist ); %Recursive process
%         end
%     end
%     SS_in = min_disp;
%     FS_in = FS * (min_disp/SS);
%     clear x y z min_disp
%     % Facet size
%         fid=fopen(sprintf('%s\\FacetSize.dat',savepath),'w');
%         fprintf(fid,sprintf('Step size = %.4f\nFacet size = %.4f',SS_in,FS_in));
%         fclose(fid);clear fid;