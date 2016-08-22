function zz_DomeProfiles_3D_Circular();

%close all


%%%% Where are the files located
    path='E:\Martin_Experiments\Sandia_Anis\AA_Paper\BT2_AllPtsForContours\';
    prefix='BT2-Stage-0-';
    %prefix={'X-Axis-Section\X-Sxn_section0_0-';'Y-Axis-Section\Y-Sxn_section1_0-'};
    last=218;
    R = 6/2 ;

    
% Add matlab extras for plane fit, circle fit
curdir=pwd;
addpath(sprintf('%s\\Matlab\\extras',curdir(1:2)));

% Columns of each section file are:
%(1) CoordX (2) CoordY (3) CoordZ (4) Mises Stn (log) (5)Index

    clear X I

    c={[238 201 0]/255,[0 201 87]/355,[0 0 1],[139 58 58]/255,[0 1 0],[238 106 167]/255,[0 1 1],[255 127 36]/255,[0 0 0],[154 50 205]/255};

% First find points in the last stage, and only plot those in stages priod
    fid = fopen(sprintf('%s\\%s%d.txt',path,prefix,last),'r');
    d=cell2mat(textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f','delimiter',',','headerlines',2));
    d(any(isnan(d),2),:)=[];
    keepers= d(:,[1 2]);
    
k = 1;
stg = [];
xdata = [];
zdata = [];
cdata = [];
maxmin = [];
figure

stages = [0:40:160 180 200 218]
stages = [120 180 218]; 
%modstn_stages = [0 40 80 120];

Xgrd = {[] [] []};
Ygrd = {[] [] []};
Zeval={[] [] []};
stneval = {[] [] []};

for i = stages
        
    clear X d dz fid
    % Load the section file
    fid = fopen(sprintf('%s\\%s%d.txt',path,prefix,i),'r');
    d=cell2mat(textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f','delimiter',',','headerlines',2));
    % Index_x  ;;  Index_y  ;;  (3,4,5) Coord-undef (x,y,z)  ;;  (6,7,8) Disp (x,y,z)  ;;  (9,10,11,12) F11 ;; F12 ;; F21 ;; F22
    d(any(isnan(d),2),:)=[];

    
    %d = d(ismember(d(:,[1 2]),keepers,'rows'),:);   % Limit ourselves to the points also in the last stage
    %any(sqrt( (d(:,3)).^2 + (d(:,4)).^2) > 1)
   % d ( sqrt( (d(:,3)).^2 + (d(:,4)).^2) > 1,:) = []; 
    
    
    ArI = d(:,1);
    Arj = d(:,2);
    X = d(:,3)+d(:,6);
    Y = d(:,4)+d(:,7);
    Z = d(:,5)+d(:,8);
    stn = zeros(length(X),1);
    
    %Z ( sqrt( X.^2 + Y.^2) > 1,:) = nan; 
    
    for j = 1:length(X)
        F = [d(j,9) d(j,10);d(j,11) d(j,12)];
        U = sqrtm(F'*F);
        evs = eig(U);
        L1 = evs(1);
        L2 = evs(2);
        e1 = log(L1);
        e2 = log(L2);
        stn(j,1) = sqrt( (2/3) * (e1^2 + e2^2 + (-e1-e2)^2) ); 
    end;
    
    %Creation of a grid in order to plot the data
    %Determination of the min distance between any two facets in the stage
    min_disp = 1000;    %Initialize to erroneously large value
    for m = [1:round(length(X)/50)] %i will go from 2 to length(X)-1 in 30 equally spaced increments
        for q= [1:round(length(X)/50)]    %j will cover all indices except i itself
            dist= ( (X(m)-X(q))^2 + (Y(m)-Y(q) )^2 + ( Z(m)-Z(q) )^2 )^0.5 ; %Recursive process
            if dist==0
                dist = 11000;
            end
            min_disp=min(min_disp,dist);
        end
    end
    
    dx=min_disp/2;
    dy=min_disp/2;    
    
    x_edge=[min(X)+0.2:dx:max(X)-0.2];
    y_edge=[min(Y)+0.2:dy:max(Y)-0.2];
    %z_edge=[floor(min(z)):dz:ceil(max(z))];
    [Xgrd{k},Ygrd{k}]=meshgrid(x_edge,y_edge);
    
    F = scatteredInterpolant(X,Y,Z);
    C = scatteredInterpolant(X,Y,stn);
    
    %LEP=griddata(x,y,LEp,Xgrd,Ygrd);  %Interpolate  
    Zeval{k} = F(Xgrd{k},Ygrd{k});
    %Zeval{k}= griddata(X,Y,Z,Xgrd{k},Xgrd{k});
    stneval{k} = C(Xgrd{k},Ygrd{k});
    
     for j = 1:size(Zeval{k},1)     % j is the row --> Y-coord (all entries in the jth row of Ygrd are identical)
         for p = 1:size(Zeval{k},2) % p is the col --> X-coord (all entries in the pth col of Xgrd are identical)
             if sqrt(Xgrd{k}(j,p)^2+Ygrd{k}(j,p)^2) > 1
                 Xgrd{k}(j,p) = nan;
                 Ygrd{k}(j,p) = nan;
                 Zeval{k}(j,p) = nan;
                 stneval{k}(j,p) = nan;
             end
         end
     end
             
    
    %mesh(Xgrd,Ygrd,Zeval);
    
    v = linspace(min(Z),max(Z),256);
    
    figure(3)
    hold on
    surf(Xgrd{k},Ygrd{k},Zeval{k},stneval{k},'linestyle','none','facecolor','interp')
    
    
    figure(4)
    hold on
    surf(Xgrd{k},Ygrd{k},Zeval{k},'linestyle','none','facecolor','interp')
    
     k = k+1;
end;

for i = [3 4];
    figure(i)
    cbar = colorbar();
    set(gcf,'color','w');
    set(gca,'Tickdir','out')
    set(gca,'linewidth',1.5,'fontsize',16)
    xlabel('X','fontsize',16,'rot',0)
    ylabel('Y','fontsize',16,'rot',0)
    zlabel('Z','fontsize',16,'rot',0)
    if i == 1
        ylabel(cbar,'e_e' ,'rot',0,'fontsize',16)
    else
        ylabel(cbar,'Z' ,'rot',0,'fontsize',16)
    end
    view( -61,  9.2644 )
    grid on
end





%{    
    % Find the maximum point
    [maxz,locz]=max(X(:,3));
    %X(~ismember(X(:,4),I),:)=[];
    % Find the distance of every poiint from this max Z point
    % dz is a column vector, same length as X(:,i)
    dz=sqrt((X(:,1)-X(locz,1)).^2+(X(:,2)-X(locz,2)).^2);

    Z = X(:,3);
    C = X(:,4);
    maxmin = [maxmin; i max(C) min(C)];   
    
    if i == 0;
        Z0 = Z;
    end;
    
    %Z = Z - min(Z0);
    
    % Create a projected distance coordinate (projected onto XY plane)
    d = [];
    for p = 1:length(X(:,1));
        if X(p,2) < X(locz,2)
            d(p)=-sqrt((X(p,1)-X(locz,1)).^2+(X(p,2)-X(locz,2)).^2);
        else
            d(p)=sqrt((X(p,1)-X(locz,1)).^2+(X(p,2)-X(locz,2)).^2);
        end;
    end;
    d = d';
    %text(0,maxz,sprintf('\n%d',i))
        
    text(.35*2,1.5*X(end,3),sprintf('   %d\n\n\n',k-1),'HorizontalAlignment','left','color','r','fontsize',14)
    
    %%% If holding, which is necessary if you want a legend:
%        patch([d' fliplr(d')], [ Z' fliplr(Z'-.01/R)],[C' fliplr(C')],'EdgeColor','interp','FaceColor','interp')
%        hold on
    %%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This little block of code spline fits the data.  The hope was that it
% would make the lines "smoother" but it seems to make no difference
% sorted=sortrows([d X(:,3) X(:,4)]);
% d = sorted(:,1);
% Z = sorted(:,2);
% C = sorted(:,3);
% 
% fitx=linspace(min(d),max(d),1000);
% fity = interp1(d,[Z C],fitx,'spline');
% clear d Z C
% d = fitx';% Z = fity(:,1);
% C = fity(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    xdata(k,:) = [d' fliplr(d')];
    zdata(k,:) = [Z' fliplr(Z'-.01/R)];
    cdata(k,:) = [C' fliplr(C')];

        %autoArrangeFigures
        
    stg = [stg; i];
        
    k = k +1;        
end;


loc = find(stages == 160);
minstn160 = min(cdata(loc,:));
cdata(1:loc-1,:) = minstn160;
maxstn = max(cdata(end,:));
%crange = linspace(minstn160,maxstn,9);      
%crange = linspace(.09,.29,9);

mycmap = [0 0 139; 0 191 255; 0 201 87; 192 255 62; 255 215 0 ; 255 97 3; 255 0 0; 255 0 255]/255;

patch(xdata',zdata',cdata','EdgeColor','interp','FaceColor','interp')
axis([-.4 .4 -.01 .35]*2)

%axis equal
%title('Bulge Profiles','fontsize',14)
%ylabel('$\displaylstyle\frac{A-A(-1)}{Y}$','interpreter','latex')
%xlab = xlabel('_                        r/R','fontsize',18)
%L = legend(strsplit(num2str(fliplr(stg'))));
cbr = colorbar;
    if exist('cbr')
          caxis([0.1 .3])
          ylabel(cbr,'e_e_q' ,'rot',0,'fontsize',18)
          set(cbr,'YTick',[.1 .15 .2 .25 .3])
          ylab = get(cbr,'YLabel');
          set(ylab,'Position',get(ylab,'Position') + [3 .07 .0])
          yt=get(cbr,'YTick');
          set(cbr,'YTickLabel',sprintf('%.2f|',yt))
    else
        colorbar
    end;
set(gcf,'color','w')
set(gca,'Tickdir','out')
set(gca,'linewidth',1.5,'fontsize',16)
set(gcf,'Units','Normalized','Outerposition',[0 0 9/11 .5]); 
a = annotation('textarrow',[.64 .71],[.023 .023]);  % X-axis arrow
b = annotation('textarrow',[.07 .07],[.57 .72]);   %Y-axis arrow
C = annotation('textarrow',[.935 .935],[.57 .72]);   %Cbar arrow

%}