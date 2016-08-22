function zz_DomeProfiles_3D_Z();

close all


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
    
    ArI = d(:,1);
    Arj = d(:,2);
    X = d(:,3)+d(:,6);
    Y = d(:,4)+d(:,7);
    Z = d(:,5)+d(:,8);
    stn = zeros(length(X),1);
    
%     for j = 1:length(X)
%         F = [d(j,9) d(j,10);d(j,11) d(j,12)];
%         U = sqrtm(F'*F);
%         evs = eig(U);
%         L1 = evs(1);
%         L2 = evs(2);
%         e1 = log(L1);
%         e2 = log(L2);
%         stn(j,1) = sqrt( (2/3) * (e1^2 + e2^2 + (-e1-e2)^2) ); 
%     end;
    
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
%     C = scatteredInterpolant(X,Y,stn);
    
    %LEP=griddata(x,y,LEp,Xgrd,Ygrd);  %Interpolate  
    Zeval{k} = F(Xgrd{k},Ygrd{k});
%     stneval{k} = C(Xgrd{k},Ygrd{k});
    
    %mesh(Xgrd,Ygrd,Zeval);
    
    v = linspace(min(Z),max(Z),256);
    
    figure(1)
    hold on
%     surf(Xgrd{k},Ygrd{k},Zeval{k},stneval{k},'linestyle','none','facecolor','interp')
    
    
    figure(2)
    hold on
    surf(Xgrd{k},Ygrd{k},Zeval{k},'linestyle','none','facecolor','interp')
    
%          for j = 1:size(Zeval{k},1)     % j is the row --> Y-coord (all entries in the jth row of Ygrd are identical)
%          for p = 1:size(Zeval{k},2) % p is the col --> X-coord (all entries in the pth col of Xgrd are identical)
%              if sqrt(Xgrd{k}(j,p)^2+Ygrd{k}(j,p)^2) > 1
%                  Xgrd{k}(j,p) = nan;
%                  Ygrd{k}(j,p) = nan;
%                  Zeval{k}(j,p) = nan;
%                  stneval{k}(j,p) = nan;
%              end
%          end
%      end
%              
%     
%     %mesh(Xgrd,Ygrd,Zeval);
%     
%     v = linspace(min(Z),max(Z),256);
%     
    figure(3)
    hold on
%     surf(Xgrd{k},Ygrd{k},Zeval{k},stneval{k},'linestyle','none','facecolor','interp')
    
    
    figure(4)
    hold on
%     surf(Xgrd{k},Ygrd{k},Zeval{k},'linestyle','none','facecolor','interp') 
        
    k = k+1;
end;
   
for i = 2;
    figure(i)
    colormap jet
    set(gcf,'position',[199         149        1330         614])
    cbar = colorbar();
    set(gcf,'color','w');
    set(gca,'Tickdir','out')
    set(gca,'Xtick',[-2 -1 0 1 2])
    set(gca,'Ytick',[-2 -1 0 1 2])
    set(gca,'Ztick',[0.6 0.8 1.0 1.2 1.4 1.6])
    set(gca,'linewidth',1.5,'fontsize',16)
    xlabel('X','fontsize',16,'rot',0)
    ylabel('Y','fontsize',16,'rot',0)
    zlabel('Z','fontsize',16,'rot',0)
    if any( i == [1 3])
        ylabel(cbar,'e_e' ,'rot',0,'fontsize',16)
    else
        %ylabel(cbar,'Z' ,'rot',0,'fontsize',16)
        axis([-2 2 -2 2 0.5 1.6])
        set(cbar,'YTick',[0.6 0.8 1.0 1.2 1.4 1.6])
        set(cbar,'Position',[0.939423558897243   0.202182410423453   0.020050125313283   0.671465798045603])
        %set(cbar,'position',
    end
    if i == 2 || i ==4
        caxis([0.6 1.6])
    end
    view( -61,  9.2644 )
    grid on
    set(gca,'GridLineStyle','-','GridColor',[0.3 0.3 0.3])
    %zlabel('Z','fontsize',16,'rot',0)
end

close 1 3 4

