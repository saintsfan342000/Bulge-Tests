function zz_DomeProfiles_PaperFigure();

close all


%%%% Where are the files located
    path='F:\Martin_Experiments\Sandia_Anis\Bulge Tests\BT-2_Results\';
    prefix='Section_Profile\Sxn_section4_0-';
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
    X=load(sprintf('%s\\%s%d.txt',path,prefix,last));
    keepers = X(:,5);
    
k = 1;
stg = [];
xdata = [];
zdata = [];
cdata = [];
maxmin = [];
figure

stages = [0:40:160 180 200 218]; 
modstn_stages = [0 40 80 120];

for i = stages
        
    clear X d dz
    % Load the section file
    X=load(sprintf('%s\\%s%d.txt',path,prefix,i));

    X = X(ismember(X(:,5),keepers),:);
    
    X(:,[1 2 3]) = X(:,[1 2 3])/R;
    
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