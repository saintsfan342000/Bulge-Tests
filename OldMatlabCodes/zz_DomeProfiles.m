function zz_DomeProfiles();

close all


%%%% Where are the files located
    path='F:\Martin_Experiments\Sandia_Anis\Bulge Tests\BT-2_Recalc';
%     prefix={'Major-Section\Maj-Sxn_section2_0-';'Minor-Section\Min-Sxn_section3_0-'};
    prefix={'Section_X\Sxn_section0_0-';'Section_Y\Sxn_section1_0-'};
    %prefix={'X-Axis-Section\X-Sxn_section0_0-';'Y-Axis-Section\Y-Sxn_section1_0-'};
    last=218;
%%%% THreshhold is the z-value in the last stage above which we select
%%%% points (the points considered to be at the apex of the bulge
    threshforradius=1.585;
    threshforstn=1.595;
    distfrommaxz = 0.5;
    
% Add matlab extras for plane fit, circle fit
curdir=pwd;
addpath(sprintf('%s\\Matlab\\extras',curdir(1:2)));

% Columns of each section file are:
%(1) CoordX (2) CoordY (3) CoordZ (4)INdex (5)StrainP1 (6)StrainP2

for n=1:2;
    clear X I

    for i=[0:25:last last];
        
        clear X d dz
        % Load the section file
        X=load(sprintf('%s\\%s%d.txt',path,prefix{n},i));
        % Find the maximum point
        [maxz,locz]=max(X(:,3));
        %X(~ismember(X(:,4),I),:)=[];
        % Find the distance of every poiint from this max Z point
        % dz is a column vector, same length as X(:,i)
        dz=sqrt((X(:,1)-X(locz,1)).^2+(X(:,2)-X(locz,2)).^2);
        
        X(dz>distfrommaxz,:)=[];
        
        % Now just create a dist along the section coorinate
        d=sqrt((X(:,1)-X(1,1)).^2+(X(:,2)-X(1,2)).^2);
        
        figure(i+1)
        subplot(1,2,n)
        plot(d,X(:,3),'.')
        %axis equal
        
        autoArrangeFigures
        
        if n==1
            title(sprintf('Major Direction - Stg %d',i))
        elseif n==2
            title('Minor Direction')
        end;
        
    end;
    
end;