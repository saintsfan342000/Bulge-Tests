%function [rad,exp]=SectionCircleFit_5();

%%%% Where are the files located
    path='E:\Martin_Experiments\Sandia_Anis\Bulge Tests\BT-6_Results';
    savepath=path;
%     prefix={'Major-Section\Maj-Sxn_section2_0-';'Minor-Section\Min-Sxn_section3_0-'};
    prefix={'Section_X\Sxn_section0_0-';'Section_Y\Sxn_section1_0-'};
    %savepath='F:\Martin_Experiments\Sandia_Anis\Bulge Tests\BT-2_FS75SS8\Radius0p3-Stn0p3'
    last=412;
%%%% THreshhold is the z-value in the last stage above which we select
%%%% points (the points considered to be at the apex of the bulge
    distfrommaxz = 1.2/2;
    distfrommaxzforstn = 0.3;
    
% Add matlab extras for plane fit, circle fit
curdir=pwd;
addpath(sprintf('%s\\Matlab\\extras',curdir(1:2)));

% Columns of each section file are:
%(1) CoordX (2) CoordY (3) CoordZ (4)INdex (5)StrainP1 (6)StrainP2
c=[];
for n=1:2;
    clear X I
    k=1;
for i=0:last;
    i
    clear X en V P E C Rt Par
    % Load the section file
    X=load(sprintf('%s\\%s%d.txt',path,prefix{n},i));
    
    %X(:,3) = X(:,3) + 0.91885186/25.4;
    
    Xcopy=X;
    X=load(sprintf('%s\\%s%d.txt',path,prefix{n},i));
    % Find the maximum point
    [~,locz]=max(X(:,3));
    % Find the distance of every poiint from this max Z point
    % dz is a column vector, same length as X(:,i)
    dz=sqrt((X(:,1)-X(locz,1)).^2+(X(:,2)-X(locz,2)).^2);
    X(dz>distfrommaxz,:)=[];    
    Xcopy(dz>distfrommaxzforstn,:)=[];
    a = length(Xcopy(:,1));
    Xcopy(Xcopy(:,5)>nanmean(Xcopy(:,5))+nanstd(Xcopy(:,5)) | Xcopy(:,5)<nanmean(Xcopy(:,5))-nanstd(Xcopy(:,5)),:)=[]; 
    b = length(Xcopy(:,1));
    e1(k,n)=nanmean(Xcopy(:,5));
    e2(k,n)=nanmean(Xcopy(:,6));
    e12(k,n)=nanmean(Xcopy(:,7));
    % Now fit a plane to these remaining points
        %en is the normal vector, a column vector; 
        %V should be 3 x 2, each column being one of the two ON basis vectors of the plane. 
        %P is a point on the plane
        [en,V,P] = affine_fit(X(:,1:3));
        E=[V en];
        O=[1 0 0;0 1 0;0 0 1];
        for p=1:3;
            for q=1:3;
                C(p,q)=dot(E(:,p),O(:,q));
            end;
        end;
        % Now project each point onto the plane and then transform to plane coordinate system carry out steps necessary to circle fit
    for j=1:length(X);
        clear R
        % R is the cartesian system coordinates of the projection of each
        % point onto the plane
        X(j,1:3)=(X(j,1:3)-P);
        R=(X(j,1:3)'-dot(X(j,1:3),en)*en);
        %R=(X(j,1:3)'-dot(X(j,1:3),en)*en)-P';  %<<--- I expected this to be
        % the procdure, but this produces a non-zero 3-component.  The
        % procedure used yields the same 1 and 2 components and this one,
        % but a zero 3 component
        % Rt is the components of the projected point in the [(V(:,1), V(:,2), en] coordinate system of the plane
        Rt(:,j)=(C*R)'; %Here this are the components projected into the plane
        % Delete the ~0 third row
    end;
    Rt(3,:)=[];
    Par = CircleFitByPratt(Rt');
    rad(k,n)=Par(3);
    k=k+1;
    c = [c ; a b];
end;
end;


exp=[rad e1 e2 e12];

%fid = fopen(sprintf('%s\\Roll-Trans-Radius.dat',savepath),'w');
%fprintf(fid,'%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f\n',exp');
%fprintf(fid,'%.8f,%.8f,%.8f,%.8f,%.8f,%.8f  \n',exp');
%fclose(fid)

figure
%subplot(2,1,1)
plot(rad(:,1),'b')
hold on 
plot(rad(:,2),'r')
legend('Roll','Trans','location','northeast')
set(gcf,'color','w')
axis([0 400 0 35])
%export_fig([savepath '/R1-R2.png'],'-nocrop','-r250')