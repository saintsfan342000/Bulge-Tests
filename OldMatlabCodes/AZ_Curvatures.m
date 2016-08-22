clear
clc
close all

savefig=1;

curdir=pwd;
addpath(sprintf('%s\\MATLAB\\extras',curdir(1:2)));


paths = {'../BT-2_Recalc_Results/','../BT-3_StgZeroDeleted_Results/','../BT-5_Results/R1p5/' '../BT-6_Results/'};
%paths = {paths{4}};
%lstring = {'BT-2','BT-3','BT-4'};

% (1)Stage (2)Pressure (3)Bulge Height (4)Rx (5)Ry (6)ex (7)ey (8)eps_xy

%path = ['../../Bulge Tests/' 'BT-6_Results/'];

for i = 1:length(paths)
    clear loc
    testres = load([paths{i} 'TestResults.dat']);
    R = testres(:,1);
    Z = testres(:,11);
    phi = testres(:,6);
    
    STLP = load( [paths{i} 'STLP.dat'] );
    P = STLP(:,4);
    
    xyrd = load( [paths{i} 'Roll-Trans-Radius.dat'] );
    Rx = xyrd(:,1);
    Ry = xyrd(:,2);
    ex = xyrd(:,3);
    ey = xyrd(:,6);
    %    exy = mean(xyrd(:,[7 8]),2);
    
    %exp = [STLP(:,1)  P  Z  Rx  Ry  ex  ey  exy R];
    
    %fid = fopen([path 'BT-6_ForKelin.dat'],'w');
    %fprintf(fid,'%.0f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f\n',exp');
    %fclose(fid);
    
    %Critical Stages
    if paths{i}(7) ~= '2'
        cstg = [.05 .1 .15 .2 .3];
        loc = [];
        
        L1=[];
        figure(i)
        set(gcf,'color','w')
        plot(Z/3,Rx,'b')
        hold on
        plot(Z/3,Ry,'r')
        title(paths{i}(4:7),'fontsize',18)
        xlabel('h/R','fontsize',18)
        ylabel('\rho','fontsize',18,'rot',0)
        
        for k = 1:length(cstg)
            [~,loc] = min(abs(Z/3 - cstg(k)));
            %loc =(1);
            rdiff = abs(Ry(loc)-Rx(loc))/mean([Ry(loc),Rx(loc)]) * 100;
            L1 = [L1;sprintf('h/R = %.2f: \\rho_x - \\rho_y = %.1f%%',cstg(k),rdiff)];
            plot(Z(loc)/3,Rx(loc),'ko','markerfacecolor','k','markersize',3)
            [cstg(k) loc abs(Rx(loc)-Ry(loc))];
        end
        axis([0,max(Z/3),min(Ry)*.95,max(Ry)*1.05])
        text(.35,30,L1,'verticalalignment','top')
        L2='';
    else    %BT-2
        cstg = [.05 .1 .15 .2 .3];
        loc = [];
        
        L1=[];
        figure(i)
        set(gcf,'color','w')
        plot(Z/3,Rx,'b')
        hold on
        plot(Z/3,Ry,'r')
        title(paths{i}(4:7),'fontsize',18)
        xlabel('h/R','fontsize',18)
        ylabel('\rho','fontsize',18,'rot',0)
        
        for k = 1:length(cstg)+1
            if k == length(cstg)+1
                loc=98;
                kolor = 'r';
                msize=6;
                rdiff = abs(Ry(loc)-Rx(loc))/mean([Ry(loc),Rx(loc)]) * 100;
%                L1 = [L1;sprintf('h/R = %.2f: \\rho_x - \\rho_y = %.1f%%',cstg(k),rdiff)];
                L2 = [sprintf('Wp=3000psi: \\Delta\\rho/\\rho = %.1f%%',rdiff)];
            else
                [~,loc] = min(abs(Z/3 - cstg(k)));
                kolor='k';
                rdiff = abs(Ry(loc)-Rx(loc))/mean([Ry(loc),Rx(loc)]) * 100;
                msize=3;
                L1 = [L1;sprintf('h/R = %.2f: \\Delta\\rho/\\rho = %.1f%%',cstg(k),rdiff)];
            end
            plot(Z(loc)/3,Rx(loc),[kolor 'o'],'markerfacecolor',kolor,'markersize',msize)
        end
        ex(loc)
    end
    
    axis([0,max(Z/3),min(Ry)*.95,max(Ry)*1.05])
    text(.35,30,L1,'verticalalignment','top')
    text(Z(98)/3*1.05,Ry(98)*1.05,L2,'color','r','verticalalignment','bottom')
end


if savefig == 1
    for i = 1:length(paths)
        figure(i)
        %set(gca,'linewidth',1)%,'fontsize',16)
        export_fig([paths{i} 'Roll-Trans-Radius.png'],'-r150','-nocrop')
        close(i)
    end
end
