%function MiscDataCompilation()

clear
clc
close all

curdir=pwd;
addpath(sprintf('%s\\MATLAB\\extras',curdir(1:2)));

savefig = 0;

paths = {'BT-2_Recalc_Results/','BT-3_StgZeroDeleted_Results/','BT-5_Results/','BT-6_Results/'};
lstring = {'BT-2','BT-3','BT-5','BT-6'};
    
c={[0 0 0],[0 0 1],[1 0 0],[0 201 87]/355,[255 130 171]/255,[156 102 31]/255,[139 87 66]/255,[255 165 79]/255,[178 58 238]/255,[0 .95 0],[0 220 255]/255};

to = .04;
Ro = 3;

%topath = '../../Bulge Tests/BT-4_Results/';

for i = 1:length(paths)

path = ['../../Bulge Tests/' paths{i}];
    
STLP = load( [path 'STLP.dat'] );
    if lstring{i} == 'BT-3'
        L = (STLP(:,3)-STLP(1,3))*10; P = STLP(:,4); time = STLP(:,2);
    else
        L = STLP(:,3)-STLP(1,3); P = STLP(:,4); time = STLP(:,2);
    end
    
keeps = true(length(L),1);    
if lstring{i}  == 'BT-5'
    pts = find(L>4.85 & P<1417);
    keeps(pts) = false;
    L = L(keeps);    P = P(keeps);    time = time(keeps);
end
    
    
D = load([path 'TestResults.dat']);
%(1)NomSts (2)MeanMajorStn (3)MeanMinorStn (4)Mean emaj/emin (5)std(emaj/emin (6)Mean principal direction
%(7)Mean ex (8)Mean ey (9)Mean ey/ex (10)std(ey/ex) (11) Dome Heigh[in] (12)Thickness (13)Change-Thickness

    R = D(keeps,1) ; Z = D(keeps,11) ; t = D(keeps,12) ; dt = D(keeps,13) ; ex = D(keeps,7) ; ey = D(keeps,8) ; phi=D(keeps,6);
    e1 = D(keeps,2) ;  e2 = D(keeps,3) ;
    
B = load([path 'Case3.dat']);
    e1p = B(keeps,2); e2p = B(keeps,3);

a=Z(isnan(e1p)== 0);a(1)

if i == 1
    Zinit = Z(1);
else
    Zdif = Z(1) - Zinit;
    Z = Z - Zdif;
end


% Plastic strain points
interppoints = [0.02;0.04;0.06;0.08];
e1e2p = interp1(e1p(isnan(e1p)==0),e2p(isnan(e1p)==0),interppoints);
E = polyfit(interppoints,e1e2p,1); 
x = linspace(min(interppoints),max(interppoints),50);
y = E(1)*x + E(2);


figure(10*i)
plot(L,P,'linewidth',1.5)
title(lstring{i},'fontsize',14)
xlabel('LVDT (in)','fontsize',14)
ylabel({'P','(psi)'},'rot',0,'fontsize',14)
set(gca,'linewidth',1.5)
set(gca,'Tickdir','out')
axis([0 6 0 1600])
set(gcf,'color','w')
if savefig == 1
    print(gcf,'-dpdf',['../' paths{i} 'P-L'])
    close
end


figure(1)
subplot(2,1,1)
plot(Z,ex,'color',c{i})
hold on
title('Roll.Total Stn','fontsize',14)
xlabel('Bulge Height, h [inches]','fontsize',14)
ylabel('e_x_''   ','fontsize',14,'rot',0)
set(gca,'linewidth',1.5)
set(gca,'Tickdir','out')
subplot(2,1,2)
plot(Z,ey,'color',c{i})
hold on
title('Transv.Total Stn','fontsize',14)
xlabel('Bulge Height, h [inches]','fontsize',14)
ylabel('e_y_''   ','fontsize',14,'rot',0)
set(gca,'linewidth',1.5)
set(gca,'Tickdir','out')

figure(2)
subplot(2,1,1)
plot(Z(isnan(e1p)==0),e1p(isnan(e1p)==0),'color',c{i})
hold on
title('Roll. Plastic Stn')
xlabel('Bulge Height, h [inches]','fontsize',14)
ylabel('e^p_x_''   ','fontsize',14,'rot',0)
set(gca,'linewidth',1.5)
set(gca,'Tickdir','out')
subplot(2,1,2)
plot(Z(isnan(e1p)==0),e2p(isnan(e1p)==0),'color',c{i})
hold on
title('Transv. Plastic Stn')
xlabel('Bulge Height, h [inches]','fontsize',14)
ylabel('e^p_y_''   ','fontsize',14,'rot',0)
set(gca,'linewidth',1.5)
set(gca,'Tickdir','out')

figure(3)
plot(Z,phi,'color',c{i})
hold on
title('Principal Angle')
xlabel('Bulge Height, h [inches]','fontsize',14)
ylabel('\phi^o','fontsize',18,'rot',0)
set(gca,'linewidth',1.5)
set(gca,'Tickdir','out')

figure(4)
set(gcf,'color','w')
yo(i) = plot(x,y,'color',c{i});
hold on
plot(interppoints,e1e2p,'o','MarkerFaceColor',c{i},'MarkerEdgeColor','none')
title('Plastic Strain Ratio')
axis([0 .1 0 .1])
%set(gca,'Xlim',[0 .1])
text(.08,e1e2p(end),sprintf('   %.3f',E(1)))
xlabel('e^p_x','fontsize',14)
ylabel('e^p_y    ','rot',0,'fontsize',14)
set(gca,'linewidth',1.5)
set(gca,'Tickdir','out')


figure(5)
plot(Z,P,'color',c{i});
hold on
title('Pressure-Height')
xlabel('Bulge Height, h [inches]','fontsize',14)
ylabel({'P','(psi)'},'rot',0,'fontsize',14)
set(gca,'linewidth',1.5)
set(gca,'Tickdir','out')

figure(6)
plot(L,P,'color',c{i});
hold on
title('Pressure-Height')
xlabel('LVDT Disp. [inches]','fontsize',14)
ylabel({'P','(psi)'},'rot',0,'fontsize',14)
set(gca,'linewidth',1.5)
set(gca,'Tickdir','out')
axis([0 6 0 1600])

end

for i = 1:2
    figure(i)
    subplot(2,1,1)
    l = legend(lstring);
    set(l,'location','northwest')
    subplot(2,1,2)
    l = legend(lstring);
    set(l,'location','northwest')
    
    set(gcf,'PaperPositionMode','auto')
    set(gcf, 'Position', [0 0 1.2*500 1.2*900])
end
for i = [3,5,6];
    figure(i)
    l = legend(lstring);
    set(l,'location','northwest')
end
for i = 4;
    figure(i)
    l = legend(yo,lstring);
    set(l,'location','northwest')
end

for i = 1:6
    set(i,'color','w')
end

pause(5)

if savefig ==1
    titles = {'TotalStn.png','PlasticStn.png','PrinAngle.png','StnRatio.png','PressHt.png','Press-LVDT.png'};
    prefix = '../../Bulge Tests/';
    for i = 1:6;
        figure(i)
        fname = [prefix titles{i}];
        export_fig(fname,'-r200','-nocrop')
        close(i)
    end
end


%{
To polyfit BT-4 up to just Z = 1:
Z=Z(isnan(e1p)==0);
e1p = e1p(isnan(e1p)==0);
e2p = e2p(isnan(e2p)==0);
%a = find(Z>1); a = a(1);
plot(e1p(1:a),e2p(1:a))
%polyfit(e1p(1:a),e2p(1:a),1)
%}


% figure
% set(gcf,'color','w')
% plot(dt/to,P,'b')
% hold on
% plot(dts/to,P,'r')
% set(gca,'xlim',[0 .3])
% ylabel(sprintf('P\n(psi)'),'fontsize',16,'rot',0)
% xlabel('\Deltat / t_o        ','rot',0,'fontsize',16)
% l=legend('Correct','Elastic neglected');
% set(l,'location','northwest','orientation','vertical')
% %export_fig('Pressure_vs_DeltaT-two.png','-nocrop','-r250')

%{

figure
set(gcf,'color','w')
plot(Z/Ro,R/Ro,'linewidth',2);
xlabel('h / R_o','fontsize',16)
ylabel('\rho/R_o   ','fontsize',16,'rot',0)
export_fig('Rho_vs_h.png','-nocrop','-r250')


F=1;
figure(F)
set(gcf,'color','w')
plot(ex,P,'r')
hold on
plot(ey,P,'b')
ylabel('Pressure (psi)','fontsize',14)
xlabel('Strain')
l=legend('e_x','e_y');
set(l,'location','northwest','orientation','vertical')
%export_fig('Pressure-vs-Strain.png','-nocrop','-r250')
close







F = 1;

figure(F)
set(gcf,'color','w')
plot(Z,P,'linewidth',2);
xlabel('h - Bulge Height (in)','fontsize',14)
ylabel('Pressure (psi)','fontsize',14)
export_fig('Pressure vs Bulge Height','-nocrop','-r250')
close
F = F+1;

figure(F)
set(gcf,'color','w')
plot(ex,P,'r')
hold on
plot(ey,P,'b')
ylabel('Pressure (psi)','fontsize',14)
xlabel('Strain')
l=legend('e_x','e_y');
set(l,'location','northwest','orientation','vertical')
export_fig('Pressure-vs-Strain.png','-nocrop','-r250')





%%%%%%%%%%%%%%%%% First Set of Rought Figs %%%%%%%%%%%%%%%%%%%%%
F = 1;

% P, LVDT, Z, R, sig, e_a_v, dt/to vs time
figure(F)
set(gcf,'color','w')
subplot(3,1,1)
plot(time, P,'b','linewidth',2)
title(' Various Parameters vs Time','fontsize',16)
ylabel('P','rot',0,'fontsize',20)
subplot(3,1,2)
plot(time,Z,'r','linewidth',2)
ylabel('h','rot',0,'fontsize',20)
subplot(3,1,3)
plot(time,L,'color',[.65 .65 .65],'linewidth',2)
ylabel('LVDT','fontsize',20)
xlabel('TIME','fontsize',20)
set(gcf,'Position',[100 100 500 900])
export_fig('All_vs_Time_1.png')
close

F = F+1;

% P, LVDT, Z, R, e_a_v, dt/to vs time
figure(F)
set(gcf,'color','w')
subplot(3,1,1)
plot(time, R,'b','linewidth',2)
title(' Various Parameters vs Time','fontsize',16)
ylabel('R','rot',0,'fontsize',20)
subplot(3,1,2)
plot(time,mean([ex ey],2),'r','linewidth',2)
ylabel('\epsilon_a_v','rot',0,'fontsize',20)
subplot(3,1,3)
plot(time,dt./to,'color',[.65 .65 .65],'linewidth',2)
ylabel('\Deltat / t_o','fontsize',20)
xlabel('TIME','fontsize',20)
set(gcf,'Position',[100 100 500 900])
export_fig('All_vs_Time_2.png')
close

F = F+1;

figure(F)
set(gcf,'color','w')
subplot (1,3,1)
plot(L,P,'linewidth',2);
set(gca,'Xlim',[-.050 5])
title('Pressure vs LVDT Disp','fontsize',18)
subplot (1,3,2)
plot(Z,P,'linewidth',2);
title('Pressure vs Bulge Height','fontsize',18)
subplot (1,3,3)
plot(dt./to,P,'linewidth',2);
set(gca,'Xlim',[-.0001 .012]/to)
title('Pressure vs \Deltat/t_o','fontsize',18)
set(gcf,'Position',[100 500 1600 450])
export_fig('Pressure.png')
close
F = F + 1;

figure(F)
set(gcf,'color','w')
subplot (1,3,1)
plot(L,sig,'linewidth',2);
set(gca,'Xlim',[-.050 5])
title('\Sigma vs LVDT Disp','fontsize',18)
subplot (1,3,2)
plot(Z,sig,'linewidth',2);
title('\Sigma  vs Bulge Height','fontsize',18)
subplot (1,3,3)
plot(dt./to,sig,'linewidth',2);
set(gca,'Xlim',[-.0001 .012]/to)
title('\Sigma vs \Deltat/t_o','fontsize',18)
set(gcf,'Position',[100 500 1600 450])
export_fig('Sigma.png')
close
F = F + 1;

figure(F)
set(gcf,'color','w')
subplot (1,2,1)
plot(mean([ex ey],2),P,'linewidth',2);
%set(gca,'Xlim',[-.050 5])
title('Pressure vs \epsilon_a_v','fontsize',18)
subplot (1,2,2)
plot(mean([ex ey],2),sig,'linewidth',2);
%set(gca,'Xlim',[-.050 5])
title('\Sigma vs \epsilon_a_v','fontsize',18)
set(gcf,'Position',[100 500 1100 450])
export_fig('Strain.png')
close
F = F + 1;


%}