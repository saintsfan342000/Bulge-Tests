function strainrate();

fp = '../BT-2_Results/Final Computation - Radius0p6-Stn0p3/TestResults.dat';
fp2 = '../BT-2_Results/Final Computation - Radius0p6-Stn0p3/STLP.dat'

D = load(fp);
STLP = load(fp2);
% (2)MeanMajorStn (3)MeanMinorStn

D = D(:,[2 3]);

epeq = sqrt( (2/3)*(D(:,1).^2 + D(:,2).^2 +(-D(:,1)-D(:,2)).^2)  );
depeq = epeq(2:end)-epeq(1:end-1);

de1 = D(2:end,1) - D(1:end-1,1);

dt = STLP(2:end,2)-STLP(1:end-1,2);

subplot(2,1,1)
plot(depeq./dt * 10^4,'b.')
set(gcf,'color','w');
ylabel({'e_e^\bullet        ', '(10^-^4/s)           '}','fontsize',16,'rot',0)

subplot(2,1,2)
plot(de1./dt*10^4,'r.')
ylabel({'e_I^\bullet        ', '(10^-^4/s)           '},'fontsize',16,'rot',0)

xlswrite('../../AA_Paper/StrainRate.xlsx',[depeq./dt de1./dt])