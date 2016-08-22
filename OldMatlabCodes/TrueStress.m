function [trunew, epeq,k]=TrueStress(P,e1,e2,E,v,R,to);

E=E*1000;   % Since E is input in KSI, we need it in PSI to be compatible with pressure
clear tnew
tnew=10000; % Initialize to erroneously high value

%Initial approximation of thickness
    ta=to*exp(-(e1+e2));

k=1;
while min(ta/tnew,tnew/ta)<.9999999999
    if k~=1;
        ta=tnew;  %Update ta to the previously calculated tnew if we're not on the first iteration
    end;
    tru=P*R/(2*ta); %Approximate true stress based on ta
    e1p=e1-(1-v)*tru/E; %e_plastic strain is e1 minus e1_elastic
    e1p=e1-(1-v)*tru/E;
    e3p=-(e1p+e2p); %Assume plastic incompressibility
    e3=e3p-2   %Then add on the elastic part of e3 (plane stress) to get e3_total
    epeq=sqrt((2/3)*(e1p^2+e2p^2+e3p^2));   %Calcu epeq 
    tnew=to*exp(e3);  %Get new thickness
    trunew=P*R/(2*tnew); %New true stress
    k=k+1;
end;