
function Rank4_RungeKutta
h=0.002;
t=0:h:40;
n=length(t); 
k3=2.5;
E_react=1;
Y(1,1)=0;
Y(2,1)=10;

% c(S) and c(ES) solving
for k=1:n-1  
    z1=concentration_calculation(t(k),Y(1:2,k));  
    z2=concentration_calculation(t(k)+h/2,Y(1:2,k)+z1*h/2);
    z3=concentration_calculation(t(k)+h/2,Y(1:2,k)+z2*h/2);
    z4=concentration_calculation(t(k)+h,Y(1:2,k)+z3*h);
    Y(1:2,k+1)=Y(1:2,k)+h*(z1+2*z2+2*z3+z4)/6; 
end
S=Y(2,:); 
ES=Y(1,:);
 
% c(E) and c(P) solving
E_react=E_react-ES;
for i=1:1:n  
    P(i)=sum(ES(1:i))*h*k3; 
end

% v(P) solving
Vp=k3.*ES; 
 
% plot map
figure(1);
plot(t,E_react,'LineWidth',3); 
legend('E');  
title('Concentration Variety of E');
xlabel('time /s');
ylabel('concentration_E /uMol');

figure(2);
plot(t,S,'LineWidth',3);
legend('S');  
title('Concentration Variety of S');
xlabel('time /s');
ylabel('concentration_S /uMol');

figure(3);
plot(t,ES,'LineWidth',3); 
legend('ES');  
title('Concentration Variety of ES');
xlabel('time /s');
ylabel('concentration_ES /uMol');

figure(4);
plot(t,P,'LineWidth',3);
legend('P');  
title('Concentration Variety of P');
xlabel('time /s');
ylabel('concentration_P /uMol');

%figure(5); 
%plot(S,Vp,'LineWidth',3);
%title('Vp change caused by concentration of S');
%xlabel('Sconcentration/uMol');
%ylabel('uMol/s');

end

