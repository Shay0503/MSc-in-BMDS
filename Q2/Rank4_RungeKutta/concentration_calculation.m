function concentration_calculation=f(t,Y)  
k2=10;
k3=2.5;
E0=1;
k1=5/3;
S=Y(2,1);
E=Y(1,1);  
% concentration calculation
f1=k1*E0*S-(k2+k3+k1*S)*E;  
f2=(k1*S+k2)*E-k1*E0*S; 
concentration_calculation=[f1;f2];
end