clc
clear 
% k=20;
Ha=0.01; 
Du=0; 
EI=0.01;
Da=10;
fi=0.0;
alpha_E=0.3;
gta=1;
bta=0.0;
alpha_e=0.0;
alpha_s=0.0;
efs_ratio=80/(6.95*10^(-10));
sigma_ratio=25000/0.1;
row_ratio=5.18/1.025;
S=0;
P=1;
i=0;
for k=10:0.1:40
  i=i+1;
 sum=0;
for n = 1: 300  
    equation = @(lamda) (cos(lamda) - bta * lamda * sin(lamda));
    lamda_guess = (2 * n - 1) * pi / 2;
    lamda = Newton_Raphson(equation, lamda_guess);
  gta_a=gta*(1+bta*k*(sinh(gta)/gta)); 
  alpha_1=alpha_e/((1+alpha_e*alpha_s)^2+alpha_e^2);
  alpha_2=(1+alpha_e*alpha_s)/((1+alpha_e*alpha_s)^2+alpha_e^2);
  efs_r=1+((2*(efs_ratio-1)*fi)/(((efs_ratio+2)-(efs_ratio-1)*fi)));
  row_r=(1-fi)+(row_ratio)*fi;
  mu_r=1/((1-fi)^2.5);
  sigma_r=1+((3*(sigma_ratio-1)*fi)/((sigma_ratio+2)-(sigma_ratio-1)*fi));
 
   c2=sqrt((sigma_r/mu_r)*alpha_2*Ha^2+(1/Da));
  
    c4=((gta_a*k*sqrt(efs_r))/(k^2-c2^2*efs_r))*(cosh(c2)*tanh(k/sqrt(efs_r)));
     c41=-(P*c4)/(mu_r*c2^2*cosh(c2)+mu_r*bta*c2^3*sinh(c2));
     c42=-(c4*alpha_1*sigma_r*Ha^2*S)/(mu_r*c2^2*cosh(c2)+mu_r*bta*c2^3*sinh(c2));
     c43=-(c4*k^2*efs_r*gta_a)/((k^2-c2^2*efs_r)*(cosh(c2)+bta*c2*sinh(c2)));
     c44=-(c4*bta*k^3*sqrt(efs_r)*gta_a*tanh(k/sqrt(efs_r)))/((k^2-c2^2*efs_r)*(cosh(c2)+bta*c2*sinh(c2)));
     c5=(c2*efs_r*gta_a*sinh(c2))/((k^2-c2^2*efs_r));
     c51=(P*c5)/(mu_r*c2^2*cosh(c2)+mu_r*bta*c2^3*sinh(c2));
     c52=(c5*alpha_1*sigma_r*Ha^2*S)/(mu_r*c2^2*cosh(c2)+mu_r*bta*c2^3*sinh(c2));
     c53=(c5*k^2*efs_r*gta_a)/((k^2-c2^2*efs_r)*(cosh(c2)+bta*c2*sinh(c2)));
     c54=(c5*bta*k^3*sqrt(efs_r)*gta_a*tanh(k/sqrt(efs_r)))/((k^2-c2^2*efs_r)*(cosh(c2)+bta*c2*sinh(c2)));
     c6=(gta_a*sqrt(efs_r)*tanh(k/sqrt(efs_r)))/(c2^2*k);
     c61=(P*c6)/mu_r;
     c62=(c6*alpha_1*sigma_r*Ha^2*S)/mu_r;
     c7=((k^2*efs_r*gta_a^2)/(2*mu_r*(k^2-c2^2*efs_r)))*(1/(cosh(k/sqrt(efs_r))^2));
     c8=(0.5*k*efs_r^1.5*gta_a^2*tanh(k/sqrt(efs_r)))/(mu_r*(k^2-c2^2*efs_r));
   E_s=(c41+c51+c61)/(alpha_E+alpha_E*Du-c42-c43-c44-c52-c53-c54-c62-c7-c8);
c3=row_r/mu_r;
c1=(P)/mu_r;    
S1=(EI*c2^2+c3)/(2*c3*EI);
S2=(c2^2+lamda^2)/(2*c3*EI);
mu_n=sqrt(2*S2-S1^2);    
D=(1/(cosh(c2)+bta*c2*sinh(c2)))*(-c1/c2^2); 
A1=(D*lamda*cosh(c2)*sin(lamda))/(lamda^2+c2^2)+D*(c2/(lamda^2+c2^2))*(sinh(c2)*cos(lamda));
end
v1 =(D/c2)*sinh(c2)+(c1/c2^2);
R1=v1;
Q_in = R1;
efficiency(i,1)=(E_s^2*alpha_E*(1+Du)*k^2)/(4*P*Q_in);
y(i,1)=k;
end
hold on, plot(y,efficiency)