clear
clc
close all
[T_inf, c_inf , p_inf, rho_inf] = atmosisa(1000.0); %irtifa kontrol et
k = 1.4000; % specific 
T_total = 300.0; %K
p_total = 4000000.0; %Pa basınç değişecek
flow_rate_L = 0.0:0.05:30.0; %L/s debi
flow_rate = flow_rate_L/1000.0; %m3/s
R = 286.9; %J/kgK 
r_inlet = (3.0/(8.0*2.0))/39.37;
A_inlet = pi()*r_inlet^2; %m

%direk girilen değişkenler
%m_dot_throat = 0.125;
%rho_inlet = p_total/(R*T_total);
%V_inlet = m_dot_throat/(rho_inlet*A_inlet);
V_inlet = flow_rate./A_inlet; %m/s
M_inlet = V_inlet./((k*R*T_total)^(0.5)); 
T_throat = T_total*(1+((k-1)/2.0))^(-1);
p_throat = p_total*(1+((k-1)/2.0))^(-k/(k-1));
rho_throat = p_throat/(R*T_throat);
v_throat = (k*R*T_throat)^(0.5);
rho_inlet =rho_throat*((k+1)/2.0)^(1.0/(k-1));

A_throat = A_inlet./((((k+1)/2.0)^(-(k+1)/(2.0*(k-1))))*((1+((k-1)/2.0).*(M_inlet.^2.0)).^((k+1)/((k+1)/(2.0*(k-1)))))./M_inlet);
m_dot_throat = A_throat*p_total/sqrt(T_total)*sqrt(k/R)*(1+(k-1)/2.0)^(-(k+1)/(2*(k-1)));

D_throat = (4*A_throat./pi()).^(0.5);

p_exit = p_inf;
M_exit = ((2.0/(k-1))*((p_total/p_exit)^((k-1)/k)-1))^(0.5);

Ae_Ai_ratio = (M_inlet./M_exit).*(((1+((k-1)/2.0).*(M_exit^2.0))./(1+((k-1)/2.0).*(M_inlet.^2.0))).^((k+1)/(2.0*(k-1))));

A_exit_1 = Ae_Ai_ratio.*A_inlet;
A_exit = (A_throat./M_exit).*((1+((k-1)/2.0)*M_exit^2.0)/((k+1)/2.0))^((k+1)/(2.0*(k-1)));
D_exit = (4*A_exit./pi()).^(0.5);

rho_exit = rho_inlet*(1+((k-1)/2.0)*M_exit^2.0)^(-1.0/(k-1));
T_exit = T_total*(1+((k-1)/2.0)*(M_exit^(2.0)))^(-1);
p_exit = p_total*(1+((k-1)/2.0)*(M_exit^(2.0)))^(-k/(k-1));
V_exit = ((2.0*k/(k-1))*R*T_total*(1-((p_exit/p_total)^((k-1)/k)))+V_inlet).^(0.5);
m_dot_exit = A_exit*p_total/sqrt(T_total)*sqrt(k/R)*M_exit*(1+(k-1)*M_exit^(2.0)/2.0)^(-(k+1)/(2*(k-1)));

thrust = m_dot_exit.*V_exit + (p_exit-p_inf).*A_exit; 

figure, clf
plot(m_dot_exit,thrust)
title('Kütlesel Debi vs İtki')
xlabel('Kütlesel Debi [kg/s]') 
ylabel('İtki [N]')
grid on
grid minor

figure, clf
plot(flow_rate_L, m_dot_throat)
title('Hacimsel Debi vs Kütlesel Debi')
xlabel('Hacimsel Debi [L/s]') 
ylabel('Kütlesel Debi [kg/s]')
grid on
grid minor

figure, clf
plot(m_dot_throat, D_throat)
title('D_t vs Kütlesel Debi')
xlabel('Kütlesel Debi [kg/s]') 
ylabel('D_t [m]')
grid on
grid minor

figure, clf
plot(m_dot_throat, D_exit)
title('D_e vs Kütlesel Debi')
xlabel('Kütlesel Debi [kg/s]') 
ylabel('D_e [m]')
grid on
grid minor

