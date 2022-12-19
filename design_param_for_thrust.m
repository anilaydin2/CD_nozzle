clear
clc

%Constants
k = 1.4000; % specific 
T_total = 300.0; %K
h = 1000.0; % m
R = 286.9; %J/kgK 

% User defined values
p_total = 4e6("Inlet Pressure [Pa]: "); %Pa
r_inlet = 0,00748("Inlet radius [m]: ");
thrust = 25("Wanted thrust [N]: ");

%Calculations
A_inlet = pi*r_inlet^2; %m
mdot = mdot_needed_for_given_thrust(thrust,p_total,k,h,T_total,R) %kg/s
[D_exit, D_throat] = nozzle_dia(mdot,p_total,T_total,k,h,R,r_inlet) %m


%% Functions

function [m_dot] = mdot_needed_for_given_thrust(thrust,P,k,h,T_total,R)
    [~, ~ , p_inf, ~] = atmosisa(h); %irtifa kontrol et
    p_exit = p_inf; %Pa
    M_exit = ((2.0/(k-1)).*((P./p_exit).^((k-1)/k)-1)).^(0.5); 
    T_exit = T_total*(1+((k-1)/2.0)*(M_exit^(2.0)))^(-1);
    V_exit = M_exit*(k*R*T_exit)^(0.5);
    m_dot = thrust/V_exit;
end

function [D_exit, D_throat] = nozzle_dia(m_dot,P,T_total,k,h,R,r_inlet)
    [~, ~ , p_inf, ~] = atmosisa(h); %irtifa kontrol et
    A_inlet = pi*r_inlet^2; %m

    %direk girilen değişkenler
    rho_inlet = P./(R*T_total);
    V_inlet = m_dot./(rho_inlet.*A_inlet);
    M_inlet = V_inlet./((k*R*T_total)^(0.5)); 
    
    A_throat = A_inlet./((((k+1)/2.0)^(-(k+1)/(2.0*(k-1))))*((1+((k-1)/2.0).*(M_inlet.^2.0)).^((k+1)/((k+1)/(2.0*(k-1)))))./M_inlet);
    D_throat = (4*A_throat./pi()).^(0.5);
    
    p_exit = p_inf;
    M_exit = ((2.0/(k-1)).*((P./p_exit).^((k-1)/k)-1)).^(0.5);
    
    A_exit = (A_throat./M_exit).*((1+((k-1)/2.0)*M_exit^2.0)/((k+1)/2.0))^((k+1)/(2.0*(k-1)));
    D_exit = (4*A_exit./pi()).^(0.5);

end
