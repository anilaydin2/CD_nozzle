function [m_dot] = mdot_needed_for_given_thrust(thrust,P,k,h,T_total,R)
    [~, ~ , p_inf, ~] = atmosisa(h); %irtifa kontrol et
    p_exit = p_inf; %Pa
    M_exit = ((2.0/(k-1)).*((P./p_exit).^((k-1)/k)-1)).^(0.5); 
    T_exit = T_total*(1+((k-1)/2.0)*(M_exit^(2.0)))^(-1);
    V_exit = M_exit*(k*R*T_exit)^(0.5);
    m_dot = thrust/V_exit;
end