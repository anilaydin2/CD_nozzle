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