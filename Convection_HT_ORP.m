function [Q_Convection] = Convection_HT_ORP(area, perimeter, TR, T_amb, WNDSPD, dyn_visc_air, rho_air)
%This function calculates the convection heat flux from the ORP

kin_visc_a = dyn_visc_air./rho_air; 
L_c = area/perimeter; 
alpha_a = 2.2*10^-5;
lamda_a = 2.6*10^-2;

Re_L = L_c.*WNDSPD/kin_visc_a; 
Pr_a = kin_visc_a/alpha_a;

    if Re_L < (3*10^5) %then the flow is laminar and use:
           Nu_L = 0.628.*(Re_L.^0.5).*(Pr_a^(1/3)); 
    elseif Re_L > (5*10^5) %then the flow is turbulent and use: 
           Nu_L = 0.035.*(Re_L.^0.8).*(Pr_a^(1/3));
    else   %Average the values sherwood values for laminar and turbulent
           Nu_L = ((0.628.*(Re_L.^0.5).*(Pr_a^(1/3)))+(0.035.*(Re_L.^0.8).*(Pr_a^(1/3))))/2;
    end 

%Caclculate the convection coefficient given the Nusselt number
h_conv = Nu_L*lamda_a/L_c; 

%Calculate the convective heat transfer
Q_Convection = h_conv.*(T_amb - TR);

% Pretty much straight from Yadala and Cremaschi, 2016. Except for the added
% correlation for laminar flow and averaging the two if in the transisition
% period.

end