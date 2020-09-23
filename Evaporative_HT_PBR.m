function [M_Evap, Q_Evap] = Evaporative_HT_PBR(depth, width, TR, T_amb, RH, WNDSPD, dyn_visc_air, rho_air) 
%This function calculates teh evaporation rate of the pond as well as the
%cooling effect due to that evaporation

kin_visc_a = dyn_visc_air./rho_air; 
L_c = (depth*width)/(2*depth+2*width); 
D_w_a = 2.4*10^-5; 
M_water = 0.018; %kg/mol
R = 8.314; %Universal gas constant Pa*m3/mol*K
hfg_water = 2.45*10^6; % j/kg

%Calculate the Reynold's Number
Re_L = L_c.*WNDSPD/kin_visc_a; 
Sch_L = kin_visc_a./D_w_a;
     
    if Re_L < (3*10^5) %then the flow is laminar and use:
           Sh_L = 0.628.*(Re_L.^0.5).*(Sch_L.^(1/3)); 
    elseif Re_L > (5*10^5) %then the flow is turbulent and use: 
           Sh_L = 0.035.*(Re_L.^0.8).*(Sch_L.^(1/3));
    else   %Average the values sherwood values for laminar and turbulent
           Sh_L = ((0.628.*(Re_L.^0.5).*(Sch_L.^(1/3)))+(0.035.*(Re_L.^0.8).*(Sch_L.^(1/3))))/2;
    end 
    
%Now with the sherwood number we can calculate mass transfer coefficient, K
K = Sh_L.*D_w_a./L_c; 

%Calculate the saturated vapor pressures at T_amb and TR
P_w = 3385.5.*exp(-8.0929 + 0.97608.*(TR + 42.607 - 273.15).^0.5); 
P_a = 3385.5.*exp(-8.0929 + 0.97608.*(T_amb + 42.607 - 273.15).^0.5); 

%Now we can calculate the evaporation rate in kg/m2*s
M_Evap = K.*((P_w./TR)-(((RH./100).*P_a)./T_amb))*M_water/R; %Evaporation rate in kg/s*m2 
% play with the 1.25
Q_Evap = -(depth*width)*M_Evap.*hfg_water; %Evaporative heat transfer in W, *-2.4 for ASU data, *-1.5 for TMY3 Data

end

        