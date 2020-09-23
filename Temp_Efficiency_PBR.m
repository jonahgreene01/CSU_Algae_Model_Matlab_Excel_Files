function [Temp_eff] =Temp_Efficiency_PBR(TR, T_opt,Tmin,Tmax)
%This function determines the temperature efficiency factor using
% a cardinal model based on  rosso et al 1993

if TR < Tmin
    Temp_eff = 0; 
end
if TR >= Tmin && TR <= Tmax
    Temp_eff = ((TR-Tmax).*(TR-Tmin).^2)./((T_opt-Tmin).*((T_opt-Tmin).*(TR-T_opt)-(T_opt-Tmax).*(T_opt+Tmin-2*TR)));
    
end
if TR > Tmax
    Temp_eff = 0;
end

if TR == T_opt
    Temp_eff=1; 
end

end
