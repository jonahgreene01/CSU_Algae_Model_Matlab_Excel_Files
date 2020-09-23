
function [decay] = Night_Respiration_PBR(solar_flux, CX_pbr, night_resp, volume_pbr, Temp_eff_pbr)


if (solar_flux < 5)
    
    decay_rate = (log(1 - night_resp))/10; %yields the decay rate 1/day as a negative decimal
    decay_specific = -(CX_pbr*(decay_rate/3600));
    decay = (decay_specific*volume_pbr)*Temp_eff_pbr;

else 
    
    decay = 0; 
    
end 



 
