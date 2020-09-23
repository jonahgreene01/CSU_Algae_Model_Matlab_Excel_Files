function [Light_eff] = Light_Efficiency_PBR(IO,Conc_eff, photoI)
%This function determines the Light-Inhibition Metabolic Efficiency 

% Io in micro mol/m2/s
% photoI in micro mol/m2/s

Iav = IO*Conc_eff;
Light_eff = (Iav/(10^-6*photoI))*exp(1-(Iav/(photoI*10^-6)));
   
end

     