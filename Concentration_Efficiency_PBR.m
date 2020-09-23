function [Conc_eff] = Concentration_Efficiency_PBR(CX_pbr, ODC_pbr, depth_pbr)
%This function determines the impact of the time-resolved algal concentration on the average light intensity hitting the culture

gpl = CX_pbr/1000; 
OD = gpl/ODC_pbr;
Conc_eff = (1-exp(-OD*depth_pbr))/(OD*depth_pbr);

end

