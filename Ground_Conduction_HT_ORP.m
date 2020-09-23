function [Q_Ground] = Ground_Conduction_HT_ORP (TR)
%This function determines the heat flux between the ground and the ORP
diff_concrete=691.70*(10^-9); 
k_concrete = 1.4;
l_ref=4400*diff_concrete^0.5;
Q_Ground=-k_concrete*(TR-290.0)/l_ref;

% This heat flux equation comes from Bechet et al., "Universal temperature
% model for shallow algal ponds provides improved accuracy" and the
% equations are found in the SI document. It is approximating a depth in
% meters (3.65m) at which soil temp is unaffected by changes in ambient
% environment. Then it calculates a simple 1-D conduction through concrete
% using the thermal conductivty of concrete, and the the two temperatures
% (pond and soil). Soil temp at reference depth is assumed to equal 290K.
end
