function [Q_Longwave_Atmo] = Longwave_Atmo_HT_ORP(T_amb)
%This function defines the heat flux from longwave atmospheric radiation
epsilon_water = 0.97; %emissivity of water under normal conditions
epsilon_air = 0.85; %emissivity of air
sigma = 5.67*(10^-8); %Stefan Boltzmann constant (W/m^2*k^4)

Q_Longwave_Atmo = epsilon_water*epsilon_air*sigma.*(T_amb.^4);

% This comes directly from Yadala and Cremaschi, 2016, however T_amb is
% substituted for T_surr in the original equation. It is an approximation
% to use T_amb


end


