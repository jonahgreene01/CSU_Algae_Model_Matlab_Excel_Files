function [Q_Rerad] = Reradiation_HT_ORP(TR)
%This function defines the heat flux from reradiation from the pond (Q_out)

sigma=5.67*10^-8;
epsilon_water = 0.97;

Q_Rerad = -1*sigma*epsilon_water.*(TR.^4); 

% This function is the pond radiation. Yadala and Cremaschi, 2016. The
% culture is approximated to have the emissivity of water.
% Stephan-Boltzmann fourth power law estimation.
end

