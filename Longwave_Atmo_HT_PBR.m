function [ atmo_q ] = Longwave_Atmo_HT_PBR( Tamb, height, space, absorp, width, transm )
% this function calculates distant infrared irradiance for the VFP
% Assumptions: no reflection from culture, therefore emissivity equals
% absorbance (Kirchoff's law) 

% Tamb = ambient temperature in K
% height = height of panel 
% space = space between panels in m
% absorp = spectral absorptivity (dimensionless) 
% transm = spectral transmissivity for material (dimensionless)
% width = width of panel

emiss_culture = absorp; % emissivity = absorption assumption
stb = 5.67e-8; % stefan boltz-mann constant
eps_atm = 0.8; % sky emissivity
F1 = (1+(height/space) - (1 + (height/space)^2)^0.5); % view factor above panels
Tsky = (eps_atm^0.25)*Tamb; % calculate sky temperature
Isky = eps_atm*stb*Tsky^4; 
atmo_q = F1*Isky*emiss_culture*transm*2*width*space; % atmospheric irradiance, in W

end

