function [ q_bub ] = Sparge_Evaporation_HT_PBR( T, Tamb, RH, v_air )
% this function calculates heat transfer due to evaporation of water into
% the sparge gas

% Inputs: T = temperature in K
% Tamb = ambient temperature in K 
% RH = relative humidity (%) 
% v_air = velocity of air in m/s

ptot = 101.315*1000; % atmospheric pressure

[density, ~, hfg, ~, ~, ~, ~ ] =Air_Properties(T);

pg = (101.325*(10^3)/760)*exp(20.36-(5132/T)); % partial pressure of water vapor
mda = v_air*density; % mass flowrate of sparge gas
om1  = 0.622*(RH/100)*pg/(ptot - pg*(RH/100)); % mass fraction of water vapor in fluid-air interface
om2 = 0.622*pg/(ptot-pg);  % mass fraction of water vapor in the ambient air

q_bub =-1000*hfg*mda*(om2-om1); % watts of evaporative cooling

if Tamb > T
    q_bub = 0; 
end

end

