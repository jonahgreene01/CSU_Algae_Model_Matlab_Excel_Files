function [ ground_temp ] = Ground_Temperature_PBR(GHI,Tgi,Tamb,space,width,dt,depth)
% this function calculates the ground temperature beneath the panels
% following the methods detailed in Enders et al. 2011 

% GHI = global horizontal irradiance in W/m2
% Tgi = initial ground temperature in K
% Tamb = ambient temperature in K
% depth = ground depth in m
% space = space between panels in m
% width = width of panel in m 
% dt = time step 
% height = height of panel in m

stb = 5.67e-8; % stefan boltz-mann constant
eps_g = .90; % ground emissivity
eps_atm =.8; % atmospheric emissivity def .8
absg = 0.82; % ground spectral absorptivity .82
g = 9.81; % gravitational acceleration 
rhog = 2400; % ground density in kg/m3
cpg = 1260; % typical for concrete in J/K/kg
k = .50; % thermal conductivity W/m*K typical for concrete

n = 52; % number of approximations for ground discretization 

b = 1/((Tgi + Tamb)/2); % assuming ideal gas behavior for air at the film temp
tg(n) = 290; % assumed distant ground temperature in K
tg(1) = Tgi; % initial temperature is set to ambient temperature for i =1 
m = (tg(n)-tg(1))/n; % step size of linear interpolation between layers
dx = 2.6/n; % depth of each layer in meters, =2.6 is depth in meters

 for j =2:(n-1) % 2:n-1 because surface temperature is fixed
    tg(j) = m*j+Tgi;
 end
 x = (Tamb+Tgi)/2;

[density, ~, ~, ~, dyn_visc, pr, ~ ] = Air_Properties(x);
  
Ral = pr*g*b*(Tgi-Tamb)*((2*space)^3)/((dyn_visc/density)^2); % rayleigh number


% Nusselt number correlation

if Ral > 0 % calculate nusselt number
    NuL = 0.540*Ral^.25;
else 
    NuL = 0.52*abs(Ral)^(1/5);
end 

hg = k*NuL/(2*space); % heat transfer convection coefficient W/m2/K

direct_ground = absg*GHI* ((2*space+depth)*width)/(2*space*width); % direct irradiance to ground in W/m2
ground_convec = -(tg(1)-Tamb)*hg; % ground convection in W/m2
ground_conduction = -k*(tg(1)-tg(2))/dx; % conduction from surface to bottom layer in W/m2
ir_ground = -eps_g*stb*tg(1)^4; % IR radiation out from ground W/m2
ir_in = eps_atm*absg*stb*Tamb^4; % incoming IR from atmosphere W/m2

ground_flux_surface = direct_ground + ground_convec + ground_conduction + ir_ground + ir_in; % W/m2
ground_temp = Tgi+dt*ground_flux_surface/(cpg*rhog*dx); % final ground temperature in Kelvin
   

end

