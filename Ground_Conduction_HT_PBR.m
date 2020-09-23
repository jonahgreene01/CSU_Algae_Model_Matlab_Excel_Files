function [ ground_flux, ground_reflected,ir_grnd ] = Ground_Conduction_HT_PBR( ground_temp, dni, DHI, hgt_p, width,space, height, absorp, depth,transm)
% this function calculates the heat transfer from the ground to the panel in W
%   Inputs: 
% ground_temp = temperature of the ground in K
% dni=direct normal irradiance in W/m2
% DHI = diffuse horizontal irradiance in W
% hgt_p = solar area due to shading
% height = height of the panel in m
% space = space between panels in m
% width = width of panel in m 
% absorp = spectral absorptivity of culture (dimensionless) 
% transm = spectral transmissivity for material (dimensionless)

alb = .18; % albedo, proportion of incident light reflected by surface 
em_g = .85;  % ground emissivity 
stb = 5.67e-8; % stefan boltz-mann constant
area = width*height;
ag = (5*space*width)-(5*depth*width); % ground area for 6 panels

if (ge(hgt_p,height)) && (le(hgt_p,5*height)) 
    AGDN = space*width*(hgt_p-height)/height;
else
    AGDN = 0; 
end

F1 = (1+(height/space) - (1 + (height/space)^2)^0.5)/2; % view factor above panels
ground_reflected = 2*((ag*DHI*alb*F1/area) + (AGDN*F1*dni*alb/area))*area*transm*absorp; % direct and diffuse reflected from ground in W
ir_grnd = 2*stb*em_g*(ground_temp^4)*F1*absorp*transm; % infrared raditation from ground W

ground_flux = ground_reflected + ir_grnd*area; % heat transfer to panels from the ground 


end

