function [ ref_sum, ref_dir, diff_ref,ref_infrar,ref_grnd ] = Reflection_HT_PBR(transm, height, space,width,absorp,hgt_p,width_p,dni,DHI,Tamb,ground_ref,ir_grnd,tc)
% this function calculates the heat transfer due to interpanel interactions

% Inputs: 
% Tamb = temperature in K
% height = height of the panel in m
% space = space between panels in m
% width = width of panel in m 
% absorp = spectral absorptivity (dimensionless) 
% transm = spectral transmissivity for material (dimensionless)
% az = azenith angle in degrees
% el = elevation angle in degrees
% north angle in radians
% dni=direct normal irradiance in W/m2
% DHI = diffuse horizontal irradiance in W
% reflected from ground 
% tc is angle of incidence in radians

eps_atm = 0.8; % sky emissivity
stb=5.67E-8; % stefan boltz-mann constant
F1 = (1+(height/space) - (1 + (height/space)^2)^0.5)/2; % view factor above panel
rho = 1- transm; % calculated panel reflectivity
IPVF1 = IPVF_PBR(space,height,width);
Tsky = (eps_atm^0.25)*Tamb;
eps_r = absorp; % grey body assumption 

% reflected direct sunlight from opposing panel
 
if hgt_p >= (height) 
    ref_dir = 0; 
elseif hgt_p <= height/2
        ref_dir = (width_p*hgt_p)*dni*abs(cos(tc))*rho*IPVF1*absorp*transm;
elseif hgt_p < height && hgt_p > height/2
            ref_dir = rho*transm*(height-hgt_p)*width_p*dni*abs(cos(tc))*IPVF1;
end

% reflected diffuse atmospheric irradiance from opposing panel
diff_ref = 2*width*space*transm*F1*rho*absorp*IPVF1*DHI; % 

% infrared irradiance reflected by the atmosphere to panel and panel to
% panel
ref_infrar = 2*eps_r*(1-eps_r)*F1*IPVF1*(height*width)*stb*eps_atm*Tsky^4;

% infrared irradiance reflected from the ground to a panel and panel to
% panel
ref_grnd  = 2*ground_ref*rho*IPVF1+ 2*ir_grnd*height*width*rho*IPVF1 ;


% reflected energy sum
ref_sum = (ref_grnd+ref_infrar+diff_ref+ref_dir); 

end






