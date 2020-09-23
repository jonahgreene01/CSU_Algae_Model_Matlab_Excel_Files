function [ rad_dir,hgt_p ] = Direct_Rad_HT_PBR(transm,absorp,dni,space,az,el,north_angle,height,width)

% this function calculates the direct normal irradiance 

% absorp = spectral absorptivity for culutre (dimensionless) 
% transm = spectral transmissivity for material (dimensionless)
% dni=direct normal irradiance in W/m2
% space = space between panels in m
% az = azimuth angle in radians
% el = elevation angle in radoians
% north angle in radians
% height = height of the panel in m
% width = width of pabel in m 

v1 = [ sin(north_angle) cos(north_angle) 0]; % unit vector of panel orthogonal vector projection onto horizontal xy plane
v2 = [ sin(az)*cos(el) cos(az)*sin(el) sin(el)]; % unit vector to the sun
tc = acos(dot(v1,v2)); % angle between a vector orthogonal to the vertical face of PBR and the sun
v3 = [ v2(1) v2(2) 0]; % unit vector of sun projection in xy plane
v4 = [ 0 v2(2) v2(3)]; % unit vector of sun projection in yz plane
txy = acos(dot(v1,v3)); % angle between panel normal vector projected to xy plane
tyz = acos(dot(v1,v4)); % angle between panel normal vector and sun projected to yz plane
INN = abs(dni*cos(tc)); % incident normal irradiance in W/m2
width_p = abs(space*tan(txy)); % solar area due to shading
hgt_p = abs(space*tan(tyz)); % solar area due to shading


if hgt_p < 0.0 || hgt_p > height % solar area check
    hgt_p = height; 
end

if width_p < 0 || width_p > width
    width_p = width; 
end
area_p = width_p*hgt_p; % unshaded area

rad_dir =INN*area_p*absorp*transm; % direct radiation in W


end

