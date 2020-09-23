function [ convection_Q ] = Convection_HT_PBR( T, Tamb, wdir, wspd, north_angle, height,width )
% this function calculates the heat transfer via convection

% Inputs: T  = temperature in kelvins
% Tamb = ambient temperature in K 
% wdir = wind direction in degrees
% wspd = wind speed in m/s 
% north angle in radians
% heigth of panel in m
% width of panel in m 

g = 9.81; % gravitational acceleration in m/s2

% correlation to adjust for height from prevailing windspeed
z = 0.95 + height/2 ; % center height measured from ground

if wspd == 0 
    na = 0; 
else 
    na = (.37-.0881*log10(abs(wspd)))/(1-.0881*log10(z/10));
end

cf_z = (z/10)^na; 
VR = [ sin(north_angle) cos(north_angle) ]; 
VW = [ sin(wdir*3.14159/180) cos(wdir*3.14159/180)]; 
WC = acos(dot(VR,VW)); 
VC = abs(cf_z*wspd*sin(WC));

% calculating the air properties at the film temperature 

[ density, conductivity, ~, ~, dyn_visc, pr, ~ ] = Air_Properties( (T+Tamb)/2 );

ReL  = VC*density*width/dyn_visc; % reynolds number 

beta = 1 / (0.5*(Tamb+T)); % assuming ideal gas behavior for air 
kin_visc = dyn_visc/density; % kinematic viscosity 
GrL = (g*beta*abs(T-Tamb)*(height^3))/(kin_visc^2); % grasshof number 
RaL = GrL*pr ; % Rayleigh  number

if GrL/(ReL^2) < 1
    NuL = 0.68*(ReL^.5)*(pr^(1/3)); % nuselt number correlation for forced convection heat transfer
else
    top = .387*RaL^(1/6);
    bottom = (1+(.492/pr)^(9/16))^(8/27);
    NuL = ((0.825 + top)/bottom)^2;% nusselt number correlation for vertical plate convecting naturally
end

hc = NuL*conductivity/height; % convection coefficient W/m2 K 
convection_Q = hc*(Tamb - T)*width*height*2 ; % Watts of convective heat transfer


end

