function [re_rad] = Infrared_HT_PBR(height,width, space, TR,absorp,transm)
% Infrared out from panel and from oppposing panels from Endres et al. 
[ F2 ] = IPVF_PBR( space,height,width ); % get inter panel view factor

eps_g = absorp; % emissivity equals absorptivity
stb = 5.67 * 10^-8;

IR_out = -2*eps_g*(height*width)*stb*TR^4; % out from panel in W
IR_in = 2*eps_g*F2*transm*stb*(height*width)*TR^4; %in from opposing panel in W
re_rad = IR_out+IR_in;% reradiation in W 
end

