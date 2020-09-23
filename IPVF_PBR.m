function [ IPVF ] = IPVF_PBR( space,height,width )
% this function calculates the inter panel view factor
% inputs: space - space between panels in m 
% height = height of panel in m 
% width = width of panel in m 

y = width/space; 
x = height/space; 
xd = sqrt(1 + x^2);
yd = sqrt(1 + y^2);
G1 = 2*log(((xd^2)*(yd^2))/(1+(x^2)+(y^2)))^.5;
G2 = 2*x*yd*atan(x/yd); 
G3 = 2*y*xd*atan(y/xd)-2*x*atan(x)-2*y*atan(y);
IPVF = (G1+G2+G3)/(pi*x*y);
end

