function [diffuse_ht] = Diffuse_Rad_HT_PBR(F1,DHI,space,width,absorp,transm)
%calculates incident diffuse irradiance in W/m2 and reradiation from panel 

stb = 5.67*10^-8;
IDF = F1*DHI; % incident diffuse irradiance in W/m2
diffuse_ht = 2*IDF*space*width*absorp*transm; % diffuse solar irradiance in W

end
