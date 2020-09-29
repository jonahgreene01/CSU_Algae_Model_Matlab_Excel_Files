%Master PBR Growth + Thermal Model for Open Source Model
clc
clear

%Import data for T_amb, GHI, RH, and WNDSPD from TMY3 Dataset
TMY3_data = xlsread('C:\CSU_Open_Source_Algae_Model\Model_Files\CSU_Algae_Model_User_Interface.xlsm', 3, 'E3:AU17522');
Tamb= TMY3_data(:,28)+273.15;
GHI= TMY3_data(:,1);
DNI= TMY3_data(:,4);
DHI= TMY3_data(:,7);
RH= TMY3_data(:,34);
WNDSPD= TMY3_data(:, 43);
WDIR = TMY3_data(:, 40);

file = 'C:\CSU_Open_Source_Algae_Model\Model_Files\CSU_Algae_Model_User_Interface.xlsm';
sheet = 3;
[W_data,date,~] = xlsread(file,sheet, 'A3:B17522');
time = datestr(W_data(:,1),'HH:MM:SS');
date_str= char(date(:,1));
UTC = strcat(date_str,{' '},time);



%Read Inputs from Excel
Inputs_Data = xlsread('C:\CSU_Open_Source_Algae_Model\Model_Files\CSU_Algae_Model_User_Interface.xlsm', 1, 'C2:B31');

%Assign values to variables within code

%PBR Geometry Paramterers
height_pbr = Inputs_Data(9,:);
width_pbr = Inputs_Data(10,:);
depth_pbr = Inputs_Data(11,:);
space_pbr = Inputs_Data(12, :);
v_air = Inputs_Data(13, :); %sparge velocity in m/s

%Strain Parameters
T_opt_PBR = Inputs_Data(14,:);
T_min = Inputs_Data(15,:); 
T_max = Inputs_Data(16,:);
ODC_pbr = Inputs_Data(17,:);
I_sat = Inputs_Data(18,:)*4.56;
night_resp= Inputs_Data(19,:);

%Harvesting Parameters
IBC_pbr = Inputs_Data(20,:);
harvest_conc_pbr = Inputs_Data(21,:);
harvest_days_pbr = Inputs_Data(22,:);

%Time Positioning
n1 = Inputs_Data(23,:);
n2 = Inputs_Data(24,:);  

%Algal Composition
Algae_N_Content = Inputs_Data(25,:); 
Algae_P_Content = Inputs_Data(26,:);
Nut_surplus = Inputs_Data(27,:); %nutrient surplus (% of requirement)

%CO2 Parameters
Mass_frac_CO2 = Inputs_Data(28,:); 
CO2_util = Inputs_Data(29,:); 
CO2_Source_Percent_CO2 = Inputs_Data(30,:);

%Lat, Long, Altitude, and north angle for solar angles code
Lat = Inputs_Data(1, :); 
Lon = Inputs_Data(2, :); 
Alt = Inputs_Data(3, :);  % altitude in km
north_angle = Inputs_Data(4, :); %north angle in radians

%Define unchanging parameters
dt = 3600; %time step of 1 hour or 3600 s
absorp =0.7; % absorptivity
transm = 0.92; % transmissivity of the panels 
rho_algae = 1000.0;  %kg/m3
cp_algae = 4184.0; 
phi = 12/9;
Conversion = 4.56 * 10^-6 ; %W/m2 to mol/m2*s
area_pbr = width_pbr*height_pbr; 
volume_pbr = width_pbr*height_pbr*depth_pbr; 
perimeter_pbr = 2*(width_pbr+height_pbr); 

%Nutrient Parameters
Ammonia_N_Content = 0.82; 
DAP_N_Content = 0.18;
DAP_P_Content = 0.20; 

%Initialize Variables
Harvest_Mass_pbr = 0.0; 
marker_pbr = 0;
pbr_harvests = 0;

%Prealocate all variables that change in the for loop to save time
% set vector size depending on number of data points
TR_pbr = zeros(n2, 1);
dT1_pbr = zeros(n2, 1);
T1_pbr = zeros(n2, 1);
dT2_pbr = zeros(n2, 1);
T2_pbr = zeros(n2, 1);
dT3_pbr = zeros(n2, 1);
T3_pbr = zeros(n2, 1);
dTdt_pbr = zeros(n2, 1);
Tg_pbr = zeros(n2,1);
Temp_eff = zeros(n2,1);
Conc_eff = zeros(n2,1);
Light_eff = zeros(n2,1);
CX_pbr = zeros(n2,1); 
Conc_at_Harvest_pbr = zeros(n2,1);
I_ave_culture = zeros(n2,1);

hgt_p = zeros(n2,1); 
IDF = zeros(n2,1); 
diffuse_ht = zeros(n2,1); 
direct_rad = zeros(n2,1);
infrared_rad = zeros(n2,1); 
q_bubble = zeros(n2,1); 
evaporation_ht = zeros(n2,1); 
convection_ht = zeros(n2,1); 
ground_q = zeros(n2,1); 
atmospheric_ht = zeros(n2,1); 
solar_flux = zeros(n2,1);
direct_top = zeros(n2,1);
decay = zeros(n2,1); 
dCXdt = zeros(n2,1); 

%Prealocate Consumption Variables
Water_consump = zeros(n2,1);
CO2_demand_hrly = zeros(n2,1);
Ammonia_demand_hrly = zeros(n2,1);
DAP_demand_hrly = zeros(n2,1);

%Compute thermal mass based on algal properties
Therm_mass = volume_pbr*rho_algae*cp_algae;
F1 = (1 + (height_pbr/space_pbr)-(1+(height_pbr/space_pbr)^2)^.5)/2; % Compute Interpanel view factor
[az, el] = SolarAzEl_PBR(UTC,Lat,Lon,Alt); % get azimuth, and elevation angles
az = az*3.14159/180.0; % convert azimuth angle to radians
el = el*3.14159/180.0; % convert elevation angle to radians

%% This is the start of the for loop 
for i = n1:(n2-1) %Only go to n-1 because program reads to i+1
 
if (i == n1)
    TR_pbr(i,1) =	Tamb(i,1);
    CX_pbr(i,1) = IBC_pbr;

else
    TR_pbr(i,1) = TR_pbr(i,1);
    CX_pbr(i,1) = CX_pbr(i,1);

end

%Compute Air Properties at this time step using the film temperature
[rho_air, ~, ~, ~, dyn_visc_air, ~, ~]= Air_Properties((TR_pbr(i,1)+Tamb(i,1))/2);

%Compute ground temperature
Tg_pbr(1,1) = Tamb(1,1);

%Calculate incidence angle and shading parameters
v1 = [ sin(north_angle) cos(north_angle) 0]; % vector of panel orthogonal vector projection onto horizontal xy plane
v2 = [ sin(az(i,1))*cos(el(i,1)) cos(az(i,1))*sin(el(i,1)) sin(el(i,1))]; % vector to the sun
tc = acos(dot(v1,v2)); % angle between a vector orthogonal to the vertical face of PBR and the sun
v3 = [ v2(1) v2(2) 0]; 
v4 = [ 0 v2(2) v2(3)];
txy = acos(dot(v1,v3)); % angle between panel normal vector projected to xy plane
tyz = acos(dot(v1,v4)); % angle between panel normal vector and sun projected to yz plane
INN = abs(DNI(i,1)*cos(tc)); % incident normal irradiance in W/m2
width_p = abs(space_pbr*tan(txy)); % solar area due to shading
hgt_p = abs(space_pbr*tan(tyz)); % solar area due to shading

if (hgt_p < 0) || (hgt_p > height_pbr) % solar area check
    hgt_p = height_pbr;
end

if width_p < 0 || width_p > width_pbr
    width_p = width_pbr;
end


%First RK4 Fluxes
[ direct_rad,~ ] = Direct_Rad_HT_PBR(transm,absorp,DNI(i,1),space_pbr,az(i,1),el(i,1),north_angle,height_pbr,width_pbr);
[ ground_q,ref_grnd, ir_grnd ] = Ground_Conduction_HT_PBR( Tg_pbr(i,1), DNI(i,1), DHI(i,1), hgt_p, width_pbr, space_pbr, height_pbr, absorp,depth_pbr,transm);
[ atmospheric_ht ] = Longwave_Atmo_HT_PBR( Tamb(i,1),height_pbr,space_pbr, absorp,width_pbr,transm);
[ ref_sum,ref_direct ,diffuse_ref,~,~ ] = Reflection_HT_PBR(transm, height_pbr, space_pbr,width_pbr,absorp,hgt_p,width_p,DNI(i,1),DHI(i,1), Tamb(i,1),ref_grnd,ir_grnd,tc);
[diffuse_ht] = Diffuse_Rad_HT_PBR(F1,DHI(i,1),space_pbr,width_pbr,absorp,transm); 
[infrared_rad] = Infrared_HT_PBR(height_pbr,width_pbr, space_pbr, TR_pbr(i,1),absorp,transm);
[M_Evap, evaporation_ht] = Evaporative_HT_PBR(depth_pbr, width_pbr, TR_pbr(i,1), Tamb(i,1), RH(i,1), WNDSPD(i,1), dyn_visc_air, rho_air);
[ convection_ht ] = Convection_HT_PBR( TR_pbr(i,1), Tamb(i,1), WDIR(i,1),WNDSPD(i,1),north_angle, height_pbr,width_pbr );
[ q_bubble ] = Sparge_Evaporation_HT_PBR( TR_pbr(i,1), Tamb(i,1), RH(i,1), v_air );
HT = infrared_rad+diffuse_ht+ref_sum+atmospheric_ht+ground_q+direct_rad+convection_ht+evaporation_ht+q_bubble + GHI(i,1)*width_pbr*depth_pbr*absorp*transm;

%RK4_1 time derivative of temperature - K1=dT1
dT1_pbr = HT/Therm_mass;
T1_pbr = TR_pbr(i,1) + 0.5*dt*dT1_pbr;

solar_flux(i,1) =   direct_rad*2+diffuse_ht+ref_direct*2+diffuse_ref;% radiation contributing to algal growth in W/m2

%Biological Growth Model
[Temp_eff(i,1)] = Temp_Efficiency_PBR(TR_pbr(i,1), T_opt_PBR,T_min,T_max);
[Conc_eff(i,1)] = Concentration_Efficiency_PBR(CX_pbr(i,1), ODC_pbr, depth_pbr);
[Light_eff(i,1)] = Light_Efficiency_PBR(0.458*solar_flux(i,1)*Conversion,Conc_eff(i,1), I_sat);
[decay(i,1)] = Night_Respiration_PBR(solar_flux(i,1), CX_pbr(i,1), night_resp, volume_pbr, Temp_eff(i,1)); % biomass loss due to dark respiration
dCXdt(i,1) = Light_eff(i,1)*Temp_eff(i,1)*phi*0.458*solar_flux(i,1)*Conversion/volume_pbr; % autotrophic growth rate in g/m3/s
CX_pbr(i+1,1) = CX_pbr(i,1)+(dCXdt(i,1)-(decay(i,1)/volume_pbr))*dt; % algae concentration g/m3


%Second RK4 Fluxes
[ direct_rad,~ ] = Direct_Rad_HT_PBR( transm,absorp,((DNI(i,1)+DNI(i+1,1))/2),space_pbr,az(i,1),el(i,1),north_angle,height_pbr,width_pbr);
[ ground_q,ref_grnd,ir_grnd ] = Ground_Conduction_HT_PBR( Tg_pbr(i,1), ((DNI(i,1)+DNI(i+1,1))/2), ((DHI(i,1)+DHI(i+1,1))/2), hgt_p, width_pbr, space_pbr, height_pbr, absorp,depth_pbr,transm);
[ atmospheric_ht ] = Longwave_Atmo_HT_PBR(((Tamb(i,1)+Tamb(i+1,1))/2),height_pbr,space_pbr, absorp, width_pbr,transm );
[  ref_sum,~ ,~,~,~  ] = Reflection_HT_PBR(transm, height_pbr, space_pbr,width_pbr,absorp,hgt_p,width_p,((DNI(i,1)+DNI(i+1,1))/2),((DHI(i,1)+DHI(i+1,1))/2),((Tamb(i,1)+Tamb(i+1,1))/2),ref_grnd,ir_grnd,tc);
[diffuse_ht] = Diffuse_Rad_HT_PBR(F1,(DHI(i,1)+DHI(i+1,1))/2,space_pbr,width_pbr,absorp,transm);   
[infrared_rad] = Infrared_HT_PBR(height_pbr,width_pbr, space_pbr, T1_pbr,absorp,transm);
[M_Evap, evaporation_ht] = Evaporative_HT_PBR(depth_pbr, width_pbr, T1_pbr, ((Tamb(i,1)+Tamb(i+1,1))/2), ((RH(i,1)+RH(i+1,1))/2), ((WNDSPD(i,1)+WNDSPD(i+1,1))/2), dyn_visc_air, rho_air);
[ convection_ht ] = Convection_HT_PBR( T1_pbr, ((Tamb(i,1)+Tamb(i+1,1))/2), ((WDIR(i,1)+WDIR(i+1,1))/2),((WNDSPD(i,1)+WNDSPD(i+1,1))/2),north_angle, height_pbr,width_pbr );
[ q_bubble ] = Sparge_Evaporation_HT_PBR( T1_pbr, ((Tamb(i,1)+Tamb(i+1,1))/2), ((RH(i,1)+RH(i+1,1))/2), v_air );
HT = infrared_rad+diffuse_ht+ref_sum+atmospheric_ht+ground_q+direct_rad+convection_ht+evaporation_ht+q_bubble + ((GHI(i,1)+GHI(i+1,1))/2)*width_pbr*depth_pbr*absorp*transm;

%RK4_2 time derivative of temperature - K2=dT2
dT2_pbr = HT/Therm_mass;
T2_pbr = TR_pbr(i,1) + 0.5*dt*dT2_pbr;

%Third RK4 Fluxes
[ direct_rad,~ ] = Direct_Rad_HT_PBR( transm, absorp,((DNI(i,1)+DNI(i+1,1))/2),space_pbr,az(i,1),el(i,1),north_angle,height_pbr,width_pbr);
[ ground_q,ref_grnd,ir_grnd ] = Ground_Conduction_HT_PBR(Tg_pbr(i,1), ((DNI(i,1)+DNI(i+1,1))/2), ((DHI(i,1)+DHI(i+1,1))/2), hgt_p, width_pbr, space_pbr, height_pbr, absorp,depth_pbr,transm);
[ atmospheric_ht ] = Longwave_Atmo_HT_PBR(((Tamb(i,1)+Tamb(i+1,1))/2),height_pbr,space_pbr, absorp,width_pbr,transm);
[  ref_sum,~ ,~,~,~] = Reflection_HT_PBR(transm, height_pbr, space_pbr,width_pbr,absorp,hgt_p,width_p,((DNI(i,1)+DNI(i+1,1))/2),((DHI(i,1)+DHI(i+1,1))/2),((Tamb(i,1)+Tamb(i+1,1))/2),ref_grnd,ir_grnd,tc);
[diffuse_ht] = Diffuse_Rad_HT_PBR(F1,(DHI(i,1)+DHI(i+1,1))/2,space_pbr,width_pbr,absorp,transm);
[infrared_rad] = Infrared_HT_PBR(height_pbr,width_pbr, space_pbr, T2_pbr,absorp,transm);
[M_Evap, evaporation_ht] = Evaporative_HT_PBR(depth_pbr, width_pbr, T2_pbr, ((Tamb(i,1)+Tamb(i+1,1))/2), ((RH(i,1)+RH(i+1,1))/2), ((WNDSPD(i,1)+WNDSPD(i+1,1))/2), dyn_visc_air, rho_air);
[ convection_ht ] = Convection_HT_PBR( T2_pbr, ((Tamb(i,1)+Tamb(i+1,1))/2), ((WDIR(i,1)+WDIR(i+1,1))/2),((WNDSPD(i,1)+WNDSPD(i+1,1))/2),north_angle, height_pbr,width_pbr );
[ q_bubble ] = Sparge_Evaporation_HT_PBR( T2_pbr, ((Tamb(i,1)+Tamb(i+1,1))/2), ((RH(i,1)+RH(i+1,1))/2), v_air );
HT = infrared_rad+diffuse_ht+ref_sum+atmospheric_ht+ground_q+direct_rad+convection_ht+evaporation_ht+q_bubble + ((GHI(i,1)+GHI(i+1,1))/2)*width_pbr*depth_pbr*absorp*transm;

%RK4_3 time derivative of temperature - K3=dT3
dT3_pbr = HT/Therm_mass;
T3_pbr = TR_pbr(i,1) + 0.5*dt*dT3_pbr;

%Fourth RK4 Fluxes
[ direct_rad,~] = Direct_Rad_HT_PBR(transm, absorp,DNI(i+1,1),space_pbr,az(i,1),el(i,1),north_angle,height_pbr,width_pbr);
[ ground_q,ref_grnd,ir_grnd ] = Ground_Conduction_HT_PBR(Tg_pbr(i,1), DNI(i+1,1), DHI(i+1,1), hgt_p, width_pbr, space_pbr, height_pbr, absorp,depth_pbr,transm);
[ atmospheric_ht ] = Longwave_Atmo_HT_PBR(Tamb(i+1,1),height_pbr,space_pbr, absorp,width_pbr,transm);
[  ref_sum,~ ,~,~,~] = Reflection_HT_PBR(transm, height_pbr, space_pbr,width_pbr,absorp,hgt_p,width_p,DNI(i+1,1),DHI(i+1,1), Tamb(i+1,1),ref_grnd,ir_grnd,tc);
[diffuse_ht] = Diffuse_Rad_HT_PBR(F1,DHI(i+1,1),space_pbr,width_pbr,absorp,transm);
[infrared_rad] = Infrared_HT_PBR(height_pbr,width_pbr, space_pbr, T3_pbr,absorp,transm);
[M_Evap, evaporation_ht] = Evaporative_HT_PBR(depth_pbr, width_pbr, T3_pbr, Tamb(i+1,1), RH(i+1,1), WNDSPD(i+1,1), dyn_visc_air, rho_air);
[ convection_ht ] = Convection_HT_PBR( T3_pbr, Tamb(i+1,1),WDIR(i+1,1),WNDSPD(i+1,1),north_angle, height_pbr,width_pbr );
[ q_bubble ] = Sparge_Evaporation_HT_PBR( T3_pbr,Tamb(i+1,1),RH(i+1,1),v_air );
HT = infrared_rad+diffuse_ht+ref_sum+atmospheric_ht+ground_q+direct_rad+convection_ht+evaporation_ht+q_bubble +(GHI(i+1,1))*width_pbr*depth_pbr*absorp*transm;

%RK4_4 time derivative of temperature - K4=dT4
dTdt_pbr = HT/Therm_mass;
TR_pbr(i+1,1)= TR_pbr(i,1) + (1/6)*dt*(dT1_pbr+2.0*(dT2_pbr+dT3_pbr+dTdt_pbr));

%Track Average light intensity with each iteration
I_ave_culture(i,1) = solar_flux(i,1)*0.458*Conc_eff(i,1); 

marker_pbr = marker_pbr + 1; 
%Harvesting sequence
    if (((mod(marker_pbr, (harvest_days_pbr*24.0))== 0) && marker_pbr ~= 0)||(CX_pbr(i+1,1) >= harvest_conc_pbr)||(i==(n2-1)))
        Harvest_Mass_pbr = Harvest_Mass_pbr + (CX_pbr(i+1,1)-IBC_pbr)*volume_pbr; % harvested mass in grams
        Conc_at_Harvest_pbr(i,1) = CX_pbr(i+1,1);% g/m3
        CX_pbr(i+1,1) = IBC_pbr;
        marker_pbr = 0;
        pbr_harvests = pbr_harvests + 1;
    end

%Update ground temperature
[ Tg_pbr(i+1,1) ] = Ground_Temperature_PBR(GHI(i,1),Tg_pbr(i,1),Tamb(i,1),space_pbr,width_pbr,dt,depth_pbr);
    
%Calculate waterloss at each time step 
Water_consump(i,1) = M_Evap*depth_pbr*width_pbr*dt; %[kg/m2*s]*[m2]*[3600 s] = [kg/hr]

    if Water_consump(i,1) < 0
       Water_consump(i,1) = 0;
    else 
       Water_consump(i,1) = Water_consump(i,1); 
    end    
        
        
%Calculate the CO2 demand kg/hr based on stoichiometric carbon balance 
    if (i == n1)
        CO2_demand_hrly(i,1) = 0;
    else
        CO2_demand_hrly(i,1) = (((CX_pbr(i+1,1)-CX_pbr(i,1))*volume_pbr/1000)*Mass_frac_CO2)/CO2_util;
    end

    if  CO2_demand_hrly(i,1) < 0
        CO2_demand_hrly(i,1) = 0;
    else   
        CO2_demand_hrly(i,1) = CO2_demand_hrly(i,1);
    end
       
    
%Calculate the Nutrient consumption at each time step
    if (i == n1)
        Ammonia_demand_hrly(i,1) = 0;
        DAP_demand_hrly(i,1) = 0;

    else
        DAP_demand_hrly(i,1) = (((CX_pbr(i+1,1)-CX_pbr(i,1))*volume_pbr/1000)*Algae_P_Content/DAP_P_Content); %g/m3 * m3 * 1 kg/1000 g * % 
        Ammonia_demand_hrly(i,1) = (((CX_pbr(i+1,1)-CX_pbr(i,1))*volume_pbr/1000)*Algae_N_Content/Ammonia_N_Content)-(DAP_demand_hrly(i,1)*DAP_N_Content);
    end

    if  Ammonia_demand_hrly(i,1) < 0
        Ammonia_demand_hrly(i,1) = 0;
    else   
        Ammonia_demand_hrly(i,1) = Ammonia_demand_hrly(i,1)*(1+ Nut_surplus);
    end
    
    if  DAP_demand_hrly(i,1) < 0
        DAP_demand_hrly(i,1) = 0;
    else   
        DAP_demand_hrly(i,1) = DAP_demand_hrly(i,1)*(1+ Nut_surplus);
    end
        
end 

%% For Loop to store thermal flux values at each iteration

for i = n1:(n2)
% Information to print heat fluxes

        % Calculate incidence angle and shading parameters
        v1 = [ sin(north_angle) cos(north_angle) 0]; % vector of panel orthogonal vector projection onto horizontal xy plane
        v2 = [ sin(az(i,1))*cos(el(i,1)) cos(az(i,1))*sin(el(i,1)) sin(el(i,1))]; % vector to the sun
        tc = acos(dot(v1,v2)); % angle between a vector orthogonal to the vertical face of PBR and the sun
        v3 = [ v2(1) v2(2) 0]; 
        v4 = [ 0 v2(2) v2(3)];
        txy = acos(dot(v1,v3)); % angle between panel normal vector projected to xy plane
        tyz = acos(dot(v1,v4)); % angle between panel normal vector and sun projected to yz plane
        INN = abs(DNI(i,1)*cos(tc)); % incident normal irradiance in W/m2
        width_p = abs(space_pbr*tan(txy)); % solar area due to shading
        hgt_p = abs(space_pbr*tan(tyz)); % solar area due to shading
        
        if (hgt_p < 0) || (hgt_p > height_pbr) % solar area check
            hgt_p = height_pbr;
        end
        
        if width_p < 0 || width_p > width_pbr
            width_p = width_pbr;
        end

        [ direct_rad(i,1),~ ] = Direct_Rad_HT_PBR(transm,absorp,DNI(i,1),space_pbr,az(i,1),el(i,1),north_angle,height_pbr,width_pbr);
        [ ground_q(i,1),ref_grnd, ir_grnd ] = Ground_Conduction_HT_PBR( Tg_pbr(i,1), DNI(i,1), DHI(i,1), hgt_p, width_pbr, space_pbr, height_pbr, absorp,depth_pbr,transm);
        [ atmospheric_ht(i,1) ] = Longwave_Atmo_HT_PBR( Tamb(i,1),height_pbr,space_pbr, absorp,width_pbr,transm);
        [ ref_sum(i,1),~ ,~,~,~ ] = Reflection_HT_PBR(transm, height_pbr, space_pbr,width_pbr,absorp,hgt_p,width_p,DNI(i,1),DHI(i,1), Tamb(i,1),ref_grnd,ir_grnd,tc);
        [diffuse_ht(i,1)] = Diffuse_Rad_HT_PBR(F1,DHI(i,1),space_pbr,width_pbr,absorp,transm); 
        [infrared_rad(i,1)] = Infrared_HT_PBR(height_pbr,width_pbr, space_pbr, TR_pbr(i,1),absorp,transm);
        [~, evaporation_ht(i,1)] = Evaporative_HT_PBR(depth_pbr, width_pbr, TR_pbr(i,1), Tamb(i,1), RH(i,1), WNDSPD(i,1), dyn_visc_air, rho_air);
        [ convection_ht(i,1) ] = Convection_HT_PBR( TR_pbr(i,1), Tamb(i,1), WDIR(i,1),WNDSPD(i,1),north_angle, height_pbr,width_pbr );
        [ q_bubble(i,1) ] = Sparge_Evaporation_HT_PBR( TR_pbr(i,1), Tamb(i,1), RH(i,1), v_air );
        direct_top(i,1)= GHI(i,1)*width_pbr*depth_pbr*absorp*transm;
    
end
%% Results and formatting

%Summary of CO2 Consumption
CO2_consump_total = sum(CO2_demand_hrly(n1:n2));
CO2_consump_daily = CO2_consump_total/((n2-n1)/24); 

%Temporal Outputs
T_amb_Plot = Tamb(n1:n2, 1);
TR_Plot=TR_pbr(n1:n2,1); 
CX_plot = CX_pbr(n1:n2, 1); 
CO2_hourly_kg = CO2_demand_hrly(n1:n2,1); 
Water_con = Water_consump(n1:n2,1);
Ammonia_con = Ammonia_demand_hrly(n1:n2,1); 
DAP_con = DAP_demand_hrly(n1:n2,1); 
I_ave_PBR = I_ave_culture(n1:n2,1); 

%Plot Temporally Resolved Outputs into Excel
Temporal_Outputs = [T_amb_Plot, TR_Plot, CX_plot, CO2_hourly_kg, Water_con, Ammonia_con, DAP_con, I_ave_PBR];
xlswrite('C:\CSU_Open_Source_Algae_Model\Model_Files\Growth_Model_Outputs.xlsx', Temporal_Outputs, 2);

%Single Number Outputs
Areal_Prod = Harvest_Mass_pbr/((depth_pbr+space_pbr)*width_pbr*((n2-n1)/(24.0)));
Vol_Prod = Harvest_Mass_pbr/(volume_pbr*1000*((n2-n1)/(24.0)));
CO2_consump_ave = CO2_consump_daily; 
Water_consump_ave = sum(Water_con)/((n2-n1)/24);
Ammonia_consump_ave = sum(Ammonia_con)/((n2-n1)/24);
DAP_consump_ave = sum(DAP_con)/((n2-n1)/24); 
Conc_at_Harvest_Plot = Conc_at_Harvest_pbr(n1:n2,1);
Conc_at_Harvest_Ave = mean(nonzeros(Conc_at_Harvest_Plot)); 
Evap_rate = Water_consump_ave*(1/1000)*(1/(width_pbr*depth_pbr))*100; 
Op_Days_per_Year = (n2-n1)/24;

%Write Single Number Outputs to Excel
Growth_Model_Static_Outputs = [Areal_Prod; Vol_Prod; Water_consump_ave; CO2_consump_ave; Ammonia_consump_ave; DAP_consump_ave; Op_Days_per_Year; Conc_at_Harvest_Ave; Evap_rate]; 
xlswrite('C:\CSU_Open_Source_Algae_Model\Model_Files\Growth_Model_Outputs.xlsx', Growth_Model_Static_Outputs, 'B3:B11')

%Plot Temporally Resolved Thermal Fluxes
Q_radiation = direct_rad + direct_top + diffuse_ht; % total combined diffuse and direct solar irradiance
Q_evaporation = evaporation_ht + q_bubble; %evaporation total

Q_Solar_Plot = Q_radiation(n1:n2,1); 
Q_Atmospheric_Plot = atmospheric_ht(n1:n2,1);
Q_Infrared_Plot = infrared_rad(n1:n2,1); 
Q_Ground_Plot = ground_q(n1:n2,1); 
Q_Reflection_Plot = ref_sum(n1:n2,1);
Q_Convection_Plot = convection_ht(n1:n2,1); 
Q_Evaporation_Plot = Q_evaporation(n1:n2,1); 


Thermal_Flux_Outputs = [Q_Solar_Plot, Q_Atmospheric_Plot, Q_Infrared_Plot, Q_Ground_Plot, Q_Reflection_Plot, Q_Convection_Plot, Q_Evaporation_Plot];
xlswrite('C:\CSU_Open_Source_Algae_Model\Model_Files\Growth_Model_Outputs.xlsx', Thermal_Flux_Outputs, 4)




