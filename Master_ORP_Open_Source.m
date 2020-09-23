%Master ORP Growth + Thermal Model for Open Source Model
clc
clear

%Import data for T_amb, GHI, RH, and WNDSPD from TMY3 Dataset
TMY3_data = xlsread('C:\CSU_Open_Source_Algae_Model\Model_Files\CSU_Algae_Model_User_Interface.xlsm', 3, 'E3:AU17522');
T_amb= (TMY3_data(:,28)+273.15);
GHI= TMY3_data(:,1);
RH= TMY3_data(:,34);
WNDSPD= TMY3_data(:, 43);

%Read Inputs from Excel
Inputs_Data = xlsread('C:\CSU_Open_Source_Algae_Model\Model_Files\CSU_Algae_Model_User_Interface.xlsm', 1, 'C2:B31');

%Assign values to variables within code

%ORP Geometry
length = Inputs_Data(5,:);
width = Inputs_Data(6,:);
depth = Inputs_Data(7,:);
num_ponds = Inputs_Data(8,:);

%Strain Parameters
T_opt = Inputs_Data(14,:);
T_min = Inputs_Data(15,:); 
T_max = Inputs_Data(16,:);
ODC = Inputs_Data(17,:);
I_sat = Inputs_Data(18,:);
night_resp= Inputs_Data(19,:);

%Harvesting Parameters
IBC = Inputs_Data(20,:);
harvest_conc = Inputs_Data(21,:);
harvest_days = Inputs_Data(22,:);

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

%Define unchanging parameters
h = 3600; %Time step is 1 hour
pi = 3.14159;
rho_algae = 1000.0;  %kg/m3
cp_algae = 4184.0; 
Conversion = 4.56 * 10^-6 ; %GHI to mol/m2*s
area = pi*0.25*(width^2) + width*(length-width);
volume = area*depth; 
perimeter = (3.14*(width/2)^2)+(2*(length-width)); 

%Nutrient Parameters
Ammonia_N_Content = 0.82; 
DAP_N_Content = 0.18;
DAP_P_Content = 0.20; 

%Initialize Variables
Harvest_Mass = 0.0; 
marker = 0;
marker_crash = 0; 
op_hours = 0;
Freezes = 0; 

%Prelocate all variables that change in the for loop to save time
TR = zeros(100000, 1);
CX = zeros(100000, 1);
Q_Atmospheric = zeros(100000, 1);
M_Evapo = zeros(100000, 1);
Q_Evapo = zeros(100000, 1);
Q_Direct_Solar = zeros(100000,1);
Q_Ground_Conduction = zeros(100000,1);
Q_Reradiation = zeros(100000,1);
Q_Convection_track = zeros(100000,1); 
Q_Makeup_Water = zeros(100000,1); 
Water_consump = zeros(100000,1); 
Temp_eff = zeros(100000,1);
Conc_eff = zeros(100000,1);
Light_eff = zeros(100000,1);
decay = zeros(100000,1);
Conc_at_Harvest = zeros(100000,1);
CO2_demand_hrly = zeros(100000,1);
Dried_Harvest = zeros(100000,1);
Harvest_Shortage = zeros(100000,1);
reliability = zeros(100000,1);
x = zeros(100000,1);
Crashes = zeros(100000,1);
Ammonia_demand_hrly = zeros(100000, 1); 
DAP_demand_hrly = zeros(100000, 1); 
Potash_demand_hrly = zeros(100000, 1); 

%% This is the start of the forloop 

for i = n1:(n2-1) %Only go to n-1 because program reads to i+1
 
    if (i == n1)
        TR(i,1) = T_amb(i,1);
        CX(i,1) = IBC; 
    elseif (CX(i,1)==0)
        CX(i,1) = IBC;
    else
        TR(i,1) = TR(i,1);
        CX(i,1) = CX(i,1);
    end

%Biological Growth Model
[Temp_eff(i,1)] = Temp_Efficiency_ORP (TR(i,1), T_opt, T_min, T_max);
[Conc_eff(i,1)] = Concentration_Efficiency_ORP(CX(i,1), ODC, depth);
[Light_eff(i,1)] = Light_Efficiency_ORP(GHI(i,1), Conc_eff(i,1), I_sat);   
[decay(i,1)] = Night_Respiration_ORP(GHI(i,1), CX(i,1), night_resp, volume, Temp_eff(i,1));
dCXdt = (12/8)*Light_eff(i,1)*Temp_eff(i,1)*GHI(i,1)*0.458*0.95*Conversion*area + decay(i,1); 
CX(i+1,1) = ((dCXdt*h)/volume + CX(i,1));


%Thermal Model - 4th Order Runge-Kutta Scheme
%Compute Air Properties at this time step using the film temperature
[rho_air, ~, ~, ~, dyn_visc_air, ~, ~]=Air_Properties((TR(i,1)+T_amb(i,1))/2);
%Compute thermal mass based on algal properties
Therm_mass = volume*rho_algae*cp_algae;

%Fist RK4 Fluxes
[Q_Solar] = Direct_Solar_HT_ORP(GHI(i,1)); 
[Q_Convection] = Convection_HT_ORP(area, perimeter, TR(i,1), T_amb(i,1), WNDSPD(i,1), dyn_visc_air, rho_air); 
[Q_Longwave_Atmo] = Longwave_Atmo_HT_ORP(T_amb(i,1)); 
[Q_Ground] = Ground_Conduction_HT_ORP(TR(i,1));
[M_Evap, Q_Evap]=Evaporative_HT_ORP(area, perimeter, TR(i,1), T_amb(i,1), RH(i,1), WNDSPD(i,1), dyn_visc_air, rho_air); 
[Q_Inflow] = Inflow_HT_ORP(M_Evap, cp_algae, TR(i,1), T_amb(i,1)); 
[Q_Rerad] = Reradiation_HT_ORP(TR(i,1)); 

%RK4_1 time derivative of temperature - K1=dT1
dT1 =(Q_Evap + Q_Convection +Q_Inflow + Q_Longwave_Atmo + Q_Solar + Q_Ground + Q_Rerad)*area/Therm_mass;
T1 = TR(i,1) + 0.5*h*dT1;

%Air Properties - Thermal mass already defined
[rho_air, ~, ~, ~, dyn_visc_air, ~, ~]=Air_Properties((T1+((T_amb(i,1)+T_amb(i+1,1))/2))/2);

%Second RK4 Fluxes
[Q_Solar] = Direct_Solar_HT_ORP((GHI(i,1)+GHI(i+1,1))/2); 
[Q_Convection] = Convection_HT_ORP(area, perimeter, T1, ((T_amb(i,1)+T_amb(i+1,1))/2), (WNDSPD(i,1)+WNDSPD(i+1,1))/2, dyn_visc_air, rho_air); 
[Q_Longwave_Atmo] = Longwave_Atmo_HT_ORP((T_amb(i,1)+T_amb(i+1,1))/2);  
[Q_Ground] = Ground_Conduction_HT_ORP(T1); 
[M_Evap, Q_Evap] = Evaporative_HT_ORP(area, perimeter, T1, ((T_amb(i,1)+T_amb(i+1))/2), (RH(i,1)+RH(i+1,1))/2, (WNDSPD(i,1)+WNDSPD(i+1,1))/2, dyn_visc_air, rho_air); 
[Q_Inflow] = Inflow_HT_ORP(M_Evap, cp_algae, T1, (T_amb(i,1)+T_amb(i+1, 1))/2);
[Q_Rerad] = Reradiation_HT_ORP(T1); 

%RK4_2 time derivative - K2 = dT2
dT2 = (Q_Evap +  Q_Convection + Q_Inflow + Q_Longwave_Atmo + Q_Solar + Q_Ground + Q_Rerad)*area/Therm_mass;
T2= TR(i,1) + 0.5*h*dT2; 

%Air Properties - Thermal mass already defined
[rho_air, ~, ~, ~, dyn_visc_air, ~, ~] = Air_Properties((T2+((T_amb(i,1)+T_amb(i+1,1))/2))/2); 

%Third RK4 Fluxes
[Q_Solar] = Direct_Solar_HT_ORP((GHI(i,1)+GHI(i+1,1))/2); 
[Q_Convection] = Convection_HT_ORP(area, perimeter, T2, (T_amb(i,1)+T_amb(i+1,1))/2, (WNDSPD(i,1)+WNDSPD(i+1,1))/2, dyn_visc_air, rho_air); 
[Q_Longwave_Atmo] = Longwave_Atmo_HT_ORP((T_amb(i,1)+T_amb(i+1,1))/2);  
[Q_Ground] = Ground_Conduction_HT_ORP(T2); 
[M_Evap, Q_Evap] = Evaporative_HT_ORP(area, perimeter, T2,(T_amb(i,1)+T_amb(i+1,1))/2,(RH(i,1)+RH(i+1,1))/2, (WNDSPD(i,1)+WNDSPD(i+1,1))/2, dyn_visc_air, rho_air); 
[Q_Inflow] = Inflow_HT_ORP(M_Evap, cp_algae, T2, (T_amb(i,1)+T_amb(i+1, 1))/2);
[Q_Rerad] = Reradiation_HT_ORP(T2);

%RK4_3 time derivative - K3 = dT3
dT3 = (Q_Evap + Q_Convection + Q_Inflow + Q_Longwave_Atmo + Q_Ground + Q_Solar + Q_Rerad)*area/Therm_mass;
T3 = TR(i,1) + h*dT3; 

%Air Properties - Thermal mass already defined
[rho_air, ~, ~, ~, dyn_visc_air, ~, ~] = Air_Properties((T3+ T_amb(i+1,1))/2); 

%Fourth RK4 Fluxes
[Q_Solar] = Direct_Solar_HT_ORP(GHI(i+1,1)); 
[Q_Convection] = Convection_HT_ORP(area, perimeter, T3, T_amb(i+1,1), WNDSPD(i+1,1), dyn_visc_air, rho_air); 
[Q_Longwave_Atmo] = Longwave_Atmo_HT_ORP(T_amb(i+1,1));  
[Q_Ground] = Ground_Conduction_HT_ORP(T3); 
[M_Evap, Q_Evap] = Evaporative_HT_ORP(area, perimeter, T3, T_amb(i+1,1), RH(i+1,1), WNDSPD(i+1,1), dyn_visc_air, rho_air); 
[Q_Inflow] = Inflow_HT_ORP(M_Evap, cp_algae, T3, T_amb(i+1, 1));
[Q_Rerad] = Reradiation_HT_ORP(T3); 

%RK4_1 time derivative  
dTdT = (Q_Evap + Q_Convection + Q_Inflow + Q_Longwave_Atmo + Q_Ground + Q_Solar + Q_Rerad)*area/Therm_mass;
TR(i+1,1)= TR(i,1) + (1/6)*h*(dT1+2.0*(dT2+dT3)+dTdT);


%Marker and Operational Hours Counting
marker = marker + 1;
op_hours = op_hours + 1;

%Harvesting Sequence
        if ((mod(marker, (harvest_days*24.0))== 0)||(CX(i+1,1) >= harvest_conc)||(i==(n2-1)))           
        Harvest_Mass = Harvest_Mass + (CX(i+1,1)-IBC)*volume;
        Conc_at_Harvest(i,1) = CX(i+1,1);
        CX(i+1,1) = IBC;
        marker = 0;
        else
        Harvest_Mass = Harvest_Mass +0; 
        Conc_at_Harvest(i,1) = 0;
        end 
        
%Track Thermal Fluxes with each Iteration
[rho_air, ~, ~, ~, dyn_visc_air, ~, ~] = Air_Properties((TR(i,1)+T_amb(i,1))/2);
Q_Atmospheric(i,1) = Longwave_Atmo_HT_ORP(T_amb(i,1));
[M_Evapo(i,1), Q_Evapo(i,1)] = Evaporative_HT_ORP(area, perimeter, TR(i,1), T_amb(i,1), RH(i,1), WNDSPD(i,1), dyn_visc_air, rho_air);
Q_Makeup_Water(i,1) = Inflow_HT_ORP(M_Evapo(i,1), cp_algae, TR(i,1), T_amb(i, 1));
Q_Direct_Solar(i,1) = Direct_Solar_HT_ORP(GHI(i,1)); 
Q_Ground_Conduction(i,1) = Ground_Conduction_HT_ORP(TR(i+1,1)); 
Q_Reradiation(i,1) = Reradiation_HT_ORP(TR(i+1,1));
Q_Convection_track(i,1) = Convection_HT_ORP(area, perimeter, TR(i+1,1), T_amb(i,1), WNDSPD(i,1), dyn_visc_air, rho_air);
        
    
%Calculate waterloss at each time step 
Water_consump(i,1) = M_Evapo(i,1)*area*h; %[kg/m2*s]*[m2]*[3600 s] = [kg/hr]

    if Water_consump(i,1) < 0
       Water_consump(i,1) = 0;
    else 
       Water_consump(i,1) = Water_consump(i,1); 
    end    
        
        
%Calculate the CO2 demand kg/hr based on stoichiometric carbon balance 
    if (i == n1)
        CO2_demand_hrly(i,1) = 0;
    else
        CO2_demand_hrly(i,1) = (((CX(i+1,1)-CX(i,1))*volume/1000)*Mass_frac_CO2)/CO2_util;
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
        DAP_demand_hrly(i,1) = (((CX(i+1,1)-CX(i,1))*volume/1000)*Algae_P_Content/DAP_P_Content); %g/m3 * m3 * 1 kg/1000 g * % 
        Ammonia_demand_hrly(i,1) = (((CX(i+1,1)-CX(i,1))*volume/1000)*Algae_N_Content/Ammonia_N_Content)-(DAP_demand_hrly(i,1)*DAP_N_Content);
    end

    if  Ammonia_demand_hrly(i,1) < 0
        Ammonia_demand_hrly(i,1) = 0;
    else   
        Ammonia_demand_hrly(i,1) = Ammonia_demand_hrly(i,1);
    end
    
    if  DAP_demand_hrly(i,1) < 0
        DAP_demand_hrly(i,1) = 0;
    else   
        DAP_demand_hrly(i,1) = DAP_demand_hrly(i,1);
    end
        
end 

%% Results and formatting

%Summary of CO2 Consumption
CO2_consump_total = sum(CO2_demand_hrly(n1:n2));
CO2_consump_daily = CO2_consump_total/((n2-n1)/24); 

%Temporal Outputs
T_amb_Plot = T_amb(n1:n2, 1);
TR_Plot=TR(n1:n2,1); 
CX_plot = CX(n1:n2, 1); 
CO2_hourly_kg = CO2_demand_hrly(n1:n2,1); 
Water_con = Water_consump(n1:n2,1);
Ammonia_con = Ammonia_demand_hrly(n1:n2,1); 
DAP_con = DAP_demand_hrly(n1:n2,1); 

%Plot Temporally Resolved Outputs into Excel
Temporal_Outputs = [T_amb_Plot, TR_Plot, CX_plot, CO2_hourly_kg, Water_con, Ammonia_con, DAP_con];
xlswrite('C:\CSU_Open_Source_Algae_Model\Model_Files\Growth_Model_Outputs.xlsx', Temporal_Outputs, 2);

%Single Number Outputs
Areal_Prod = Harvest_Mass/(area*((n2-n1)/(24.0)));
Vol_Prod = Harvest_Mass/(volume*((n2-n1)/(24.0))*1000);
CO2_consump_ave = CO2_consump_daily; 
Water_consump_ave = sum(Water_con)/((n2-n1)/24); %kg per pond per day
Ammonia_consump_ave = sum(Ammonia_con)/((n2-n1)/24);
DAP_consump_ave = sum(DAP_con)/((n2-n1)/24); 
Conc_at_Harvest_Plot = Conc_at_Harvest(n1:n2,1);
Conc_at_Harvest_Ave = mean(nonzeros(Conc_at_Harvest_Plot)); 
Evap_rate = Water_consump_ave*(1/1000)*(1/area)*100; 
Op_Days_per_Year = (n2-n1)/24;

%Write Single Number Outputs to Excel
Growth_Model_Static_Outputs = [Areal_Prod; Vol_Prod; Water_consump_ave; CO2_consump_ave; Ammonia_consump_ave; DAP_consump_ave; Op_Days_per_Year; Conc_at_Harvest_Ave; Evap_rate]; 
xlswrite('C:\CSU_Open_Source_Algae_Model\Model_Files\Growth_Model_Outputs.xlsx', Growth_Model_Static_Outputs, 'B3:B11')

%Plot Temporally Resolved Thermal Fluxes
Q_Atmospheric_Plot = Q_Atmospheric(n1:n2-1,1);
Q_Direct_Solar_Plot = Q_Direct_Solar(n1:n2-1,1);
Q_Reradiation_Plot = Q_Reradiation(n1:n2-1,1);
Q_Ground_Conduction_Plot = Q_Ground_Conduction(n1:n2-1,1); 
Q_Evaporative_Plot = Q_Evapo(n1:n2-1,1);
Q_Convection_Plot = Q_Convection_track(n1:n2-1,1);
Q_Makeup_Water_Plot = Q_Makeup_Water(n1:n2-1,1);

Thermal_Flux_Outputs = [Q_Atmospheric_Plot, Q_Direct_Solar_Plot, Q_Reradiation_Plot, Q_Ground_Conduction_Plot, Q_Evaporative_Plot, Q_Convection_Plot, Q_Makeup_Water_Plot];
xlswrite('C:\CSU_Open_Source_Algae_Model\Model_Files\Growth_Model_Outputs.xlsx', Thermal_Flux_Outputs, 3)


