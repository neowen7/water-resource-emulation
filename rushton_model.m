%%%%%%%%%%%%%%%%% RIVER FLOW MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% References: Rushton, K.R., Eilers, V.H.M. and Carter, R.C., 2006. %%%%%
%%% Improved soil moisture balance methodology for recharge estimation. %%%
%%% Journal of Hydrology 318, 379-399. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 05/09/2017 version %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = rushton_model(input)

%% LOAD INPUT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data=load('basin_data.txt');           % basins characteristics 
weather=load('weather_data.txt');      % weather data
crop=load('crop_data.txt');            % growth stages for crop 
param=load('parameters.txt');          % best set of parameters from model calibration

%% SELECT WEATHER DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prec=weather(:,[1:3 4]);       % precipitation [mm]
temp=weather(:,[1:3 5:6]);     % maximum and minimum temperatures [°C]
sun=weather(:,[1:3 7]);        % total daily solar radiation [MJ/m2]
rh=weather(:,[1:3 8]);         % mean daily relative humidity [%]
ws=weather(:,[1:3 9]);         % mean daily wind speed [m/s]

%% EXTRACT WEATHER SERIES FOR MODEL SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
first_year=data(1,10);         % first year of simulation
last_year=data(1,11);          % last year of simulation
indp=find(prec(:,1)==1 & prec(:,2)==1 & prec(:,3)==first_year);
indpF=find(prec(:,1)==31 & prec(:,2)==12 & prec(:,3)==last_year);
prec=prec(indp:indpF,:);   
indt=find(temp(:,1)==1 & temp(:,2)==1 & temp(:,3)==first_year);
indtF=find(temp(:,1)==31 & temp(:,2)==12 & temp(:,3)==last_year);
temp=temp(indt:indtF,:);       
date=prec(:,1:3);
year=[first_year:last_year]';
nyear=size(unique(date(:,3)),1);
ndays=size(date,1);
P=prec(:,4);                    % precipitation series [mm]
Tmax=temp(:,4);                 % maximum temperature series [°C]
Tmin=temp(:,5);                 % minimum temperature series [°C]
sun=sun(indp:end,4);            % solar radiation series [MJ/m2]
rh=rh(indp:end,4);              % relative humidity series [-]
ws=ws(indp:end,4);              % wind speed series [m/s]
x1=[1:ndays]';                  % number of days of simulation

%% CATCHMENT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Area=data(1,2);        % catchment area [m2]
lat=data(1,3);         % latitude 
Z=data(1,4);           % elevation [m AOD]

%% LAND USE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Woodlands = input(1);               % percentage of woodlands area in the catchment
Arable = input(2);                  % percentage of arable area in the catchment
Pasture = input(3);                 % percentage of pasture area in the catchment
Impermeable = input(4);             % percentage of impermeable area in the catchment

%% GROUNDWATER AND SURFACE WATER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ga=0;    % groundwater abstraction [mm/day]
gd=0;    % groundwater discharges [mm/day]
swa=0;   % surface water abstractions [m3/s]
swd=0;   % surface water discharges [m3/s]
   
%% MODEL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fc=param(1,1);             % field capacity [m3/m3]
wp=param(1,2);             % wilting point  [m3/m3]
pRAW=param(1,3);           % depletion factor for RAW calculation [-]
rdg=param(1,4);            % root depth grass [m]
storestart=param(1,5);     % store start [mm]
sm=param(1,6);             % store max  [mm]
td=param(1,7);             % discharge time [day]
storage=param(1,8);        % storage parameter [mm]
gls=param(1,9);            % groundwater level start [m AOD]
bf=param(1,10);            % bypass fraction [-]
fracstor=param(1,11);      % factor for Near Surface Storage (NSS) [-]

%% POTENTIAL EVAPOTRANSPIRATION CALCULATION [PENMAN-MONTEITH] %%%%%%%%%%%%%%
Tmean=(Tmax+Tmin)/2;                    % mean temperature [°C]
Pa=101.325;                             % atmospheric pressure  [kPa]
Cp=1.005;                               % specific heat at constant pressure [kJ/kg°C]
lambda=2.51-2.361*(10^-3)*Tmean;        % latent heat of vaporization [kJ/kg],
albedo=0.15;                            % albedo 
sigma= 4.903*10^-9;                     % Stefan-Boltzman constant [MJ*m2*K^-4*giorno^-1]
kelvin=Tmean+273.15;                    % mean temperature [Kelvin]
kelvinTmax=Tmax+273.15;                 % maximum temperature [Kelvin]
kelvinTmin=Tmin+273.15;                 % minimum temperature [Kelvin]
dimW=size(Tmean,1);
G(1,1)=0;                               % soil heat flux density [MJ/m2*day]
gamma=Pa*0.000655;                      % psychrometric constant [kPa/°C],
J=[1:ndays]';                           % number of the day in the year between 1 (1 January) and 365 or 366 (31 December)
Dr=1+0.033*cos((2*pi*J/365));           % inverse relative distance Earth-Sun [rad]
deltab=0.409*sin((2*pi/365).*J-1.39);   % solar declination [rad]
phi=lat*pi/180;                         % latitude [rad]
omegas=acos(-tan(phi)*tan(deltab));     % sunset hour angle [rad]

for i=1:ndays
    delta(i,1)=4098*(0.6108*exp(17.27*Tmean(i,1)/(Tmean(i,1)+237.3)))/((237.3+Tmean(i,1))^2);
    esmax(i,1)=0.6108*(exp(17.27*Tmax(i,1)/(Tmax(i,1)+237.3)));
    esmin(i,1)=0.6108*(exp(17.27*Tmin(i,1)/(Tmin(i,1)+237.3)));
    es(i,1)=(esmax(i,1)+esmin(i,1))/2;      % saturation vapour pressure [kPa]
    ea(i,1)=es(i,1)*(rh(i,1));              % actual vapour pressure [kPa]
    Ra(i,1)=(24*60/pi)*0.082*(Dr(i,1))*(omegas(i,1)*sin(phi)*sin(deltab(i,1))+cos(phi)*cos(deltab(i,1))*sin(omegas(i,1))); % extraterrestrial radiation [MJ/m2*day]
    Rso(i,1)=(0.75+(2*10^-5)*Z)*Ra(i,1);
    Rns(i,1)=(1-albedo).*sun(i,1);
    Rnl(i,1)=(sigma*(kelvinTmax(i,1)^4+kelvinTmin(i,1)^4)/2)*(0.34-0.14*sqrt(ea(i,1)))*(1.35*(sun(i,1)/Rso(i,1))-0.35);
    Rn(i,1)= Rns(i,1)-Rnl(i,1);             % net radiation at the crop surface [MJ/m2*day]
end

for i=2:ndays
    G(i,1)=(0.38/lambda(i,1))*(kelvin(i,1)-kelvin(i-1,1));  % soil heat flux density [MJ/m2*day]
 end

for i=1:ndays
    % reference potential evapotranspiration [mm]
    ET0(i,1)=(0.408*delta(i,1)*(Rn(i,1)-G(i,1))+gamma*(900/(Tmean(i,1)+273))*ws(i,1)*(es(i,1)-ea(i,1)))/(delta(i,1)+gamma*(1+0.34*ws(i,1)));
    if ET0(i,1)<0
       ET0(i,1)=0;
    end
end

%% GRASS & CROP EVAPOTRANSPIRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ETgrass(:,1)=ET0;              % potential evapotranspiration for grass [mm]
Kc_crop(:,1:3)=date(:,1:3);    % crop coefficient for FAO actual evapotranspiration formula [-]
rd_crop(:,1:3)=date(:,1:3);    % root depth of crop [m]

for kk=1:size(crop)
        ind_k=find(Kc_crop(:,1)==crop(kk,1) & Kc_crop(:,2)==crop(kk,2));
        Kc_crop(ind_k,4)=crop(kk,3);   
        rd_crop(ind_k,4)=crop(kk,4);   
 end
 
 for k=1:ndays
     ETcrop(k,1)=Kc_crop(k,4)*ET0(k,1);  % potential evapotranspiration for crop [mm]
 end
  
%% MODEL INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SMDactual_grass=zeros(ndays+1,1);   % soil moisture deficit 
runCoeff=zeros(ndays,1);            % runoff coefficient at first time step [-]
SurfStor=zeros(ndays+1,1);          % surface storage
abstractionGW=ga*ones(ndays,1);     % abstraction from groundwater
dischargeGW=gd*ones(ndays,1);       % discharge to groundwater
abstractionSW=swa*ones(ndays,1);    % abstraction from surface water
dischargeSW=swd*ones(ndays,1);      % discharge to surface water
storeinit=zeros(ndays+1,1);         % store at first time step [mm]
   
%% GRASS WATER BALANCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TAW=(1000*(fc-wp)*rdg)*ones(size(date,1),1);    % total available water calculation 
RAW=pRAW.*TAW;                                  % readily available water calculation
   
for i=1:ndays
    runoff_grass(i,1)=P(i,1)*runCoeff(i,1);
    bypass(i,1)=max(0,(P(i,1)-runoff_grass(i,1))*bf);
    InfSoilZone(i,1)=P(i,1)-runoff_grass(i,1)-bypass(i,1)+SurfStor(i,1);
    if InfSoilZone(i,1)>ETgrass(i,1)
       SurfStor(i+1,1)=fracstor*(InfSoilZone(i,1)-ETgrass(i,1));
    else
      InfSoilZone(i+1,1)=0;
    end
    Ks(i,1)=(TAW(i,1)-SMDactual_grass(i,1))/(TAW(i,1)-RAW(i,1));
            
    % actual evapotranspiration calculation
    if SMDactual_grass(i)<RAW(i,1)
        AEgrass(i,1)=ETgrass(i,1);
    else
        if InfSoilZone(i,1)>ETgrass(i,1)
            AEgrass(i,1)=ETgrass(i,1);
        else
            if SMDactual_grass(i)>=RAW(i,1) && SMDactual_grass(i)<=TAW(i,1) && InfSoilZone(i,1)<ETgrass(i,1)
                AEgrass(i,1)=ETgrass(i,1)+Ks(i,1)*(ETgrass(i,1)-InfSoilZone(i,1));
            else
                if SMDactual_grass(i)>=TAW(i,1) && InfSoilZone(i,1)<ETgrass(i,1)
                    AEgrass(i,1)=InfSoilZone(i,1);
                else
                    AEgrass(i,1)=0;
                end
            end
        end
    end          
    
    SMDtemp_grass(i,1)=SMDactual_grass(i,1)-InfSoilZone(i,1)+SurfStor(i,1)+AEgrass(i,1);
    SMDactual_grass(i+1,1)=max(0,SMDtemp_grass(i,1));
    
% Rainfall and soil moisture deficit exceeds
if SMDactual_grass(i+1,1)>=0
    
    if P(i,1)>=20 && SMDactual_grass(i+1,1)<25
        runCoeff(i,1)=0.01;
    else
        if P(i,1)>=50 && SMDactual_grass(i+1,1)<25
            runCoeff(i,1)=0.05;
        else
            if P(i,1)>50 && SMDactual_grass(i+1,1)>=25
                runCoeff(i,1)=0.02;
            else
                runCoeff(i,1)=0;
            end
        end
    end
    
else
    runCoeff(i,1)=0;
end

if SMDtemp_grass(i,1)<=0
    recharge_grass(i,1)=-SMDtemp_grass(i,1)+bypass(i,1);
else
    recharge_grass(i,1)=bypass(i,1);
end

cumrecharge=cumsum(recharge_grass);

end

%% CROP WATER BALANCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
TAWcrop=1000*(fc-wp)*rd_crop(:,4);  % total available water calculation
RAWcrop=pRAW.*TAWcrop;              % readily available water calculation

% Runoff coefficient
SMDactual_crop=zeros(ndays+1,1);
runCoeff_crop=zeros(ndays,1);       
SurfStor_crop=zeros(ndays+1,1);

for i=1:ndays
    runoff_crop(i,1)=P(i,1)*runCoeff_crop(i,1);
    bypass_crop(i,1)=max(0,(P(i,1)-runoff_crop(i,1))*bf);
    InfSoilZone_crop(i,1)=P(i,1)-runoff_crop(i,1)-bypass_crop(i,1)+SurfStor_crop(i,1);
    if InfSoilZone_crop(i,1)>ETcrop(i,1)
       SurfStor_crop(i+1,1)=fracstor*(InfSoilZone_crop(i,1)-ETcrop(i,1));
    else
       InfSoilZone_crop(i+1,1)=0;
    end
    Ks_crop(i,1)=(TAWcrop(i,1)-SMDactual_crop(i,1))/(TAWcrop(i,1)-RAWcrop(i,1));  % evaporation stress coefficient
            
    % actual evapotranspiration calculation
    if SMDactual_crop(i)<RAWcrop(i,1)
        AEcrop(i,1)=ETcrop(i,1);
    else
        if InfSoilZone_crop(i,1)>ETcrop(i,1)
            AEcrop(i,1)=ETcrop(i,1);
        else
            if SMDactual_crop(i)>=RAWcrop(i,1) && SMDactual_crop(i)<=TAWcrop(i,1) && InfSoilZone_crop(i,1)<ETcrop(i,1)
                AEcrop(i,1)=ETcrop(i,1)+Ks_crop(i,1)*(ETcrop(i,1)-InfSoilZone_crop(i,1));
            else
                if SMDactual_crop(i)>=TAWcrop(i,1) && InfSoilZone_crop(i,1)<ETcrop(i,1)
                    AEcrop(i,1)=InfSoilZone_crop(i,1);
                else
                    AEcrop(i,1)=0;
                end
            end
        end
    end          
    
    SMDtemp_crop(i,1)=SMDactual_crop(i,1)-InfSoilZone_crop(i,1)+SurfStor_crop(i,1)+AEcrop(i,1);
    SMDactual_crop(i+1,1)=max(0,SMDtemp_crop(i,1));
    
% Rainfall and soil moisture deficit exceeds
if SMDactual_crop(i+1,1)>=0
    
    if P(i,1)>=20 && SMDactual_crop(i+1,1)<25
        runCoeff_crop(i,1)=0.01;
    else
        if P(i,1)>=50 && SMDactual_crop(i+1,1)<25
            runCoeff_crop(i,1)=0.05;
        else
            if P(i,1)>50 && SMDactual_crop(i+1,1)>=25
                runCoeff_crop(i,1)=0.02;
            else
                runCoeff_crop(i,1)=0;
            end
        end
    end
    
else
    runCoeff_crop(i,1)=0;
end       

if SMDtemp_crop(i,1)<=0
    recharge_crop(i,1)=-SMDtemp_crop(i,1)+bypass_crop(i,1);
else
    recharge_crop(i,1)=bypass_crop(i,1);
end

cumrecharge_crop=cumsum(recharge_crop);

end

%% DAILY WATER BALANCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
recharge=(Arable/100)*recharge_crop+(Pasture/100)*recharge_grass+(Woodlands/100)*recharge_grass+(Impermeable/100)*P;    % groundwater recharge 
runoff=(Arable/100)*runoff_crop+(Pasture/100)*runoff_grass+(Woodlands/100)*runoff_grass;                                % runoff
storeinit(1,1)=storestart;
baseflow(1,1)=max(0,(storeinit(1,1)/td));
if (storeinit(1,1)-baseflow(1,1)-abstractionGW(1,1))>sm
    overtop(1,1)=storeinit(1,1)-baseflow(1,1)-abstractionGW(1,1)-sm;
else
    overtop(1,1)=0;
end
storend(1,1)=recharge(1,1)+storeinit(1,1)+dischargeGW(1,1)-abstractionGW(1,1)-baseflow(1,1)-overtop(1,1);
changeGWlevel(1,1)=(storend(1,1)-storeinit(1,1))/storage;
GWlevel(1,1)=gls+changeGWlevel(1,1)/1000;
naturalised_flow(1,1)=runoff(1,1)+baseflow(1,1)+overtop(1,1);
total_river_flow(1,1)=naturalised_flow(1,1)+abstractionSW(1,1)+dischargeSW(1,1);
riverflow(1,1)=(total_river_flow(1,1)/1000)*(Area*1000000)/86400;

for w=2:ndays
    storeinit(w,1)=storend(w-1,1);
    baseflow(w,1)=max(0,(storeinit(w,1)/td));
    if (storeinit(w,1)-baseflow(w,1)-abstractionGW(w,1))>sm
        overtop(w,1)=storeinit(w,1)-baseflow(w,1)-abstractionGW(w,1)-sm;
    else
        overtop(w,1)=0;
    end
    storend(w,1)=recharge(w,1)+storeinit(w,1)+dischargeGW(w,1)-abstractionGW(w,1)-baseflow(w,1)-overtop(w,1);
    changeGWlevel(w,1)=(storend(w,1)-storeinit(w,1))/storage;
    GWlevel(w,1)=gls+changeGWlevel(w,1)/1000;
    naturalised_flow(w,1)=runoff(w,1)+baseflow(w,1)+overtop(w,1);
    total_river_flow(w,1)=naturalised_flow(w,1)+abstractionSW(w,1)+dischargeSW(w,1);
    riverflow(w,1)=(total_river_flow(w,1)/1000)*(Area*1000000)/86400;      % modelled river flow series [m3/s]
end

out=riverflow;  % modelled river flow series [m3/s]

