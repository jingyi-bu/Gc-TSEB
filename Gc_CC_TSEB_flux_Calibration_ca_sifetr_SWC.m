function [LE,H,Rn,G_dry,LEc_dry,LEs_dry,Qc,Tc_dry,Ts_dry,Tsurf_new,Rnc_dry,Rns_dry,interation_total,ra,rs] = Gc_CC_TSEB_flux_Calibration_ca_sifetr_SWC(g1,etr_a,apm1,bpm1,alphag1,omi)
%Soil water content was taken into consideration
%SIF
%% Paramters
rou=1.225;%Air density
cp=1005;% Atmospheric specific heat
a=6.11; % unit:hpa
b=17.502;%Basic parameters of saturated water vapor press
c=240.97;
gama=0.667; %Psychrometric constant 0.667hpa/K
sigma=5.67e-8;   %Stefan-Bolzmann constant
emis_c=0.98;  %canopy emissivity
emis_s=0.95;  %soil emissivity



% Data input
% date_temp=evalin('base','date_temp');


lai=evalin('base','lai');
lai_plus_stem=evalin('base','lai_plus_stem');
albedo=evalin('base','albedo');
 
air_temp=evalin('base','air_temp');
air_pres=evalin('base','air_pres');
air_humi=evalin('base','air_humi');
wind=evalin('base','wind');
incoming_short_radiation=evalin('base','incoming_short_radiation');
incoming_long_radiation=evalin('base','incoming_long_radiation');
Da=evalin('base','Da');
% LE_mea=evalin('base','LE_mea');
% H_mea=evalin('base','H_mea');
% Tsurf_mea=evalin('base','Tsurf_mea');

sif=evalin('base','sif');
SWC=evalin('base','SWC');
f_dry=evalin('base','f_dry');

hc=evalin('base','hc');        % canopy height
d0=evalin('base','d0');        % zero plane displacement
zm0=evalin('base','zm0');      % the  length of momentum transfer roughness
zh0=evalin('base','zh0');      % the  length of heat transfer roughness
zm=evalin('base','zm');        % the observation height of wind speed
zh=evalin('base','zh');        % the observation height of air temperature
J=evalin('base','J');        %date index
lontitude=evalin('base','lontitude');
lat=evalin('base','lat');
time=evalin('base','time');  % hour
 

Ca=evalin('base','Ca');
Qf=evalin('base','Qf');
Qw=evalin('base','Qw');
c4=evalin('base','c4');

 
n=length(air_temp);
alphag=alphag1.*hc./hc;
apm=apm1.*hc./hc;
bpm=bpm1.*hc./hc;
 
%% Attenuation index for net radiation
KA=0*air_temp-9999.0;
index_lai=find((lai<1) & (lai>0));
KA(index_lai)=0.8;
index_lai=find((lai<2) & (lai>=1));
KA(index_lai)=0.6;
index_lai=find(lai>=2);
KA(index_lai)=0.45;
 
 
%% Calculation of solar zenith angle
% J=0*air_temp;%day of year
% % length(J)
% for i=1:1:length(J)
%     test=num2str(date_temp(i,1)); 
%     J(i)=str2num(test(5:7));
% end
 
%The sine of the declination for the sun
sin_declination=0.39785*sin((278.97+0.9856.*J+1.9165*sin((356.6+0.9856.*J)*pi/180))*pi/180);% sine of solar declination, solar declination varies from -23.45(summer) to 23.45(winter)
 
f=(279.575+0.9856*J)*pi/180;
ET=596.2*sin(2*f)-104.7*sin(f)+4.3*sin(3*f)-12.7*sin(4*f)-429.3*cos(f)-2*cos(2*f)+19.3*cos(3*f);
ET=ET/3600;
LC=lontitude/15;
t0=-LC-ET+12;
%
 
% time=0*air_temp;
% for i=1:n
%     test=num2str(date_temp(i));
%     if length(test)>10
%         time(i)=str2num(test(9:10))+str2num(test(11))/6; 
%
%     else
%         time(i)=str2num(test(9:10));
%     end
% end
 
cos_solar_zenith=0*air_temp;   %cos_solar_zenith
cos_solar_zenith=sin(lat).*sin_declination+cos(lat).*abs(sqrt(1-sin_declination.^2)).*cos(15.*(time-t0)*pi/180);
fc=1-exp(-KA.*lai./sqrt(2*cos_solar_zenith)); %fc is used to partition the net radiation (calculating the distribution of net radiation between canopy and soil)
fr=1-exp(-0.5.*lai);%fr is used to synthesize LST
fwet=0*fc;
index_wet=find(air_humi>70);
fwet(index_wet)=(air_humi(index_wet)/100).^4;
fwet=0*fc;
 
%% Determination of resistances

%aerodynamic resistance
ra=(log((zm-d0)./zm0)).*(log((zh-d0)./zh0))./(0.16*wind);
%subcanopy resistance
ux=0.4*wind./(log((zm-d0)./zm0));     %friction wind speed
cs_bare=3.0769.*(0.01.*ux./(1.5/100000)).^(-0.45);   %3.0769=k/a
cs_dense=0.004+0*air_temp;
cs=cs_bare.*exp(-1.*lai-omi.*hc)+cs_dense.*(1-exp(-1.*lai-omi.*hc));
rs=1./(cs.*ux);
 
%boundary layer resistance //(L+S)
rb=20*(ux).^(-0.5)./lai_plus_stem;    % ux is the friction wind speed, and lai_plus_stem is the sum of the leaf area index and the stem area index.
 
%canopy conductance Qc
 
es=a*exp(b*(air_temp-273.15)./(air_temp-273.15+c)); %es is the saturated vapor pressure(unit hpa).
ea=es.*air_humi;%atmospheric water vapor pressure(unit hpa), air_humi is relative humidity!
delta=b*c*es./(c+air_temp-273.15).^2;
delta_gama=delta./(delta+0.667);
% Da=(es-ea)/10;  %Da is the water vapor pressure saturation difference(unit kpa)

%calculate Ci
air_t_c=air_temp-273.15;
o_2=2.09*10^4; %Pa
psi_cf=700;   %mol mol-1
omi=5.*o_2./(2600.*(0.57.^((air_t_c-25)./10))); %ppm
% DD=Da./(air_pres./100);  %%BB
DD=Da./(air_pres./10);
Ci=(1-sqrt(1.6.*DD.*(Ca-omi)./(psi_cf.*(Ca.^2)))).*Ca;

% etr_a=500;
ETR=etr_a.*sif;
EUE_c4=1/5;  %C4
EUE_c3=(1/5.8).*((Ci-omi)./(Ci+2.*omi));   %C3
EUE=EUE_c4.*c4+EUE_c3.*(1-c4);

GPP=ETR.*EUE;  

fSW=(SWC-Qw)./(Qf-Qw);   
fSW(find(SWC<Qw))=0;
fSW(find(SWC>Qf))=1;
Gc=1.6.*(1+g1.*fSW./((Da).^0.5)).*GPP./Ca;

Qc=Gc./41.03;
rc=1./Qc;
emis_land=0.98*fc+0.95*(1-fc); %The land surface emissivity is estimated roughly, because the surface emissivity is not very important in this model
emis_air=1.24*(ea./air_temp).^(1/7);%air emissivity based on ea(hpa) and air_temp(k)
 
trans_long=exp(-0.95.*lai); %Transmission of longwave radiation through the canopy
 
%%
Tsurf=0*air_temp-9999.0;
Tc_dry=0*air_temp-9999.0;
Ts_dry=0*air_temp-9999.0;
Tc_wet=0*air_temp;
Ts_wet=0*air_temp;
 
Hc_dry=0*air_temp-9999.0;
Hs_dry=0*air_temp-9999.0;
Hc_wet=0*air_temp;
Hs_wet=0*air_temp;
H=0*air_temp-9999.0;
LE=0*air_temp-9999.0;
LEc_dry=0*air_temp-9999.0;
LEs_dry=0*air_temp-9999.0;
LEc_wet=0*air_temp;
LEs_wet=0*air_temp;
 
LEp_mt=0*air_temp-9999.0; %the potential evaporation of the soil is calculated based on mass transfer
LEp_pm=0*air_temp-9999.0; %the potential evaporation of the soil is calculated based on PM
fpx=0*air_temp-9999.0;
fp=0*air_temp-9999.0;
Rnc_dry=0*air_temp-9999.0;
Rns_dry=0*air_temp-9999.0;
Rnc_wet=0*air_temp;
Rns_wet=0*air_temp;
Rn=0*air_temp-9999.0;
G_dry=0*air_temp-9999.0;
G_wet=0*air_temp;
 
etc_dry=0*air_temp-9999.0; %the saturated water vapor pressure for canopy
ets_dry=0*air_temp-9999.0; %the saturated water vapor pressure for soil, calculated by Ts
Tsurf_old=air_temp-2.1; %initial value of LST
Tsurf_new=air_temp; %initial value of LST
energy_error_soil=0*air_temp+100;
 
 
interation1=0;  %important!!! The interation1 must be 0!!!
 
%index is the point that satisfies the calculation condition
index=find((incoming_short_radiation>0) & (incoming_long_radiation>0) & (air_temp>0) & (air_humi>0) & (air_pres>0) & (wind>0)  & (albedo>0));
length(index);
temp_plus=zeros(length(air_temp),1);
temp_plus(index)=1;%the point that satisfies the calculation condition is 1, and the point that does not meet the calculation condition is 0.
index_wet=find((temp_plus>0) & (fwet>0)); %There are points that exist intercept evaporation.
 
%The model does not take into account the case of wet (that is, the relevant term of wet is 0).
 
%% The decomposition of net radiation
while length(index)>0 & max(abs(Tsurf_old(index)-Tsurf_new(index)))>0.1 & interation1<20
    interation1=interation1+1;
    
    Tsurf_old(index)=Tsurf_new(index);
    
    if interation1==1
        Rn(index)=(1-albedo(index)).*incoming_short_radiation(index)+emis_land(index).*incoming_long_radiation(index)-sigma*emis_land(index).*Tsurf_old(index).^4;
        Rns_dry(index)=Rn(index).*(1-fc(index)).*(1-fwet(index));
        Rnc_dry(index)=Rn(index).*fc(index).*(1-fwet(index));
        
    else
        [Rnc_dry(index),Rns_dry(index)] = radiation_Campbell(cos_solar_zenith(index),lai(index),incoming_long_radiation(index),incoming_short_radiation(index),Tc_dry(index),Ts_dry(index));
        Rn(index)=Rnc_dry(index)+Rns_dry(index);
        
        
    end
    
    Rns_wet(index)=Rn(index).*(1-fc(index)).*fwet(index);
    Rnc_wet(index)=Rn(index).*fc(index).*fwet(index);
    G_dry(index)=alphag(index).*Rns_dry(index);
    G_wet(index)=alphag(index).*Rns_wet(index);
    
    %% Calculation of the sensible and latent heat flux
    
    %% The energy balance of canopy
    
    Tc_dry(index)=Tsurf_old(index);
    f_Tc_dry=0*air_temp-9999.0;
    derivitive_Tc_dry=0*air_temp-9999.0;
    index_f_Tc_dry=find((temp_plus>0) & abs(f_Tc_dry)>5); %The points where the basic calculation condition is met and the air temperature is normal (not 9999).
    interation3=1;
    
    while length(index_f_Tc_dry)>0 & interation3<20
        interation3=interation3+1;
        %Newton's method
        %f(x)
        f_Tc_dry(index_f_Tc_dry)=rou*cp*(Tc_dry(index_f_Tc_dry)-air_temp(index_f_Tc_dry))./(ra(index_f_Tc_dry)+rb(index_f_Tc_dry))+rou*cp*(a*exp(b*(Tc_dry(index_f_Tc_dry)-273.15)./(Tc_dry(index_f_Tc_dry)-273.15+c))-ea(index_f_Tc_dry))./gama./(ra(index_f_Tc_dry)+(rc(index_f_Tc_dry)))-Rnc_dry(index_f_Tc_dry);
        %f'(x)
        
        if interation1==1
            derivitive_Tc_dry(index_f_Tc_dry)=rou*cp./(ra(index_f_Tc_dry)+rb(index_f_Tc_dry))+rou*cp*a*b*c*exp(b*(Tc_dry(index_f_Tc_dry)-273.15)./(Tc_dry(index_f_Tc_dry)-273.15+c))./(gama*(ra(index_f_Tc_dry)+(rc(index_f_Tc_dry))).*(c+Tc_dry(index_f_Tc_dry)-273.15).^2);
            
        else
            derivitive_Tc_dry(index_f_Tc_dry)=rou*cp./(ra(index_f_Tc_dry)+rb(index_f_Tc_dry))+rou*cp*a*b*c*exp(b*(Tc_dry(index_f_Tc_dry)-273.15)./(Tc_dry(index_f_Tc_dry)-273.15+c))./(gama*(ra(index_f_Tc_dry)+(rc(index_f_Tc_dry))).*(c+Tc_dry(index_f_Tc_dry)-273.15).^2)-(1-trans_long(index_f_Tc_dry)).*(-8.*emis_c.*sigma.*Tc_dry(index_f_Tc_dry).^3);
            
        end
        %x=x-f(x)/f'(x)
        Tc_dry(index_f_Tc_dry)=Tc_dry(index_f_Tc_dry)-f_Tc_dry(index_f_Tc_dry)./derivitive_Tc_dry(index_f_Tc_dry);
        
        if interation1==1
            
        else
            [Rnc_dry(index),Rns_dry(index)] = radiation_Campbell(cos_solar_zenith(index),lai(index),incoming_long_radiation(index),incoming_short_radiation(index),Tc_dry(index),Ts_dry(index));
        end
        
        f_Tc_dry(index_f_Tc_dry)=rou*cp*(Tc_dry(index_f_Tc_dry)-air_temp(index_f_Tc_dry))./(ra(index_f_Tc_dry)+rb(index_f_Tc_dry))+rou*cp*(a*exp(b*(Tc_dry(index_f_Tc_dry)-273.15)./(Tc_dry(index_f_Tc_dry)-273.15+c))-ea(index_f_Tc_dry))./gama./(ra(index_f_Tc_dry)+(rc(index_f_Tc_dry)))-Rnc_dry(index_f_Tc_dry);
        
        
        
        % Re-iterative calculation of points with excessive closure difference
        index_f_Tc_dry=find((temp_plus>0) & abs(f_Tc_dry)>5);
        
    end
    
    
    LEc_wet(index_wet)=1.26*delta_gama(index_wet).*Rnc_wet(index_wet);
    Hc_wet(index_wet)=Rnc_wet(index_wet)-LEc_wet(index_wet);
    Tc_wet(index_wet)=Hc_wet(index_wet).*(ra(index_wet)+rb(index_wet))./rou./cp+air_temp(index_wet);
    
    
    %% The energy balance of soil
    Ts_dry(index)=Tsurf_old(index);
    f_Ts_dry=0*air_temp-9999.0;
    derivitive_Ts_dry=0*air_temp-9999.0;
    index_f_Ts_dry=find((temp_plus>0) & abs(f_Ts_dry)>5);
    interation4=1;

            while length(index_f_Ts_dry)>0 & interation4<20
        interation4=interation4+1;
        ets_dry(index_f_Ts_dry)=a*exp(b*(Ts_dry(index_f_Ts_dry)-273.15)./(Ts_dry(index_f_Ts_dry)-273.15+c)); % saturated water vapor pressure(hpa)
        LEp_mt(index_f_Ts_dry)=rou*cp*(ets_dry(index_f_Ts_dry)-ea(index_f_Ts_dry))./(ra(index_f_Ts_dry)+rs(index_f_Ts_dry))./gama;
        LEp_pm(index_f_Ts_dry)=(delta(index_f_Ts_dry).*((1-alphag(index_f_Ts_dry)).*Rns_dry(index_f_Ts_dry))+rou*cp*(es(index_f_Ts_dry)-ea(index_f_Ts_dry))./(ra(index_f_Ts_dry)+rs(index_f_Ts_dry)))./delta_gama(index_f_Ts_dry);
        % If LEp_mt<LEp_pm, LEp_mt is replaced by LEp_pm
        index_small=find(LEp_mt<LEp_pm);
        LEp_mt(index_small)=LEp_pm(index_small);
        
        fpx(index_f_Ts_dry)=LEp_mt(index_f_Ts_dry)./(LEp_mt(index_f_Ts_dry)+ (1-alphag(index_f_Ts_dry)).*Rns_dry(index_f_Ts_dry));
        fp(index_f_Ts_dry)=1./(1+apm(index_f_Ts_dry).*exp(bpm(index_f_Ts_dry).*fpx(index_f_Ts_dry)));
        
        %The equation of soil energy balance
        f_Ts_dry(index_f_Ts_dry)=rou*cp*(Ts_dry(index_f_Ts_dry)-air_temp(index_f_Ts_dry))./(ra(index_f_Ts_dry)+rs(index_f_Ts_dry))+fp(index_f_Ts_dry).*f_dry(index_f_Ts_dry).*LEp_pm(index_f_Ts_dry)-(1-alphag(index_f_Ts_dry)).*Rns_dry(index_f_Ts_dry);
        
        %The calculation of derivative by difference method.f'(x)=(f(x+0.001)-f(x))/0.001
        step=0*f_Ts_dry(index_f_Ts_dry)+0.001;
        LEp_mt_step=rou*cp*(a*exp(b*(Ts_dry(index_f_Ts_dry)+step-273.15)./(Ts_dry(index_f_Ts_dry)+step-273.15+c))-ea(index_f_Ts_dry))./(ra(index_f_Ts_dry)+rs(index_f_Ts_dry))./gama;
        
        if interation1==1
            
        else
            [Rnc_dry(index), Rns_dry(index)] = radiation_Campbell(cos_solar_zenith(index),lai(index),incoming_long_radiation(index),incoming_short_radiation(index),Tc_dry(index),(Ts_dry(index) + 0.001));
        end
        
        fpx_step=LEp_mt_step./(LEp_mt_step+(1-alphag(index_f_Ts_dry)).*Rns_dry(index_f_Ts_dry));
        
        fp_step=1./(1+apm(index_f_Ts_dry).*exp(bpm(index_f_Ts_dry).*fpx_step));
        
        LEp_pm_step=(delta(index_f_Ts_dry).*((1-alphag(index_f_Ts_dry)).*Rns_dry(index_f_Ts_dry))+rou*cp*(es(index_f_Ts_dry)-ea(index_f_Ts_dry))./(ra(index_f_Ts_dry)+rs(index_f_Ts_dry)))./delta_gama(index_f_Ts_dry);
        
        f_Ts_dry_step=rou.*cp.*((Ts_dry(index_f_Ts_dry)+0.001)-air_temp(index_f_Ts_dry))./(ra(index_f_Ts_dry)+rs(index_f_Ts_dry))+fp_step.*f_dry(index_f_Ts_dry).*LEp_pm_step-(1-alphag(index_f_Ts_dry)).*Rns_dry(index_f_Ts_dry);
        derivitive_Ts_dry(index_f_Ts_dry)=(f_Ts_dry_step-f_Ts_dry(index_f_Ts_dry))./0.001;
        
        %If LEp_mt<LEp_pm, LEp_mt is replaced by LEp_pm
        if interation1==1  %Rns is independent of Ts in intial calculation
            derivitive_Ts_dry(index_small)=rou*cp./(ra(index_small)+rs(index_small));
        else
            LEp_pm_step=(delta(index_small).*((1-alphag(index_small)).*Rns_dry(index_small))+rou*cp*(es(index_small)-ea(index_small))./(ra(index_small)+rs(index_small)))./delta_gama(index_small);
            
            fpx_step2= LEp_pm_step./( LEp_pm_step+(1-alphag(index_small)).*Rns_dry(index_small));
            fp_step2=1./(1+apm(index_small).*exp(bpm(index_small).*fpx_step2));
            f_Ts_dry_step2=rou*cp*((Ts_dry(index_small)+0.001)-air_temp(index_small))./(ra(index_small)+rs(index_small))+fp_step2.*f_dry(index_small).*LEp_pm_step-(1-alphag(index_small)).*Rns_dry(index_small);
            
            derivitive_Ts_dry(index_small)=(f_Ts_dry_step2-f_Ts_dry(index_small))./0.001;
        end
        
        % The iteration (Newton's method).
        Ts_dry(index_f_Ts_dry)=Ts_dry(index_f_Ts_dry)-f_Ts_dry(index_f_Ts_dry)./derivitive_Ts_dry(index_f_Ts_dry);
        
        %calculate closure error
        ets_dry(index_f_Ts_dry)=a*exp(b*(Ts_dry(index_f_Ts_dry)-273.15)./(Ts_dry(index_f_Ts_dry)-273.15+c)); %es?????? ??hpa
        LEp_mt(index_f_Ts_dry)=rou*cp*(ets_dry(index_f_Ts_dry)-ea(index_f_Ts_dry))./(ra(index_f_Ts_dry)+rs(index_f_Ts_dry))./gama;
        
        if interation1==1
            
        else
            [Rnc_dry(index), Rns_dry(index)] = radiation_Campbell(cos_solar_zenith(index),lai(index),incoming_long_radiation(index),incoming_short_radiation(index),Tc_dry(index),(Ts_dry(index)));
        end
        
        LEp_pm(index_f_Ts_dry)=(delta(index_f_Ts_dry).*((1-alphag(index_f_Ts_dry)).*Rns_dry(index_f_Ts_dry))+rou*cp*(es(index_f_Ts_dry)-ea(index_f_Ts_dry))./(ra(index_f_Ts_dry)+rs(index_f_Ts_dry)))./delta_gama(index_f_Ts_dry);
        
        
        
        
        
        index_small=find(LEp_mt<LEp_pm);
        LEp_mt(index_small)=LEp_pm(index_small);
        
        fpx(index_f_Ts_dry)=LEp_mt(index_f_Ts_dry)./(LEp_mt(index_f_Ts_dry)+ (1-alphag(index_f_Ts_dry)).*Rns_dry(index_f_Ts_dry));
        fp(index_f_Ts_dry)=1./(1+apm(index_f_Ts_dry).*exp(bpm(index_f_Ts_dry).*fpx(index_f_Ts_dry)));
        f_Ts_dry(index_f_Ts_dry)=rou*cp*(Ts_dry(index_f_Ts_dry)-air_temp(index_f_Ts_dry))./(ra(index_f_Ts_dry)+rs(index_f_Ts_dry))+fp(index_f_Ts_dry).*f_dry(index_f_Ts_dry).*LEp_pm(index_f_Ts_dry)-(1-alphag(index_f_Ts_dry)).*Rns_dry(index_f_Ts_dry);
        
        % Re-iterative calculation of points with excessive closure difference
        index_f_Ts_dry=find((temp_plus>0) & abs(f_Ts_dry)>5);
    end
    
    
    LEs_wet(index_wet)=1.26*delta_gama(index_wet).*(Rns_wet(index_wet)-G_wet(index_wet));
    Hs_wet(index_wet)=(Rns_wet(index_wet)-G_wet(index_wet))-LEs_wet(index_wet);
    Ts_wet(index_wet)=Hs_wet(index_wet).*(ra(index_wet)+rs(index_wet))./rou./cp+air_temp(index_wet);
    
    Hc_dry(index)=rou*cp*(Tc_dry(index)-air_temp(index))./(ra(index)+rb(index));
    Hs_dry(index)=rou*cp*(Ts_dry(index)-air_temp(index))./(ra(index)+rs(index));
    H(index)=Hc_dry(index)+Hs_dry(index)+Hc_wet(index)+Hs_wet(index);
    
    LEc_dry(index)=Rnc_dry(index)-Hc_dry(index);
    LEs_dry(index)=(1-alphag(index)).*Rns_dry(index)-Hs_dry(index);
    LE(index)=LEc_dry(index)+LEs_dry(index)+LEc_wet(index)+LEs_wet(index);
    
    Tsurf_new(index)=((1-fwet(index)).*fr(index).*Tc_dry(index).^4+fwet(index).*fr(index).*Tc_wet(index).^4+(1-fwet(index)).*(1-fr(index)).*Ts_dry(index).^4+fwet(index).*(1-fr(index)).*Ts_wet(index).^4).^0.25;
    
end
 
%% calculate the objective function value Z
temp_evaluate=zeros(length(air_temp),1);
temp_evaluate(index)=1;
index_evaluate=find((temp_evaluate>0)); %those points that meet the basic calculation conditions but require forward correction
 
 
% Tsurf(index_evaluate)=((1-fwet(index_evaluate)).*fr(index_evaluate).*Tc_dry(index_evaluate).^4+fwet(index_evaluate).*fr(index_evaluate).*Tc_wet(index_evaluate).^4+(1-fwet(index_evaluate)).*(1-fr(index_evaluate)).*Ts_dry(index_evaluate).^4+fwet(index_evaluate).*(1-fr(index_evaluate)).*Ts_wet(index_evaluate).^4).^0.25;
length(index_evaluate);
 
Hc_dry(index_evaluate)=rou*cp*(Tc_dry(index_evaluate)-air_temp(index_evaluate))./(ra(index_evaluate)+rb(index_evaluate));
Hs_dry(index_evaluate)=rou*cp*(Ts_dry(index_evaluate)-air_temp(index_evaluate))./(ra(index_evaluate)+rs(index_evaluate));
H(index_evaluate)=Hc_dry(index_evaluate)+Hs_dry(index_evaluate);
 
LEc_dry(index_evaluate)=Rnc_dry(index_evaluate)-Hc_dry(index_evaluate);
LEs_dry(index_evaluate)=(1-alphag(index_evaluate)).*Rns_dry(index_evaluate)-Hs_dry(index_evaluate);
LE(index_evaluate)=LEc_dry(index_evaluate)+LEs_dry(index_evaluate);
 
 
 
% %calculate the RMSE between the simulated value and the measured value, the aim of annealing algorithm is to minimize the Z value.
% z=sqrt(sum((Tsurf_new(index_evaluate)-Tsurf_mea(index_evaluate)).^2)/length(index_evaluate));
%  
% %If the number of simulated values that meet the decision condition is too small, then increase the Z value.
% if length(index_evaluate)<0.9*length(index)
%     z=z+10;
% else
%     z=z;
% end
%  
% % if (sum(LEc_dry(index_evaluate))<sum(0.85*delta_gama(index_evaluate).*Rnc_dry(index_evaluate)) | sum(LEc_dry(index_evaluate))>sum(1.65*delta_gama(index_evaluate).*Rnc_dry(index_evaluate)))
% %     z=z+20;
% % else
% %     z=z;
% % end
%  
% % If the input value that meets the original condition is 0, the Z value is 0
% if length(index)==0
%     z=0;
% else
%     z=z;
% end
 
interation_total=interation1.*hc./hc;

end


