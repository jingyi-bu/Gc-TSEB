function [Rnc_dry,Rns_dry] = radiation_Kool(cos_solar_zenith,lai,incoming_long_radiation,incoming_short_radiation,Tc_dry,Ts_dry)
% Parameters
% 	    ----------
% 	    lai: float
%            leaf area index (-)
%       incoming_long_radiation : float
% 	        Incoming longwave radiation (W m-2).
%       incoming_short_radiation : float
% 	        Incoming shortwave radiation (W m-2).
% 	    cos_solar_zenith : float
% 	        Cosine of the Solar Zenith Angle(-)
%       Tc_dry: float
%           Canopy temperature(K)
%       Ts_dry: float
%           Soil temperature(K)
%  	
% 	    Returns
% 	    -------
% 	    Rnc_dry : float
% 	        Canopy net radiation (W m-2).
% 	    Rns_dry : float
% 	        Soil net radiation (W m-2).
%  
% 	    References
% 	    ----------
%       D. Kool, W.P. Kustas, A. Ben-Gal, N. Agam, Energy partitioning between plant canopy and soil, performance of the two-source energy balance model in a vineyard,
%       Agricultural and Forest Meteorology,Volume 300,2021,https://doi.org/10.1016/j.agrformet.2021.108328.

sigma=5.67e-8;
emis_c=0.98;  %canopy emissivity
emis_s=0.95;  %soil emissivity

tan_solar_zenith2=1./(cos_solar_zenith.^2)-1;   %Tangent of solar zenith

%transmittance of direct radiation
kb=sqrt(1+tan_solar_zenith2)./(1+1.744.*(1+1.182).^(-0.773));    %Eq.(A.3)

%The extinction coefficient for diffuse radiation
if lai<0.5
    kd=0.9;
elseif lai>=0.5&lai<2
    kd=0.8;
else 
    kd=0.7;
end

%leaf radiation absorption for direct (subscript b) and diffuse (subscript d) shortwave radiation
a_b_par=0.85;
a_d_par=0.85;
a_b_nir=0.15;
a_d_nir=0.15;


%plant canopy reflectance
ref_c_b_par=2.*kb.*(1-sqrt(a_b_par))./(1+kb)./(1+sqrt(a_b_par));     %Eq.(A.4a) 
ref_c_b_nir=2.*kb.*(1-sqrt(a_b_nir))./(1+kb)./(1+sqrt(a_b_nir));     %Eq.(A.4a)
ref_c_d_par=2.*kd.*(1-sqrt(a_d_par))./(1+kd)./(1+sqrt(a_d_par));     %Eq.(A.4b)
ref_c_d_nir=2.*kd.*(1-sqrt(a_d_nir))./(1+kd)./(1+sqrt(a_d_nir));     %Eq.(A.4b)

% soil albedo
alb_s=0.2;
alb_s_par=0.15;
alb_s_nir=0.25;

%canopy albedo
alb_c_b_par_1=ref_c_b_par+(ref_c_b_par-alb_s_par)./(ref_c_b_par.*alb_s_par-1).*exp(-2.*sqrt(a_b_par).*kb.*lai);
alb_c_b_par=alb_c_b_par_1./(1+alb_c_b_par_1);    %Eq.(A.5a) 

alb_c_b_nir_1=ref_c_b_nir+(ref_c_b_nir-alb_s_nir)./(ref_c_b_nir.*alb_s_nir-1).*exp(-2.*sqrt(a_b_nir).*kb.*lai);
alb_c_b_nir=alb_c_b_nir_1./(1+alb_c_b_nir_1);    %Eq.(A.5a) 

alb_c_d_par_1=ref_c_d_par+(ref_c_d_par-alb_s_par)./(ref_c_d_par.*alb_s_par-1).*exp(-2.*sqrt(a_d_par).*kd.*lai);
alb_c_d_par=alb_c_d_par_1./(1+alb_c_d_par_1);   %Eq.(A.5b) 

alb_c_d_nir_1=ref_c_d_nir+(ref_c_d_nir-alb_s_nir)./(ref_c_d_nir.*alb_s_nir-1).*exp(-2.*sqrt(a_d_nir).*kd.*lai);
alb_c_d_nir=alb_c_d_nir_1./(1+alb_c_d_nir_1);   %Eq.(A.5b) 

% canopy transmission
trans_c_b_par=(ref_c_b_par.^2-1).*exp(-2.*sqrt(a_b_par).*kb.*lai)./((ref_c_b_par.*alb_s_par-1)+ref_c_b_par.*(ref_c_b_par-alb_s_par).*exp(-2.*sqrt(a_b_par).*kb.*lai));  %Eq.(A.6a)
trans_c_b_nir=(ref_c_b_nir.^2-1).*exp(-2.*sqrt(a_b_nir).*kb.*lai)./((ref_c_b_nir.*alb_s_nir-1)+ref_c_b_nir.*(ref_c_b_nir-alb_s_nir).*exp(-2.*sqrt(a_b_nir).*kb.*lai));  %Eq.(A.6a)

trans_c_d_par=(ref_c_d_par.^2-1).*exp(-2.*sqrt(a_d_par).*kd.*lai)./((ref_c_d_par.*alb_s_par-1)+ref_c_d_par.*(ref_c_d_par-alb_s_par).*exp(-2.*sqrt(a_d_par).*kd.*lai));  %Eq.(A.6b)
trans_c_d_nir=(ref_c_d_nir.^2-1).*exp(-2.*sqrt(a_d_nir).*kd.*lai)./((ref_c_d_nir.*alb_s_nir-1)+ref_c_d_nir.*(ref_c_d_nir-alb_s_nir).*exp(-2.*sqrt(a_d_nir).*kd.*lai));  %Eq.(A.6b)


alb_c=0.5.*(0.2.*alb_c_d_par+0.8.*alb_c_b_par)+0.5.*(0.1.*alb_c_d_nir+0.9.*alb_c_b_par);  %Eq.(A.7)

trans_short=0.5.*(0.2.*trans_c_d_par+0.8.*trans_c_b_par)+0.5.*(0.1.*trans_c_d_nir+0.9.*trans_c_b_par);  %Eq.(A.8)


% Transmission of longwave radiation through the canopy
trans_long=exp(-0.95.*lai);    %Eq.(A.9)


Rnc_dry=(1-trans_long).*(incoming_long_radiation+emis_s.*sigma.*Ts_dry.^4-2.*emis_c.*sigma.*Tc_dry.^4)...
    +(1-trans_short).*(1-alb_c).*incoming_short_radiation;    %Eq.(A.10)

Rns_dry=(trans_long.*incoming_long_radiation+(1-trans_long).*emis_c.*sigma.*Tc_dry.^4-emis_s.*sigma.*Ts_dry.^4)...
    +trans_short.*(1-alb_s).*incoming_short_radiation;        %Eq.(A.11)
    
end