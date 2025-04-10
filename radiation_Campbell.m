function [Rnc_dry,Rns_dry,Rn] = radiation_Campbell(cos_solar_zenith,lai,incoming_long_radiation,incoming_short_radiation,Tc_dry,Ts_dry)

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
% 	    Rn : float
% 	        Net radiation (W m-2).
%  
% 	    References
% 	    ----------
% 	    .. [Campbell1998] Campbell, G. S. & Norman, J. M. (1998), An introduction to environmental
% 	        biophysics. Springer, New York
% 	        https://archive.org/details/AnIntroductionToEnvironmentalBiophysics.
% 	    """
%       D. Kool, W.P. Kustas, A. Ben-Gal, N. Agam, Energy partitioning between plant canopy and soil, performance of the two-source energy balance model in a vineyard,
%       Agricultural and Forest Meteorology,Volume 300,2021,https://doi.org/10.1016/j.agrformet.2021.108328.


sigma=5.67e-8;
T_C=Tc_dry;
T_S=Ts_dry;
L_dn=incoming_long_radiation;


emis_c=0.98;  %canopy emissivity
emis_s=0.95;  %soil emissivity

% tan_solar_zenith2=1./(cos_solar_zenith.^2)-1;

% Compute Short Radiation
[difvis,difnir,fvis,fnir]=calc_difuse_ratio(incoming_short_radiation,cos_solar_zenith); % fraction of difuse and PAR/NIR radiation from shortwave irradiance (W m-2, solar zenith angle, atmospheric pressure and precipitable water vapour )
Skyl=difvis.*fvis+difnir.*fnir; % broadband difuse fraction
Sdn_dir=incoming_short_radiation.*(1.0-Skyl);
Sdn_dif=incoming_short_radiation.*Skyl;

%%leaf radiation absorption for shortwave radiation
a_vis=0.85;
a_nir=0.15;

%soil albedo
% alb_s=0.2;
alb_s_vis=0.15;
alb_s_nir=0.25;

% 	    albb : float or array_like
% 	        Beam (black sky) canopy albedo
% 	    albd : float or array_like
% 	        Diffuse (white sky) canopy albedo
% 	    taubt : float or array_like
% 	        Beam (black sky) canopy transmittance
% 	    taudt : float or array_like
% 	        Beam (white sky) canopy transmittance
[albb_vis, albd_vis, taubt_vis, taudt_vis]=calc_spectra_Cambpell(lai, cos_solar_zenith, a_vis, alb_s_vis);
[albb_nir, albd_nir, taubt_nir, taudt_nir]=calc_spectra_Cambpell(lai, cos_solar_zenith, a_nir, alb_s_nir);

    Sn_C = ((1.0 - taubt_vis).* (1.0- albb_vis).* Sdn_dir.*fvis...
            + (1.0 - taubt_nir).* (1.0- albb_nir).* Sdn_dir.*fnir...
            + (1.0 - taudt_vis).* (1.0- albd_vis).* Sdn_dif.*fvis...
            + (1.0 - taudt_nir).* (1.0- albd_nir).* Sdn_dif.*fnir);
            
    Sn_S = (taubt_vis.* (1.0 - alb_s_vis).* Sdn_dir.*fvis...
            + taubt_nir.* (1.0 - alb_s_nir).* Sdn_dir.*fnir...
            + taudt_vis.* (1.0 - alb_s_vis).* Sdn_dif.*fvis...
            + taudt_nir.* (1.0 - alb_s_nir).* Sdn_dif.*fnir);


% Compute Long Radiation


% 	    # calculate long wave emissions from canopy, soil and sky
 
	    L_C = emis_c.*sigma.*T_C.^4;
% 	    L_C[np.isnan(L_C)] = 0
	    L_S = emis_s.*sigma.*T_S.^4;
% 	    L_S[np.isnan(L_S)] = 0
% 	    # Calculate the canopy spectral properties
% 	    _, albl, _, taudl = calc_spectra_Cambpell(lai,np.zeros(emisVeg.shape),1.0 - emisVeg,np.zeros(emisVeg.shape),1.0 - emisGrd,
%         [albb, albl, taubt, taudl]=calc_spectra_Cambpell(lai, cos_solar_zenith, (1.0-emisVeg), cos_solar_zenith, (1.0-emisGrd));
% 	    clear albb taubt albl

% 	    # calculate net longwave radiation divergence of the soil
% 	    L_nS = emisGrd.* taudl.* L_dn + emisGrd.* (1.0 - taudl).* L_C - L_S;
% 	    L_nC = (1 - albl).* (1.0 - taudl).* (L_dn + L_S) - 2.0.* (1.0 - taudl).* L_C;
        
        trans_long=exp(-0.95.*lai);   %Transmission of longwave radiation through the canopy D. Kool et al.,2021, AFM

        Ln_S = trans_long.* L_dn + (1.0 - trans_long).* L_C - L_S;
        Ln_C = (1.0 - trans_long).* (L_dn + L_S - 2.0.* L_C);
        
% 	    # Compute Net Radiation
	    Rns_dry = Sn_S + Ln_S;
	    Rnc_dry = Sn_C + Ln_C;
	    Rn = Rns_dry + Rnc_dry;

    
end