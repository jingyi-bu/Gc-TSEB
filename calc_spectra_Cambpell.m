function[albb, albd, taubt, taudt]=calc_spectra_Cambpell(lai, cos_solar_zenith, amean, rho_soil)



% 	    """ Canopy spectra
% 	
% 	    Estimate canopy spectral using the [Campbell1998]_
% 	    Radiative Transfer Model
% 	
% 	    Parameters
% 	    ----------
% 	    lai : float
% 	        Effective Leaf (Plant) Area Index.
% 	    sza : float
% 	        Sun Zenith Angle (degrees).
%       amean: float, or array_like
%           Leaf bihemispherical aborprtivity
% 	    rho_soil : float
% 	        Soil bihemispherical reflectance
% 	    x_lad : float,  optional
% 	        x parameter for the ellipsoildal Leaf Angle Distribution function of
% 	        Campbell 1988 [default=1, spherical LIDF].
% 	    lai_eff : float or None, optional
% 	        if set, its value is the directional effective LAI
% 	        to be used in the beam radiation, if set to None we assume homogeneous canopies.
% 	
% 	    Returns
% 	    -------
% 	    albb : float or array_like
% 	        Beam (black sky) canopy albedo
% 	    albd : float or array_like
% 	        Diffuse (white sky) canopy albedo
% 	    taubt : float or array_like
% 	        Beam (black sky) canopy transmittance
% 	    taudt : float or array_like
% 	        Beam (white sky) canopy transmittance
% 	
% 	    References
% 	    ----------
% 	    .. [Campbell1998] Campbell, G. S. & Norman, J. M. (1998), An introduction to environmental
% 	        biophysics. Springer, New York
% 	        https://archive.org/details/AnIntroductionToEnvironmentalBiophysics.
% 	    """
%       D. Kool, W.P. Kustas, A. Ben-Gal, N. Agam, Energy partitioning between plant canopy and soil, performance of the two-source energy balance model in a vineyard,
%       Agricultural and Forest Meteorology,Volume 300,2021,https://doi.org/10.1016/j.agrformet.2021.108328.
	
        x_lad=1;
% 	    # calculate aborprtivity
	    amean_sqrt = sqrt(amean);
	    clear   amean 
	

% 	    # Calculate canopy beam extinction coefficient
% 	    # Modification to include other LADs

	        lai_eff = lai;

	

% 	    # D I F F U S E   C O M P O N E N T S
% 	    # Integrate to get the diffuse transmitance
% 	    taud = calc_taud(x_lad, lai)
% 	  taud = 0.*lai;
%     
%     for angle=0:5:90     
%         
%         angle=angle.*pi./180;
%         akd = sqrt(x_lad.^2 + tan(angle).^2)./(x_lad + 1.774.*(x_lad + 1.182).^(-0.733));
%         taub = exp(-akd.* lai);
%         taud = taud+taub.* cos(angle).* sin(angle).*(5.*pi./180);
%         taud= 2.0.* taud;
%     end
%     
%   
%     clear akd
    
% 	    # Diffuse light canopy reflection coefficients  for a deep canopy
%       The extinction coefficient for diffuse radiation (Kd) is independent of solar zenith, 
%       and has been approximated to equal 0.9 for LAI <0.5, 0.8 for 0.5
%       <LAI <2, and 0.7 for LAI >2. D. Kool et al.,2021, AFM
% 	    akd = -log(taud)./ lai;   %???
    
        if lai<0.5
            akd=0.9;
        elseif lai>=0.5&lai<2
             akd=0.8;
        else 
            akd=0.7;
        end
        
        
	    rcpy= (1.0 - amean_sqrt)./ (1.0 + amean_sqrt);  %# Eq 15.7
	    rdcpy = 2.0.* akd.* rcpy./ (akd + 1.0);  %# Eq 15.8
	

% 	    # Diffuse canopy transmission and albedo coeff for a generic canopy (visible)
	    expfac = amean_sqrt.* akd.* lai;
% 	    clear akd
	    neg_exp = exp(-expfac);
        d_neg_exp =exp(-2.0.* expfac);
	    xnum = (rdcpy.* rdcpy - 1.0).* neg_exp;
	    xden = (rdcpy.* rho_soil - 1.0) + rdcpy.* (rdcpy - rho_soil).* d_neg_exp;
	    taudt = xnum./ xden;  %# Eq 15.11
% 	    clear xnum  xden
	    fact = ((rdcpy - rho_soil)./ (rdcpy.* rho_soil - 1.0)).* d_neg_exp;
	    albd = (rdcpy + fact)./ (1.0 + rdcpy.* fact);  %# Eq 15.9
% 	    clear rdcpy fact 
	

% 	    # B E A M   C O M P O N E N T S
% 	    # Direct beam extinction coeff (spher. LAD)
%       Calculates the beam extinction coefficient based on [Campbell1998]_ ellipsoidal  %# Eq. 15.4
        tan_solar_zenith2=1./(cos_solar_zenith.^2)-1;
        akb = (sqrt(x_lad.^2 + tan_solar_zenith2)./(x_lad + 1.774.*(x_lad + 1.182).^(-0.733)));   
	

% 	    # Direct beam canopy reflection coefficients for a deep canopy
	    rbcpy = 2.0.* akb.* rcpy./ (akb + 1.0);  %# Eq 15.8
	    clear rcpy x_lad 
% 	    # Beam canopy transmission and albedo coeff for a generic canopy (visible)
	    expfac = amean_sqrt.* akb.* lai_eff;
	    neg_exp = exp(-expfac);
        d_neg_exp = exp(-2.0.* expfac);
	    clear amean_sqrt akb  lai_eff 
	    xnum = (rbcpy.* rbcpy - 1.0).* neg_exp;
	    xden = (rbcpy.* rho_soil - 1.0) + rbcpy.* (rbcpy - rho_soil).* d_neg_exp;
	    taubt = xnum./ xden;  %# Eq 15.11
	    clear xnum  xden 
	    fact = ((rbcpy - rho_soil)./ (rbcpy.* rho_soil - 1.0)).* d_neg_exp;
	    clear expfac
	    albb = (rbcpy + fact)./ (1.0 + rbcpy.* fact);  %# Eq 15.9
	    clear rbcpy fact 
	
    i=find(isnan(taubt));
    taubt(i)=1;
    clear i

    
    i=find(isnan(taudt));
    taudt(i)=1;
    clear i
 
    
    i=find(isnan(albb));
    albb(i) = rho_soil;
    clear i
    
    
    i= isnan(albd);
    albd(i) = rho_soil;
    clear i
    
	

% 	    return albb, albd, taubt, taudt

        end