function[Rdirvis, Rdifvis, Rdirnir, Rdifnir] = calc_potential_irradiance_weiss(cos_solar_zenith)

% cos_solar_zenith=[0.1  0.5 0.9];
        
        press=1013.25;
        SOLAR_CONSTANT=1320;
        fnir_ini=0.5455;
%     ''' Estimates the potential visible and NIR irradiance at the surface
%     Parameters
%     ----------
%     cos_solar_zenith: the cosine of solar Zenith Angle
% 
%     Returns
%     -------
%     Rdirvis : float
%         Potential direct visible irradiance at the surface (W m-2)
%     Rdifvis : float
%         Potential diffuse visible irradiance at the surface (W m-2)
%     Rdirnir : float
%         Potential direct NIR irradiance at the surface (W m-2)
%     Rdifnir : float
%         Potential diffuse NIR irradiance at the surface (W m-2)
%     based on Weiss & Normat 1985, following same strategy in Cupid's RADIN4 subroutine.
%     '''

%     # Convert input scalars to numpy arrays
%     sza, press = map(np.asarray, (sza, press))
% 
%     # Set defaout ouput values
%     Rdirvis, Rdifvis, Rdirnir, Rdifnir, w = [
%         np.zeros(sza.shape) for i in range(5)]

    coszen = cos_solar_zenith;
%     # Calculate potential (clear-sky) visible and NIR solar components
%     # Weiss & Norman 1985
%     # Correct for curvature of atmos in airmas (Kasten and Young,1989)
%     i = sza < 90

    i=find(coszen==0);
    coszen(i)=coszen(i)+0.1;
    clear i
    airmas = 1.0 ./ coszen;
%     # Visible PAR/NIR direct beam radiation
    Sco_vis = SOLAR_CONSTANT .* (1.0 - fnir_ini);
    Sco_nir = SOLAR_CONSTANT .* fnir_ini;
%     # Directional trasnmissivity
%     # Calculate water vapour absorbance (Wang et al 1976)
%     # A=10**(-1.195+.4459*np.log10(1)-.0345*np.log10(1)**2)
%     # opticalDepth=np.log(10.)*A
%     # T=np.exp(-opticalDepth/coszen)
%     # Asssume that most absortion of WV is at the NIR
    Rdirvis = (Sco_vis.* exp(-0.185 .* (press/ 1313.25) .* airmas)) .* coszen;
%                   # Modified Eq1 assuming water vapor absorption
%     # Rdirvis=(Sco_vis*exp(-.185*(press/1313.25)*airmas))*coszen
%     # #Eq. 1
    i=find(Rdirvis<0);
    Rdirvis(i)=0;
    clear i
    
%     # Potential diffuse radiation
%     # Eq 3                                      #Eq. 3
    Rdifvis = 0.4 .* (Sco_vis .* coszen - Rdirvis);
%     Rdifvis = np.maximum(0, Rdifvis)
    
    i=find(Rdifvis<0);
    Rdifvis(i)=0;
    clear i
    
%     # Same for NIR
%     # w=SOLAR_CONSTANT*(1.0-T)
    w = SOLAR_CONSTANT .* 10.^(-1.195 + 0.4459.* log10(coszen) - 0.0345.*(log10(coszen)).^2);  % Eq. .6
    Rdirnir = (Sco_nir .* exp(-0.06 .* (press./ 1313.25).* airmas) - w).* coszen;  % Eq. 4
%     Rdirnir = np.maximum(0, Rdirnir)

    i=find(Rdirnir<0);
    Rdirnir(i)=0;
    clear i
%     # Potential diffuse radiation
    Rdifnir = 0.6.* (Sco_nir .* coszen - Rdirvis - w);  % Eq. 5
%     Rdifnir = np.maximum(0, Rdifnir)
    i=find(Rdifnir<0);
    Rdifnir(i)=0;
    

end