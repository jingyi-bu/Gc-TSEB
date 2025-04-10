function[difvis,difnir,fvis,fnir]=calc_difuse_ratio(incoming_short_radiation,cos_solar_zenith)
 
% def calc_difuse_ratio(S_dn, sza, press=1013.25, SOLAR_CONSTANT=1320):
% 	    Fraction of difuse shortwave radiation.
	
% 	    Partitions the incoming solar radiation into PAR and non-PR and
% 	    diffuse and direct beam component of the solar spectrum.
% 	
% 	    Parameters
% 	    ----------
% 	    S_dn : float
% 	        Incoming shortwave radiation (W m-2).
% 	    cos_solar_zenith : float
% 	        Cosine of the Solar Zenith Angle.
% 	    Wv : float, optional
% 	        Total column precipitable water vapour (g cm-2), default 1 g cm-2.
% 	    press : float, optional
% 	        atmospheric pressure (mb), default at sea level (1013mb).
% 	
% 	    Returns
% 	    -------
% 	    difvis : float
% 	        diffuse fraction in the visible region.
% 	    difnir : float
% 	        diffuse fraction in the NIR region.
% 	    fvis : float
% 	        fration of total visible radiation.
% 	    fnir : float
% 	        fraction of total NIR radiation.
% 	
% 	    References
% 	    ----------
% 	    .. [Weiss1985] Weiss and Norman (1985) Partitioning solar radiation into direct and diffuse,
% 	        visible and near-infrared components, Agricultural and Forest Meteorology,
% 	        Volume 34, Issue 2, Pages 205-213,
% 	        http://dx.doi.org/10.1016/0168-1923(85)90020-6.
 	    

	
      
       Wv=1;
       press=1013.25;
       SOLAR_CONSTANT=1320;
       S_dn=incoming_short_radiation; 
        
	    %Convert input scalars to numpy arrays
% 	    S_dn, sza, press = map(np.asarray, (S_dn, sza, press))
% 	    difvis, difnir, fvis, fnir = [np.zeros(S_dn.shape) for i in range(4)]
% 	    fvis = fvis + 0.6
% 	    fnir = fnir + 0.4
	

% 	    Calculate potential (clear-sky) visible and NIR solar components
% 	    Weiss & Norman 1985
	    [Rdirvis, Rdifvis, Rdirnir, Rdifnir] = calc_potential_irradiance_weiss(cos_solar_zenith);
	

% 	    # Potential total solar radiation
	    potvis = Rdirvis + Rdifvis;
        i=find(potvis <= 0);
        potvis(i)=1e-6;
         clear i
% 	    potvis[potvis <= 0] = 1e-6
	    potnir = Rdirnir + Rdifnir;
% 	    potnir[potnir <= 0] = 1e-6
       i=find(potnir <= 0);
        potnir(i)=1e-6;
         clear i
	    fclear = S_dn ./ (potvis + potnir);
        i=find(fclear >= 1);
        fclear(i)=1;
         clear i
% 	    fclear = np.minimum(1.0, fclear)
	

% 	    # Partition S_dn into VIS and NIR
	    fvis = potvis ./ (potvis + potnir);  %# Eq. 7
	    fnir = potnir ./ (potvis + potnir);  %# Eq. 8
% 	    fvis = np.clip(fvis, 0.0, 1.0)
        i=find(fvis >= 1);
        fvis(i)=1;
         clear i
        i=find(fvis <= 0);
        fvis(i)=0;
         clear i
	    fnir = 1.0 - fvis;
	

% 	    # Estimate direct beam and diffuse fractions in VIS and NIR wavebands
	    ratiox = fclear;
	    ratiox(fclear > 0.9) = 0.9;
	    dirvis = (Rdirvis./ potvis).*(1- ((0.9 - ratiox)./0.7).^0.6667);  %# Eq. 11
	    ratiox = fclear;
	    ratiox(fclear > 0.88) = 0.88;
	    dirnir = (Rdirnir./ potnir).*(1- ((0.88 - ratiox)./ 0.68).^0.6667);  %# Eq. 12
	

% 	    dirvis = np.clip(dirvis, 0.0, 1.0)
        i=find(dirvis >= 1);
        dirvis(i)=1;
         clear i
        i=find(dirvis <= 0);
        dirvis(i)=0;
         clear i
% 	    dirnir = np.clip(dirnir, 0.0, 1.0)
        i=find(dirnir >= 1);
        dirnir(i)=1;
         clear i
        i=find(dirnir <= 0);
        dirnir(i)=0;
         clear i
	    difvis = 1.0 - dirvis;
	    difnir = 1.0 - dirnir;
	


        
end