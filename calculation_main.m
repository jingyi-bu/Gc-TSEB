clear
load Meteorological_data.mat
load RS_data_365.mat


datetime=readtable('day.xlsx');
co2_s=readtable('co2.xlsx');
Ca=co2_s.co2;

g1_filled=geotiffread('./parameters/filled_g1.tif');
etr_a_filled=geotiffread('./parameters/filled_etr_a.tif');
apm_filled=geotiffread('./parameters/filled_apm.tif');
bpm_filled=geotiffread('./parameters/filled_bpm.tif');
alphag_filled=geotiffread('./parameters/filled_alphag.tif');
omi_filled=geotiffread('./parameters/filled_omi.tif');


variables = who;

for idx = 1:length(variables)
    varName = variables{idx};
    varData = eval(varName);
    

    if isnumeric(varData)

        varData(varData == 32767 | varData == -9999 | varData == -32767) = NaN;
        assignin('base', varName, varData);
    end
end



ex=stem_005(:, :, 1);
[dim1, dim2, dim3] = size(ex);

LE_sim_total=zeros(dim1, dim2, 8760)-9999;
LE_sim_c=zeros(dim1, dim2, 8760)-9999;

for ii=148:dim1
    ii
    
    for jj=1:dim2
        
% for ii=73     %US-Var
%     for jj=50

% for ii=207    %US-Wkg
%     for jj=271
        
% for ii=205    %US-SRM 
%     for jj=252       
        
        if ~isnan(ex(ii,jj))

            lon_1 = long(ii,jj);
            lat_1 = latitude(ii,jj);


            if lon_1 < -112.5
            J=datetime.day_UTC_8;
            date_temp=datetime.datetime_UTC_8;
            else
            J=datetime.day_UTC_7;
            date_temp=datetime.datetime_UTC_7;
            end

            

            % lai
            lai = zeros(365*24, 1);
            for k = 1:365
                lai((k-1)*24+1:k*24) = repmat(lai_005(ii, jj, k), [24, 1]);
            end

            % stem
            stem = zeros(365*24, 1);
            for k = 1:365
                stem((k-1)*24+1:k*24) = repmat(stem_005(ii, jj, k), [24, 1]);
            end

            lai_plus_stem=lai+stem;

            % albedo
            albedo = zeros(365*24, 1);
            for k = 1:365
                albedo((k-1)*24+1:k*24) = repmat(albedo_005(ii, jj, k), [24, 1]);
            end        


            air_temp=squeeze(t2m(ii, jj, :));
            air_humi=squeeze(rh(ii, jj, :))*0.01;  
            air_pres=squeeze(sp(ii, jj, :))*0.01;  
            wind=squeeze(wind005(ii, jj, :)); 
            incoming_short_radiation=squeeze(ssrd(ii, jj, :)); 
            incoming_long_radiation=squeeze(strd(ii, jj, :)); 
            Da=squeeze(vpd(ii, jj, :))./10;  %kpa  

            % sif
            sif = zeros(365*24, 1);
            for k = 1:365
                 sif((k-1)*24+1:k*24) = repmat(CSIF(ii, jj, k), [24, 1]);
            end            

            SWC=squeeze(sw(ii, jj, :)); 

            % f_dry
            f_dry = zeros(365*24, 1);
            for k = 1:365
                 f_dry((k-1)*24+1:k*24) = repmat(f_drying(ii, jj, k), [24, 1]);
            end  


            hc=zeros(365*24, 1)+hc_005(ii, jj);

            % d0
            % Further testing is required for trees with excessive height,
            % and a average value was taken here.
            % d0=hc*0.65;
            d0=0.33500000834465+0*hc;   

            % zm0
            zm0 = zeros(365*24, 1);
            for k = 1:365
                 zm0((k-1)*24+1:k*24) = repmat(zm0_005(ii, jj, k), [24, 1]);
            end          



            zh0=zm0./7;      



            zm=10+0*hc; 
            zh=2+0*hc; 
            lontitude=lon_1*(pi./180)+0*hc;
            lat=lat_1*(pi./180)+0*hc; 

            c4=zeros(365*24, 1)+c4_005(ii, jj)*0.01;


            if lucc(ii,jj) == 6 || lucc(ii,jj) ==7
                Qf=0.261;
                Qw=0.059;
            elseif lucc(ii,jj) == 8 || lucc(ii,jj) == 9
                Qf=0.538;
                Qw=0.07;
            elseif lucc(ii,jj) == 10
                Qf=0.175;
                Qw=0.016;
            end


            for kk=1:length(date_temp)

                test=num2str(date_temp(kk));
                time(kk)=str2num(test(9:10)); %time?hour

            end

            g1=g1_filled(ii,jj);
            etr_a=etr_a_filled(ii,jj);
            apm=apm_filled(ii,jj);
            bpm=bpm_filled(ii,jj);
            alphag=alphag_filled(ii,jj);
            omi=omi_filled(ii,jj);

            [LE_sim1,H_sim1,Rn,G_dry,LEc_dry,LEs_dry,Qc,Tc_dry]= Gc_CC_TSEB_flux_Calibration_ca_sifetr_SWC(g1,etr_a,apm,bpm,alphag,omi);
            
            virtual_indices = imag(LE_sim1) ~= 0;
            LE_sim1(virtual_indices) = -9999;
            LE_sim_total(ii, jj, :)=LE_sim1;

            virtual_indices = imag(LEc_dry) ~= 0;
            LEc_dry(virtual_indices) = -9999;          
            LE_sim_c(ii, jj, :)=LEc_dry;
        
        end
    end
end

save('LE_LEc_results_ori_148_215.mat', 'LE_sim_total', 'LE_sim_c', '-v7.3');

% Replace all imaginary elements in the array
% virtual_indices = imag(LE_sim_total) ~= 0;
% LE_sim_total(virtual_indices) = -9999;
% 
% nan_indices = isnan(LE_sim_total);
% LE_sim_total(nan_indices) = -9999;
% 
% virtual_indices = imag(LE_sim_c) ~= 0;
% LE_sim_c(virtual_indices) = -9999;
% 
% nan_indices = isnan(LE_sim_c);
% LE_sim_c(nan_indices) = -9999;
% 
% save('LE_LEc_results.mat', 'LE_sim_total', 'LE_sim_c', '-v7.3');