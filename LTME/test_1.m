

clearvars temp
for i = 1:size(lon,1)
    for j = 1:size(lon,1)
            % temp
            temp_t_s=squeeze(temp_sur_clim_g_f(i,j,:)); % pick time series of one point
            temp_t_b=squeeze(temp_bot_clim_g_f(i,j,:)); % pick time series of one point
            
            % salt
            temp_s_s=squeeze(salt_sur_clim_g_f(i,j,:)); % pick time series of one point
            temp_s_b=squeeze(salt_bot_clim_g_f(i,j,:)); % pick time series of one point
            
            % do
            temp_d_s=squeeze(do_sur_clim_g_f(i,j,:)); % pick time series of one point
            temp_d_b=squeeze(do_bot_clim_g_f(i,j,:)); % pick time series of one point
            
            % nh4
            temp_nh_s=squeeze(nh4_sur_clim_g_f(i,j,:)); % pick time series of one point
            temp_nh_b=squeeze(nh4_bot_clim_g_f(i,j,:)); % pick time series of one point
            
            % no3
            temp_no_s=squeeze(no3_sur_clim_g_f(i,j,:)); % pick time series of one point
            temp_no_b=squeeze(no3_bot_clim_g_f(i,j,:)); % pick time series of one point
            
            % chl
            temp_ch_s=squeeze(chl_sur_clim_g_f(i,j,:)); % pick time series of one point
            temp_ch_b=squeeze(chl_bot_clim_g_f(i,j,:)); % pick time series of one point

     
     temp_t_s = [temp_t_s(end); NaN; temp_t_s; NaN; NaN; temp_t_s(2);];
     temp_t_b = [temp_t_b(end); NaN; temp_t_b; NaN; NaN; temp_t_b(2);];
     
     temp_s_s = [temp_s_s(end); NaN; temp_s_s; NaN; NaN; temp_s_s(2);];
     temp_s_b = [temp_s_b(end); NaN; temp_s_b; NaN; NaN; temp_s_b(2);];
     
     temp_d_s = [temp_d_s(end); NaN; temp_d_s; NaN; NaN; temp_d_s(2);];
     temp_d_b = [temp_d_b(end); NaN; temp_d_b; NaN; NaN; temp_d_b(2);];
     
     temp_nh_s = [temp_nh_s(end); NaN; temp_nh_s; NaN; NaN; temp_nh_s(2);];
     temp_nh_b = [temp_nh_b(end); NaN; temp_nh_b; NaN; NaN; temp_nh_b(2);];
     
     temp_no_s = [temp_no_s(end); NaN; temp_no_s; NaN; NaN; temp_no_s(2);];
     temp_no_b = [temp_no_b(end); NaN; temp_no_b; NaN; NaN; temp_no_b(2);];
     
     temp_ch_s = [temp_ch_s(end); NaN; temp_ch_s; NaN; NaN; temp_ch_s(2);];
     temp_ch_b = [temp_ch_b(end); NaN; temp_ch_b; NaN; NaN; temp_ch_b(2);];
     
     % filling NaNs on the time order
     t = 1:length(temp_t_s); %general time
     temp_t_s(isnan(temp_t_s)) = interp1(t(~isnan(temp_t_s)), temp_t_s(~isnan(temp_t_s)), t(isnan(temp_t_s)) ); 
     temp_t_b(isnan(temp_t_b)) = interp1(t(~isnan(temp_t_b)), temp_t_b(~isnan(temp_t_b)), t(isnan(temp_t_b)) ); 
     temp_t_s(1:2)=[]; temp_t_s(end-1:end)=[]; % 11:12, 1:2 mth removing
     temp_t_b(1:2)=[]; temp_t_b(end-1:end)=[]; % 11:12, 1:2 mth removing
     
     temp_s_s(isnan(temp_s_s)) = interp1(t(~isnan(temp_s_s)), temp_s_s(~isnan(temp_s_s)), t(isnan(temp_s_s)) ); 
     temp_s_b(isnan(temp_s_b)) = interp1(t(~isnan(temp_s_b)), temp_s_b(~isnan(temp_s_b)), t(isnan(temp_s_b)) );
     temp_s_s(1:2)=[]; temp_s_s(end-1:end)=[]; % 11:12, 1:2 mth removing
     temp_s_b(1:2)=[]; temp_s_b(end-1:end)=[]; % 11:12, 1:2 mth removing
     
     temp_d_s(isnan(temp_d_s)) = interp1(t(~isnan(temp_d_s)), temp_d_s(~isnan(temp_d_s)), t(isnan(temp_d_s)) ); 
     temp_d_b(isnan(temp_d_b)) = interp1(t(~isnan(temp_d_b)), temp_d_b(~isnan(temp_d_b)), t(isnan(temp_d_b)) );
     temp_d_s(1:2)=[]; temp_d_s(end-1:end)=[]; % 11:12, 1:2 mth removing
     temp_d_b(1:2)=[]; temp_d_b(end-1:end)=[]; % 11:12, 1:2 mth removing
     
     temp_nh_s(isnan(temp_nh_s)) = interp1(t(~isnan(temp_nh_s)), temp_nh_s(~isnan(temp_nh_s)), t(isnan(temp_nh_s)) ); 
     temp_nh_b(isnan(temp_nh_b)) = interp1(t(~isnan(temp_nh_b)), temp_nh_b(~isnan(temp_nh_b)), t(isnan(temp_nh_b)) );
     temp_nh_s(1:2)=[]; temp_nh_s(end-1:end)=[]; % 11:12, 1:2 mth removing
     temp_nh_b(1:2)=[]; temp_nh_b(end-1:end)=[]; % 11:12, 1:2 mth removing
     
     temp_no_s(isnan(temp_no_s)) = interp1(t(~isnan(temp_no_s)), temp_no_s(~isnan(temp_no_s)), t(isnan(temp_no_s)) ); 
     temp_no_b(isnan(temp_no_b)) = interp1(t(~isnan(temp_no_b)), temp_no_b(~isnan(temp_no_b)), t(isnan(temp_no_b)) );
     temp_no_s(1:2)=[]; temp_no_s(end-1:end)=[]; % 11:12, 1:2 mth removing
     temp_no_b(1:2)=[]; temp_no_b(end-1:end)=[]; % 11:12, 1:2 mth removing
     
     temp_ch_s(isnan(temp_ch_s)) = interp1(t(~isnan(temp_ch_s)), temp_ch_s(~isnan(temp_ch_s)), t(isnan(temp_ch_s)) ); 
     temp_ch_b(isnan(temp_ch_b)) = interp1(t(~isnan(temp_ch_b)), temp_ch_b(~isnan(temp_ch_b)), t(isnan(temp_ch_b)) );
     temp_ch_s(1:2)=[]; temp_ch_s(end-1:end)=[]; % 11:12, 1:2 mth removing
     temp_ch_b(1:2)=[]; temp_ch_b(end-1:end)=[]; % 11:12, 1:2 mth removing
     
      % interp on the layer order  

     prepretemp = [NaN(20,12)];
     prepretemp(20,1:12) = temp_t_s'; prepretemp(1,1:12) = temp_t_b';
     
     prepresalt = [NaN(20,12)];
     prepresalt(20,1:12) = temp_s_s'; prepresalt(1,1:12) = temp_s_b';
     
     prepredo = [NaN(20,12)];
     prepredo(20,1:12) = temp_d_s'; prepredo(1,1:12) = temp_d_b';
     
     preprenh4 = [NaN(20,12)];
     preprenh4(20,1:12) = temp_nh_s'; preprenh4(1,1:12) = temp_nh_b';
     
     prepreno3 = [NaN(20,12)];
     prepreno3(20,1:12) = temp_no_s'; prepreno3(1,1:12) = temp_no_b';
     
     preprech = [NaN(20,12)];
     preprech(20,1:12) = temp_ch_s'; preprech(1,1:12) = temp_ch_b';
     
          for k = 1:12 %time
              temp_t = squeeze(prepretemp(:,k));
              temp_s = squeeze(prepresalt(:,k));
              temp_do =  squeeze(prepredo(:,k));
              temp_nh4 =  squeeze(prepreno4(:,k));
              temp_no3 =  squeeze(prepreno3(:,k));
              temp_ch =  squeeze(preprech(:,k));
              
              g_t = 1:length(temp_t_s); %general layer num.
              temp_t(isnan(temp_t))=interp1(g_t(~isnan(temp_t)), temp_t(~isnan(temp_t)), g_t(isnan(temp_t)) );
              temp_s(isnan(temp_s))=interp1(g_t(~isnan(temp_s)), temp_s(~isnan(temp_s)), g_t(isnan(temp_s)) );
              temp_do(isnan(temp_do))=interp1(g_t(~isnan(temp_do)), temp_do(~isnan(temp_do)), g_t(isnan(temp_do)) );
              temp_nh4(isnan(temp_nh4))=interp1(g_t(~isnan(temp_nh4)), temp_nh4(~isnan(temp_nh4)), g_t(isnan(temp_nh4)) );
              temp_no3(isnan(temp_no3))=interp1(g_t(~isnan(temp_no3)), temp_no3(~isnan(temp_no3)), g_t(isnan(temp_no3)) );
              temp_ch(isnan(temp_ch))=interp1(g_t(~isnan(temp_ch)), temp_ch(~isnan(temp_ch)), g_t(isnan(temp_ch)) );
              
              pretemp(:,k) = temp_t;
              presalt(:,k) = temp_s;
              predo(:,k) = temp_do;
              preno4(:,k) = temp_nh4;
              preno3(:,k) = temp_no3;
              prech(:,k) = temp_ch;
              clervars temp_t temp_s temp_do temp_nh4 temp_no3 temp_ch g_t
          end
               
          % completet
          temperature(i,j,:,:) = pretemp;
          salt(i,j,:,:) = presalt;
          do(i,j,:,:) = predo;
          nh4(i,j,:,:) = preno4;
          no3(i,j,:,:) = preno3;
          chla(i,j,:,:) = prech;
          
          clearvars temp_* pre*
    end
end