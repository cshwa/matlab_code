close all; clear all; close all;

% from : plot_KOEM_st_seasonal_climate_3regime_v5_spatial_sigma_extract_10y
% koem_1sig=load(['koem_climate_10yr_v5_1sig_gy_each_st.mat']); % full (9point) data (extract_sigma_for_each_st)
% koem_2sig=load(['koem_climate_10yr_v5_2sig_gy_each_st.mat']); % full (9point) data (extract_sigma_for_each_st)
% koem_3sig=load(['koem_climate_10yr_v5_3sig_gy_each_st.mat']); % full (9point) data (extract_sigma_for_each_st)

load model_daily_9point_gy_1997to2010.mat % from plot_KOEM_st_GY_model_case_compare_po4_sewer_det_10yr_clim_v3_10yrclim_1st_vs_2nd.m

%% in bay only left on model results
% % station name tag
% % (for preventing the miss picking using tag from 1 col's name on excel) 
[res L]=sort([1,5,4,3,2,6,8,9,7]);
% % sorting confirm
nn=['±¤¾çÇ×'; '±¤¾ç4'; '±¤¾ç3'; '±¤¾ç2'; '±¤¾ç1'; '±¤¾ç5'; '¿©¼ö2'; '¿©¼ö3'; '¿©¼ö1']; % order of origin koem
nn(L,:) % order of yoonjas koem
% sp_gy=c1997.sp_gy;
% name_tag=c1997.name_tag;
% name_tag{sp_gy(I)} % means nn(I,:)
% 
% sp_gy(I(5:9))
% sp_gy(I(1))
% 
% yoonja_indx= I;


%% row : st. point, col : month 
%  eachst_std_gy_no3_ft_1(jst,i)
%  eachst_std_gy_no3_b_ft_1(jst,i)
%  eachst_clim_gy_no3_ft_1(jst,i)
%  eachst_clim_gy_no3_b_ft_1(jst,i)

[res I]=sort([1,5,4,3,2,6,8,9,7]);
% sorting confirm
nn=['±¤¾çÇ×'; '±¤¾ç4'; '±¤¾ç3'; '±¤¾ç2'; '±¤¾ç1'; '±¤¾ç5'; '¿©¼ö2'; '¿©¼ö3'; '¿©¼ö1'];
nn(I,:)
% sp_gy(I(j))

P_MW = 30.973762;
N_MW = 14.006720;

% make eom_d
k=0
for i = 1997:2000
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(1997,j); % 365 days only (model output dont have 366)
    end
end

% make 1980~present
k=0
for i = 2001:2010
            k=k+1;
    for j = 1:12
        eom_d_2nd(k,j) = eomday(1997,j); % 365 days only (model output dont have 366)
    end
end

clearvars ref_*
k=0; m=0;
for i = 1:length(1997:2000)
    l=0       
    for n = 1:12
        m = m+1;
        ref_yymm{m}=[num2str(i+1979) '-' num2str(n,'%02d')];
    for j = 1:eom_d(i,n)
        l=l+1; % index for 1~31(or 30 or 29)
        k=k+1; % index for 1~366
        ref_yymmdd{k}=[num2str(i+1979) '-' num2str(n,'%02d') '-'  num2str(j,'%02d')];
        ref_mmdd{k}=[num2str(n,'%02d') '-'  num2str(j,'%02d')];
        ref_mm{k}=[num2str(n,'%02d')];
        ref_yy{k}=[num2str(i+1979)];
    end
    end
end

k=0; m=0;
for i = 1:length(2001:2010)
    l=0       
    for n = 1:12
        m = m+1;
        ref_yymm_2nd{m}=[num2str(i+2000) '-' num2str(n,'%02d')];
    for j = 1:eom_d_2nd(i,n)
        l=l+1; % index for 1~31(or 30 or 29)
        k=k+1; % index for 1~366
        ref_yymmdd_2nd{k}=[num2str(i+2000) '-' num2str(n,'%02d') '-'  num2str(j,'%02d')];
        ref_mmdd_2nd{k}=[num2str(n,'%02d') '-'  num2str(j,'%02d')];
        ref_mm_2nd{k}=[num2str(n,'%02d')];
        ref_yy_2nd{k}=[num2str(i+2000)];
    end
    end
end



% matched date 'yymm' form on surf & bottom
% for j = 1:length(name_tag_1) % st. axis
for i = 1:12 %month
   id_mn{i} = find(strcmp(ref_mm, num2str(i,'%02d')) == 1);
end

for i = 1:12 %month
   id_mn_2nd{i} = find(strcmp(ref_mm_2nd, num2str(i,'%02d')) == 1);
end


% sptatial value

% len_temp=size(obs_gy_no3_ft_3(:,i:12:end),1)*size(obs_gy_no3_ft_3(:,i:12:end),2);
clearvars mod_mn_*
%% monthly climate mean
for i = 1:size(mer_sp_no3_1,1)
    for j=1:12 %month
    %% 1st    
    mod_mn_no3_1(i,j) = mean(mer_sp_no3_1(L(i),id_mn{j}),2);
    mod_mn_nh4_1(i,j) = mean(mer_sp_nh4_1(L(i),id_mn{j}),2);
    mod_mn_chl_1(i,j) = mean(mer_sp_chl_1(L(i),id_mn{j}),2);
    mod_mn_temp_1(i,j) = mean(mer_sp_temp_1(L(i),id_mn{j}),2);
    mod_mn_salt_1(i,j) = mean(mer_sp_salt_1(L(i),id_mn{j}),2);
    mod_mn_po4_1(i,j) = mean(mer_sp_po4_1(L(i),id_mn{j}),2);

    %bot merged
    mod_mn_no3_b_1(i,j) = mean(mer_sp_no3_b_1(L(i),id_mn{j}),2);
    mod_mn_nh4_b_1(i,j) = mean(mer_sp_nh4_b_1(L(i),id_mn{j}),2);
    mod_mn_chl_b_1(i,j) = mean(mer_sp_chl_b_1(L(i),id_mn{j}),2);
    mod_mn_temp_b_1(i,j) = mean(mer_sp_temp_b_1(L(i),id_mn{j}),2);
    mod_mn_salt_b_11(i,j) = mean(mer_sp_salt_b_1(L(i),id_mn{j}),2);
    mod_mn_po4_b_1(i,j) = mean(mer_sp_po4_b_1(L(i),id_mn{j}),2);

    %% 2nd 
    mod_mn_no3_2(i,j) = mean(mer_sp_no3_2(L(i),id_mn{j}),2);
    mod_mn_nh4_2(i,j) = mean(mer_sp_nh4_2(L(i),id_mn{j}),2);
    mod_mn_chl_2 (i,j) = mean(mer_sp_chl_2 (L(i),id_mn{j}),2);
    mod_mn_temp_2(i,j) = mean(mer_sp_temp_2(L(i),id_mn{j}),2);
    mod_mn_salt_2(i,j) = mean(mer_sp_salt_2(L(i),id_mn{j}),2);
    mod_mn_po4_2(i,j) = mean(mer_sp_po4_2(L(i),id_mn{j}),2);

    %bot merged
    mod_mn_no3_2(i,j) = mean(mer_sp_no3_b_2(L(i),id_mn{j}),2);
    mod_mn_nh4_b_2(i,j) = mean(mer_sp_nh4_b_2(L(i),id_mn{j}),2);
    mod_mn_chl_b_2(i,j) = mean(mer_sp_chl_b_2(L(i),id_mn{j}),2);
    mod_mn_temp_b_2(i,j) = mean(mer_sp_temp_b_2(L(i),id_mn{j}),2);
    mod_mn_salt_b_2(i,j) = mean(mer_sp_salt_b_2(L(i),id_mn{j}),2);
    mod_mn_po4_b_2(i,j) = mean(mer_sp_po4_b_2(L(i),id_mn{j}),2);

end
end

%% monthly climate std






for i= 1:12
 %% each st. seansonal climate compare (monthly)
        %each st time std
for jst = 1:size(obs_gy_no3_ft_1,1)  %st num
        eachst_std_gy_no3_ft_1(jst,i)=nanstd(obs_gy_no3_ft_1(jst,i:12:end));
        eachst_std_gy_no3_b_ft_1(jst,i)=nanstd(obs_gy_no3_b_ft_1(jst,i:12:end));
        eachst_std_gy_nh4_ft_1(jst,i)=nanstd(obs_gy_nh4_ft_1(jst,i:12:end));
        eachst_std_gy_nh4_b_ft_1(jst,i)=nanstd(obs_gy_nh4_b_ft_1(jst,i:12:end));
        eachst_std_gy_chl_ft_1(jst,i)=nanstd(obs_gy_chl_ft_1(jst,i:12:end));
        eachst_std_gy_chl_b_ft_1(jst,i)=nanstd(obs_gy_chl_b_ft_1(jst,i:12:end));
        eachst_std_gy_temp_ft_1(jst,i)=nanstd(obs_gy_temp_ft_1(jst,i:12:end));
        eachst_std_gy_temp_b_ft_1(jst,i)=nanstd(obs_gy_temp_b_ft_1(jst,i:12:end));
        eachst_std_gy_salt_ft_1(jst,i)=nanstd(obs_gy_salt_ft_1(jst,i:12:end));
        eachst_std_gy_salt_b_ft_1(jst,i)=nanstd(obs_gy_salt_b_ft_1(jst,i:12:end));
        eachst_std_gy_do_ft_1(jst,i)=nanstd(obs_gy_do_ft_1(jst,i:12:end));
        eachst_std_gy_do_b_ft_1(jst,i)=nanstd(obs_gy_do_b_ft_1(jst,i:12:end));
        eachst_std_gy_po4_ft_1(jst,i)=nanstd(obs_gy_po4_ft_1(jst,i:12:end));
        eachst_std_gy_po4_b_ft_1(jst,i)=nanstd(obs_gy_po4_b_ft_1(jst,i:12:end));

        eachst_std_gy_no3_ft_2(jst,i)=nanstd(obs_gy_no3_ft_2(jst,i:12:end));
        eachst_std_gy_no3_b_ft_2(jst,i)=nanstd(obs_gy_no3_b_ft_2(jst,i:12:end));
        eachst_std_gy_nh4_ft_2(jst,i)=nanstd(obs_gy_nh4_ft_2(jst,i:12:end));
        eachst_std_gy_nh4_b_ft_2(jst,i)=nanstd(obs_gy_nh4_b_ft_2(jst,i:12:end));
        eachst_std_gy_chl_ft_2(jst,i)=nanstd(obs_gy_chl_ft_2(jst,i:12:end));
        eachst_std_gy_chl_b_ft_2(jst,i)=nanstd(obs_gy_chl_b_ft_2(jst,i:12:end));
        eachst_std_gy_temp_ft_2(jst,i)=nanstd(obs_gy_temp_ft_2(jst,i:12:end));
        eachst_std_gy_temp_b_ft_2(jst,i)=nanstd(obs_gy_temp_b_ft_2(jst,i:12:end));
        eachst_std_gy_salt_ft_2(jst,i)=nanstd(obs_gy_salt_ft_2(jst,i:12:end));
        eachst_std_gy_salt_b_ft_2(jst,i)=nanstd(obs_gy_salt_b_ft_2(jst,i:12:end));
        eachst_std_gy_do_ft_2(jst,i)=nanstd(obs_gy_do_ft_2(jst,i:12:end));
        eachst_std_gy_do_b_ft_2(jst,i)=nanstd(obs_gy_do_b_ft_2(jst,i:12:end));        
        eachst_std_gy_po4_ft_2(jst,i)=nanstd(obs_gy_po4_ft_2(jst,i:12:end));
        eachst_std_gy_po4_b_ft_2(jst,i)=nanstd(obs_gy_po4_b_ft_2(jst,i:12:end));
        
        eachst_std_gy_no3_ft_3(jst,i)=nanstd(obs_gy_no3_ft_3(jst,i:12:end));
        eachst_std_gy_no3_b_ft_3(jst,i)=nanstd(obs_gy_no3_b_ft_3(jst,i:12:end));
        eachst_std_gy_nh4_ft_3(jst,i)=nanstd(obs_gy_nh4_ft_3(jst,i:12:end));
        eachst_std_gy_nh4_b_ft_3(jst,i)=nanstd(obs_gy_nh4_b_ft_3(jst,i:12:end));
        eachst_std_gy_chl_ft_3(jst,i)=nanstd(obs_gy_chl_ft_3(jst,i:12:end));
        eachst_std_gy_chl_b_ft_3(jst,i)=nanstd(obs_gy_chl_b_ft_3(jst,i:12:end));
        eachst_std_gy_temp_ft_3(jst,i)=nanstd(obs_gy_temp_ft_3(jst,i:12:end));
        eachst_std_gy_temp_b_ft_3(jst,i)=nanstd(obs_gy_temp_b_ft_3(jst,i:12:end));
        eachst_std_gy_salt_ft_3(jst,i)=nanstd(obs_gy_salt_ft_3(jst,i:12:end));
        eachst_std_gy_salt_b_ft_3(jst,i)=nanstd(obs_gy_salt_b_ft_3(jst,i:12:end));
        eachst_std_gy_do_ft_3(jst,i)=nanstd(obs_gy_do_ft_3(jst,i:12:end));
        eachst_std_gy_do_b_ft_3(jst,i)=nanstd(obs_gy_do_b_ft_3(jst,i:12:end));        
        eachst_std_gy_po4_ft_3(jst,i)=nanstd(obs_gy_po4_ft_3(jst,i:12:end));
        eachst_std_gy_po4_b_ft_3(jst,i)=nanstd(obs_gy_po4_b_ft_3(jst,i:12:end));

        %each st monthly mean climate
        eachst_clim_gy_no3_ft_1(jst,i)=nanmean(obs_gy_no3_ft_1(jst,i:12:end));
        eachst_clim_gy_no3_b_ft_1(jst,i)=nanmean(obs_gy_no3_b_ft_1(jst,i:12:end));
        eachst_clim_gy_nh4_ft_1(jst,i)=nanmean(obs_gy_nh4_ft_1(jst,i:12:end));
        eachst_clim_gy_nh4_b_ft_1(jst,i)=nanmean(obs_gy_nh4_b_ft_1(jst,i:12:end));
        eachst_clim_gy_chl_ft_1(jst,i)=nanmean(obs_gy_chl_ft_1(jst,i:12:end));
        eachst_clim_gy_chl_b_ft_1(jst,i)=nanmean(obs_gy_chl_b_ft_1(jst,i:12:end));
        eachst_clim_gy_temp_ft_1(jst,i)=nanmean(obs_gy_temp_ft_1(jst,i:12:end));
        eachst_clim_gy_temp_b_ft_1(jst,i)=nanmean(obs_gy_temp_b_ft_1(jst,i:12:end));
        eachst_clim_gy_salt_ft_1(jst,i)=nanmean(obs_gy_salt_ft_1(jst,i:12:end));
        eachst_clim_gy_salt_b_ft_1(jst,i)=nanmean(obs_gy_salt_b_ft_1(jst,i:12:end));
        eachst_clim_gy_do_ft_1(jst,i)=nanmean(obs_gy_do_ft_1(jst,i:12:end));
        eachst_clim_gy_do_b_ft_1(jst,i)=nanmean(obs_gy_do_b_ft_1(jst,i:12:end));
        eachst_clim_gy_po4_ft_1(jst,i)=nanmean(obs_gy_po4_ft_1(jst,i:12:end));
        eachst_clim_gy_po4_b_ft_1(jst,i)=nanmean(obs_gy_po4_b_ft_1(jst,i:12:end));

        eachst_clim_gy_no3_ft_2(jst,i)=nanmean(obs_gy_no3_ft_2(jst,i:12:end));
        eachst_clim_gy_no3_b_ft_2(jst,i)=nanmean(obs_gy_no3_b_ft_2(jst,i:12:end));
        eachst_clim_gy_nh4_ft_2(jst,i)=nanmean(obs_gy_nh4_ft_2(jst,i:12:end));
        eachst_clim_gy_nh4_b_ft_2(jst,i)=nanmean(obs_gy_nh4_b_ft_2(jst,i:12:end));
        eachst_clim_gy_chl_ft_2(jst,i)=nanmean(obs_gy_chl_ft_2(jst,i:12:end));
        eachst_clim_gy_chl_b_ft_2(jst,i)=nanmean(obs_gy_chl_b_ft_2(jst,i:12:end));
        eachst_clim_gy_temp_ft_2(jst,i)=nanmean(obs_gy_temp_ft_2(jst,i:12:end));
        eachst_clim_gy_temp_b_ft_2(jst,i)=nanmean(obs_gy_temp_b_ft_2(jst,i:12:end));
        eachst_clim_gy_salt_ft_2(jst,i)=nanmean(obs_gy_salt_ft_2(jst,i:12:end));
        eachst_clim_gy_salt_b_ft_2(jst,i)=nanmean(obs_gy_salt_b_ft_2(jst,i:12:end));
        eachst_clim_gy_do_ft_2(jst,i)=nanmean(obs_gy_do_ft_2(jst,i:12:end));
        eachst_clim_gy_do_b_ft_2(jst,i)=nanmean(obs_gy_do_b_ft_2(jst,i:12:end));        
        eachst_clim_gy_po4_ft_2(jst,i)=nanmean(obs_gy_po4_ft_2(jst,i:12:end));
        eachst_clim_gy_po4_b_ft_2(jst,i)=nanmean(obs_gy_po4_b_ft_2(jst,i:12:end));  
        
        eachst_clim_gy_no3_ft_3(jst,i)=nanmean(obs_gy_no3_ft_3(jst,i:12:end));
        eachst_clim_gy_no3_b_ft_3(jst,i)=nanmean(obs_gy_no3_b_ft_3(jst,i:12:end));
        eachst_clim_gy_nh4_ft_3(jst,i)=nanmean(obs_gy_nh4_ft_3(jst,i:12:end));
        eachst_clim_gy_nh4_b_ft_3(jst,i)=nanmean(obs_gy_nh4_b_ft_3(jst,i:12:end));
        eachst_clim_gy_chl_ft_3(jst,i)=nanmean(obs_gy_chl_ft_3(jst,i:12:end));
        eachst_clim_gy_chl_b_ft_3(jst,i)=nanmean(obs_gy_chl_b_ft_3(jst,i:12:end));
        eachst_clim_gy_temp_ft_3(jst,i)=nanmean(obs_gy_temp_ft_3(jst,i:12:end));
        eachst_clim_gy_temp_b_ft_3(jst,i)=nanmean(obs_gy_temp_b_ft_3(jst,i:12:end));
        eachst_clim_gy_salt_ft_3(jst,i)=nanmean(obs_gy_salt_ft_3(jst,i:12:end));
        eachst_clim_gy_salt_b_ft_3(jst,i)=nanmean(obs_gy_salt_b_ft_3(jst,i:12:end));
        eachst_clim_gy_do_ft_3(jst,i)=nanmean(obs_gy_do_ft_3(jst,i:12:end));
        eachst_clim_gy_do_b_ft_3(jst,i)=nanmean(obs_gy_do_b_ft_3(jst,i:12:end));        
        eachst_clim_gy_po4_ft_3(jst,i)=nanmean(obs_gy_po4_ft_3(jst,i:12:end));
        eachst_clim_gy_po4_b_ft_3(jst,i)=nanmean(obs_gy_po4_b_ft_3(jst,i:12:end));
end
end





cd D:\Àå±â»ýÅÂ\Dynamic\KOEM\compare_sewer\eachst_climate\3sig
for i = 1:size(koem_3sig.eachst_clim_gy_no3_ft_1,1)
    
%% surface
    %% NO3
    if i==1
        fig=figure;hold on;
        plot(2:6:12,koem_3sig.eachst_clim_gy_no3_ft_1(i,2:6:12)./N_MW,'ro'); 
        plot(2:6:12,koem_3sig.eachst_clim_gy_no3_ft_2(i,2:6:12)./N_MW,'go');
        plot(2:6:12,koem_3sig.eachst_clim_gy_no3_ft_3(i,2:6:12)./N_MW,'bo'); 
        plot(2:6:12,koem_3sig.eachst_clim_gy_no3_ft_1(i,2:6:12)./N_MW,'r'); 
        plot(2:6:12,koem_3sig.eachst_clim_gy_no3_ft_2(i,2:6:12)./N_MW,'g');
        plot(2:6:12,koem_3sig.eachst_clim_gy_no3_ft_3(i,2:6:12)./N_MW,'b'); 
        grid on;
            errorbar(1:12,koem_3sig.eachst_clim_gy_no3_ft_1(i,:)./N_MW,koem_3sig.eachst_std_gy_no3_ft_1(i,:)./N_MW,'r'); 
            errorbar(1:12,koem_3sig.eachst_clim_gy_no3_ft_2(i,:)./N_MW,koem_3sig.eachst_std_gy_no3_ft_2(i,:)./N_MW,'g');
            errorbar(1:12,koem_3sig.eachst_clim_gy_no3_ft_3(i,:)./N_MW,koem_3sig.eachst_std_gy_no3_ft_3(i,:)./N_MW,'b'); 
        alpha(0.3); 
        xlim([1 12]); ylim([0 24])
        title(['NO3- 3sig - ',nn(i,:),' st. seasonal climate']);
        xlabel('month')
        ylabel('mmol/m3')
%         legend('1997-2000','2001-2010','2011-2019');
        print(fig,['NO3_3sig_',nn(i,:),'_st._seasonal_climate.png'],'-dpng')
        
    else
        fig=figure;hold on;
        plot(2:3:12,koem_3sig.eachst_clim_gy_no3_ft_1(i,2:3:12)./N_MW,'ro'); 
        plot(2:3:12,koem_3sig.eachst_clim_gy_no3_ft_2(i,2:3:12)./N_MW,'go');
        plot(2:3:12,koem_3sig.eachst_clim_gy_no3_ft_3(i,2:3:12)./N_MW,'bo'); 
        plot(2:3:12,koem_3sig.eachst_clim_gy_no3_ft_1(i,2:3:12)./N_MW,'r'); 
        plot(2:3:12,koem_3sig.eachst_clim_gy_no3_ft_2(i,2:3:12)./N_MW,'g');
        plot(2:3:12,koem_3sig.eachst_clim_gy_no3_ft_3(i,2:3:12)./N_MW,'b'); 
        grid on;
            errorbar(1:12,koem_3sig.eachst_clim_gy_no3_ft_1(i,:)./N_MW, koem_3sig.eachst_std_gy_no3_ft_1(i,:)./N_MW,'r'); 
            errorbar(1:12,koem_3sig.eachst_clim_gy_no3_ft_2(i,:)./N_MW, koem_3sig.eachst_std_gy_no3_ft_2(i,:)./N_MW,'g');
            errorbar(1:12,koem_3sig.eachst_clim_gy_no3_ft_3(i,:)./N_MW, koem_3sig.eachst_std_gy_no3_ft_3(i,:)./N_MW,'b'); 
        alpha(0.3)
        xlim([1 12]); ylim([0 24])
        title(['NO3- 3sig - ',nn(i,:),' st. seasonal climate']);
        xlabel('month')
        ylabel('mmol/m3')
        print(fig,['NO3_3sig_',nn(i,:),'_st._seasonal_climate.png'],'-dpng')
    end
    
        %% NH4
    if i==1
        fig=figure;hold on;
        plot(2:6:12,koem_3sig.eachst_clim_gy_nh4_ft_1(i,2:6:12)./N_MW,'ro'); 
        plot(2:6:12,koem_3sig.eachst_clim_gy_nh4_ft_2(i,2:6:12)./N_MW,'go');
        plot(2:6:12,koem_3sig.eachst_clim_gy_nh4_ft_3(i,2:6:12)./N_MW,'bo'); 
        plot(2:6:12,koem_3sig.eachst_clim_gy_nh4_ft_1(i,2:6:12)./N_MW,'r'); 
        plot(2:6:12,koem_3sig.eachst_clim_gy_nh4_ft_2(i,2:6:12)./N_MW,'g');
        plot(2:6:12,koem_3sig.eachst_clim_gy_nh4_ft_3(i,2:6:12)./N_MW,'b'); 
        grid on;
            errorbar(1:12,koem_3sig.eachst_clim_gy_nh4_ft_1(i,:)./N_MW,koem_3sig.eachst_std_gy_nh4_ft_1(i,:)./N_MW,'r'); 
            errorbar(1:12,koem_3sig.eachst_clim_gy_nh4_ft_2(i,:)./N_MW,koem_3sig.eachst_std_gy_nh4_ft_2(i,:)./N_MW,'g');
            errorbar(1:12,koem_3sig.eachst_clim_gy_nh4_ft_3(i,:)./N_MW,koem_3sig.eachst_std_gy_nh4_ft_3(i,:)./N_MW,'b'); 
        alpha(0.3); 
        xlim([1 12]); ylim([0 12])
        title(['NH4- 3sig - ',nn(i,:),' st. seasonal climate']);
        xlabel('month')
        ylabel('mmol/m3')
        print(fig,['NH4_3sig_',nn(i,:),'_st._seasonal_climate.png'],'-dpng')
    else
        fig=figure;hold on;
        plot(2:3:12,koem_3sig.eachst_clim_gy_nh4_ft_1(i,2:3:12)./N_MW,'ro'); 
        plot(2:3:12,koem_3sig.eachst_clim_gy_nh4_ft_2(i,2:3:12)./N_MW,'go');
        plot(2:3:12,koem_3sig.eachst_clim_gy_nh4_ft_3(i,2:3:12)./N_MW,'bo'); 
        plot(2:3:12,koem_3sig.eachst_clim_gy_nh4_ft_1(i,2:3:12)./N_MW,'r'); 
        plot(2:3:12,koem_3sig.eachst_clim_gy_nh4_ft_2(i,2:3:12)./N_MW,'g');
        plot(2:3:12,koem_3sig.eachst_clim_gy_nh4_ft_3(i,2:3:12)./N_MW,'b'); 
        grid on;
            errorbar(1:12,koem_3sig.eachst_clim_gy_nh4_ft_1(i,:)./N_MW,koem_3sig.eachst_std_gy_nh4_ft_1(i,:)./N_MW,'r'); 
            errorbar(1:12,koem_3sig.eachst_clim_gy_nh4_ft_2(i,:)./N_MW,koem_3sig.eachst_std_gy_nh4_ft_2(i,:)./N_MW,'g');
            errorbar(1:12,koem_3sig.eachst_clim_gy_nh4_ft_3(i,:)./N_MW,koem_3sig.eachst_std_gy_nh4_ft_3(i,:)./N_MW,'b'); 
        alpha(0.3)
        xlim([1 12]); ylim([0 12])
        title(['NH4- 3sig - ',nn(i,:),' st. seasonal climate']);
        xlabel('month')
        ylabel('mmol/m3')
        print(fig,['NH4_3sig_',nn(i,:),'_st._seasonal_climate.png'],'-dpng')
    end
    
        %% PO4   
    if i==1
        fig=figure;hold on;
        plot(2:6:12,koem_3sig.eachst_clim_gy_po4_ft_1(i,2:6:12)./P_MW,'ro'); 
        plot(2:6:12,koem_3sig.eachst_clim_gy_po4_ft_2(i,2:6:12)./P_MW,'go');
        plot(2:6:12,koem_3sig.eachst_clim_gy_po4_ft_3(i,2:6:12)./P_MW,'bo'); 
        plot(2:6:12,koem_3sig.eachst_clim_gy_po4_ft_1(i,2:6:12)./P_MW,'r'); 
        plot(2:6:12,koem_3sig.eachst_clim_gy_po4_ft_2(i,2:6:12)./P_MW,'g');
        plot(2:6:12,koem_3sig.eachst_clim_gy_po4_ft_3(i,2:6:12)./P_MW,'b'); 
        grid on;
            errorbar(1:12,koem_3sig.eachst_clim_gy_po4_ft_1(i,:)./P_MW,koem_3sig.eachst_std_gy_po4_ft_1(i,:)./P_MW,'r'); 
            errorbar(1:12,koem_3sig.eachst_clim_gy_po4_ft_2(i,:)./P_MW,koem_3sig.eachst_std_gy_po4_ft_2(i,:)./P_MW,'g');
            errorbar(1:12,koem_3sig.eachst_clim_gy_po4_ft_3(i,:)./P_MW,koem_3sig.eachst_std_gy_po4_ft_3(i,:)./P_MW,'b'); 
        alpha(0.3); 
        xlim([1 12]); ylim([0 180./P_MW])
        title(['PO4- 3sig - ',nn(i,:),' st. seasonal climate']);
        xlabel('month')
        ylabel('mmol/m3')
        print(fig,['PO4_3sig_',nn(i,:),'_st._seasonal_climate.png'],'-dpng')
    else
        fig=figure;hold on;
        plot(2:3:12,koem_3sig.eachst_clim_gy_po4_ft_1(i,2:3:12)./P_MW,'ro'); 
        plot(2:3:12,koem_3sig.eachst_clim_gy_po4_ft_2(i,2:3:12)./P_MW,'go');
        plot(2:3:12,koem_3sig.eachst_clim_gy_po4_ft_3(i,2:3:12)./P_MW,'bo'); 
        plot(2:3:12,koem_3sig.eachst_clim_gy_po4_ft_1(i,2:3:12)./P_MW,'r'); 
        plot(2:3:12,koem_3sig.eachst_clim_gy_po4_ft_2(i,2:3:12)./P_MW,'g');
        plot(2:3:12,koem_3sig.eachst_clim_gy_po4_ft_3(i,2:3:12)./P_MW,'b'); 
        grid on;
            errorbar(1:12,koem_3sig.eachst_clim_gy_po4_ft_1(i,:)./P_MW,koem_3sig.eachst_std_gy_po4_ft_1(i,:)./P_MW,'r'); 
            errorbar(1:12,koem_3sig.eachst_clim_gy_po4_ft_2(i,:)./P_MW,koem_3sig.eachst_std_gy_po4_ft_2(i,:)./P_MW,'g');
            errorbar(1:12,koem_3sig.eachst_clim_gy_po4_ft_3(i,:)./P_MW,koem_3sig.eachst_std_gy_po4_ft_3(i,:)./P_MW,'b'); 
        alpha(0.3)
        xlim([1 12]); ylim([0 180./P_MW])
        title(['PO4- 3sig - ',nn(i,:),' st. seasonal climate']);
        xlabel('month')
        ylabel('mmol/m3')
        print(fig,['PO4_3sig_',nn(i,:),'_st._seasonal_climate.png'],'-dpng')
    end
    
        %% Chl    
    if i==1
        fig=figure;hold on;
        plot(2:6:12,koem_3sig.eachst_clim_gy_chl_ft_1(i,2:6:12),'ro'); 
        plot(2:6:12,koem_3sig.eachst_clim_gy_chl_ft_2(i,2:6:12),'go');
        plot(2:6:12,koem_3sig.eachst_clim_gy_chl_ft_3(i,2:6:12),'bo'); 
        plot(2:6:12,koem_3sig.eachst_clim_gy_chl_ft_1(i,2:6:12),'r'); 
        plot(2:6:12,koem_3sig.eachst_clim_gy_chl_ft_2(i,2:6:12),'g');
        plot(2:6:12,koem_3sig.eachst_clim_gy_chl_ft_3(i,2:6:12),'b'); 
        grid on;
            errorbar(1:12,koem_3sig.eachst_clim_gy_chl_ft_1(i,:),koem_3sig.eachst_std_gy_chl_ft_1(i,:),'r'); 
            errorbar(1:12,koem_3sig.eachst_clim_gy_chl_ft_2(i,:),koem_3sig.eachst_std_gy_chl_ft_2(i,:),'g');
            errorbar(1:12,koem_3sig.eachst_clim_gy_chl_ft_3(i,:),koem_3sig.eachst_std_gy_chl_ft_3(i,:),'b'); 
        alpha(0.3); 
        xlim([1 12]); ylim([0 14])
        title(['chl- 3sig - ',nn(i,:),' st. seasonal climate']);
        xlabel('month')
        ylabel('ug/L')
        print(fig,['chl_3sig_',nn(i,:),'_st._seasonal_climate.png'],'-dpng')
    else
        fig=figure;hold on;
        plot(2:3:12,koem_3sig.eachst_clim_gy_chl_ft_1(i,2:3:12),'ro'); 
        plot(2:3:12,koem_3sig.eachst_clim_gy_chl_ft_2(i,2:3:12),'go');
        plot(2:3:12,koem_3sig.eachst_clim_gy_chl_ft_3(i,2:3:12),'bo'); 
        plot(2:3:12,koem_3sig.eachst_clim_gy_chl_ft_1(i,2:3:12),'r'); 
        plot(2:3:12,koem_3sig.eachst_clim_gy_chl_ft_2(i,2:3:12),'g');
        plot(2:3:12,koem_3sig.eachst_clim_gy_chl_ft_3(i,2:3:12),'b'); 
        grid on;
            errorbar(1:12,koem_3sig.eachst_clim_gy_chl_ft_1(i,:),koem_3sig.eachst_std_gy_chl_ft_1(i,:),'r'); 
            errorbar(1:12,koem_3sig.eachst_clim_gy_chl_ft_2(i,:),koem_3sig.eachst_std_gy_chl_ft_2(i,:),'g');
            errorbar(1:12,koem_3sig.eachst_clim_gy_chl_ft_3(i,:),koem_3sig.eachst_std_gy_chl_ft_3(i,:),'b'); 
        alpha(0.3)
        xlim([1 12]); ylim([0 14])
        title(['chl- 3sig - ',nn(i,:),' st. seasonal climate']);
        xlabel('month')
        ylabel('ug/L')
        print(fig,['chl_3sig_',nn(i,:),'_st._seasonal_climate.png'],'-dpng')
    end
    
%% Bottom
    %% NO3
    if i==1
        fig=figure;hold on;
        plot(2:6:12,koem_3sig.eachst_clim_gy_no3_b_ft_1(i,2:6:12)./N_MW,'ro'); 
        plot(2:6:12,koem_3sig.eachst_clim_gy_no3_b_ft_2(i,2:6:12)./N_MW,'go');
        plot(2:6:12,koem_3sig.eachst_clim_gy_no3_b_ft_3(i,2:6:12)./N_MW,'bo'); 
        plot(2:6:12,koem_3sig.eachst_clim_gy_no3_b_ft_1(i,2:6:12)./N_MW,'r'); 
        plot(2:6:12,koem_3sig.eachst_clim_gy_no3_b_ft_2(i,2:6:12)./N_MW,'g');
        plot(2:6:12,koem_3sig.eachst_clim_gy_no3_b_ft_3(i,2:6:12)./N_MW,'b'); 
        grid on;
            errorbar(1:12,koem_3sig.eachst_clim_gy_no3_b_ft_1(i,:)./N_MW,koem_3sig.eachst_std_gy_no3_b_ft_1(i,:)./N_MW,'r'); 
            errorbar(1:12,koem_3sig.eachst_clim_gy_no3_b_ft_2(i,:)./N_MW,koem_3sig.eachst_std_gy_no3_b_ft_2(i,:)./N_MW,'g');
            errorbar(1:12,koem_3sig.eachst_clim_gy_no3_b_ft_3(i,:)./N_MW,koem_3sig.eachst_std_gy_no3_b_ft_3(i,:)./N_MW,'b'); 
        alpha(0.3); 
        xlim([1 12]); ylim([0 24])
        title(['NO3- 3sig - ',nn(i,:),' st. seasonal climate bot']);
        xlabel('month')
        ylabel('mmol/m3')
        print(fig,['NO3_3sig_',nn(i,:),'_st._seasonal_climate_bot.png'],'-dpng')
        
    else
        fig=figure;hold on;
        plot(2:3:12,koem_3sig.eachst_clim_gy_no3_b_ft_1(i,2:3:12)./N_MW,'ro'); 
        plot(2:3:12,koem_3sig.eachst_clim_gy_no3_b_ft_2(i,2:3:12)./N_MW,'go');
        plot(2:3:12,koem_3sig.eachst_clim_gy_no3_b_ft_3(i,2:3:12)./N_MW,'bo'); 
        plot(2:3:12,koem_3sig.eachst_clim_gy_no3_b_ft_1(i,2:3:12)./N_MW,'r'); 
        plot(2:3:12,koem_3sig.eachst_clim_gy_no3_b_ft_2(i,2:3:12)./N_MW,'g');
        plot(2:3:12,koem_3sig.eachst_clim_gy_no3_b_ft_3(i,2:3:12)./N_MW,'b'); 
        grid on;
            errorbar(1:12,koem_3sig.eachst_clim_gy_no3_b_ft_1(i,:)./N_MW, koem_3sig.eachst_std_gy_no3_b_ft_1(i,:)./N_MW,'r'); 
            errorbar(1:12,koem_3sig.eachst_clim_gy_no3_b_ft_2(i,:)./N_MW, koem_3sig.eachst_std_gy_no3_b_ft_2(i,:)./N_MW,'g');
            errorbar(1:12,koem_3sig.eachst_clim_gy_no3_b_ft_3(i,:)./N_MW, koem_3sig.eachst_std_gy_no3_b_ft_3(i,:)./N_MW,'b'); 
        alpha(0.3)
        xlim([1 12]); ylim([0 24])
        title(['NO3- 3sig - ',nn(i,:),' st. seasonal climate bot']);
        xlabel('month')
        ylabel('mmol/m3')
        print(fig,['NO3_3sig_',nn(i,:),'_st._seasonal_climate_bot.png'],'-dpng')
    end
    
        %% NH4
    if i==1
        fig=figure;hold on;
        plot(2:6:12,koem_3sig.eachst_clim_gy_nh4_b_ft_1(i,2:6:12)./N_MW,'ro'); 
        plot(2:6:12,koem_3sig.eachst_clim_gy_nh4_b_ft_2(i,2:6:12)./N_MW,'go');
        plot(2:6:12,koem_3sig.eachst_clim_gy_nh4_b_ft_3(i,2:6:12)./N_MW,'bo'); 
        plot(2:6:12,koem_3sig.eachst_clim_gy_nh4_b_ft_1(i,2:6:12)./N_MW,'r'); 
        plot(2:6:12,koem_3sig.eachst_clim_gy_nh4_b_ft_2(i,2:6:12)./N_MW,'g');
        plot(2:6:12,koem_3sig.eachst_clim_gy_nh4_b_ft_3(i,2:6:12)./N_MW,'b'); 
        grid on;
            errorbar(1:12,koem_3sig.eachst_clim_gy_nh4_b_ft_1(i,:)./N_MW,koem_3sig.eachst_std_gy_nh4_b_ft_1(i,:)./N_MW,'r'); 
            errorbar(1:12,koem_3sig.eachst_clim_gy_nh4_b_ft_2(i,:)./N_MW,koem_3sig.eachst_std_gy_nh4_b_ft_2(i,:)./N_MW,'g');
            errorbar(1:12,koem_3sig.eachst_clim_gy_nh4_b_ft_3(i,:)./N_MW,koem_3sig.eachst_std_gy_nh4_b_ft_3(i,:)./N_MW,'b'); 
        alpha(0.3); 
        xlim([1 12]); ylim([0 12])
        title(['NH4- 3sig - ',nn(i,:),' st. seasonal climate bot']);
        xlabel('month')
        ylabel('mmol/m3')
        print(fig,['NH4_3sig_',nn(i,:),'_st._seasonal_climate_bot.png'],'-dpng')
    else
        fig=figure;hold on;
        plot(2:3:12,koem_3sig.eachst_clim_gy_nh4_b_ft_1(i,2:3:12)./N_MW,'ro'); 
        plot(2:3:12,koem_3sig.eachst_clim_gy_nh4_b_ft_2(i,2:3:12)./N_MW,'go');
        plot(2:3:12,koem_3sig.eachst_clim_gy_nh4_b_ft_3(i,2:3:12)./N_MW,'bo'); 
        plot(2:3:12,koem_3sig.eachst_clim_gy_nh4_b_ft_1(i,2:3:12)./N_MW,'r'); 
        plot(2:3:12,koem_3sig.eachst_clim_gy_nh4_b_ft_2(i,2:3:12)./N_MW,'g');
        plot(2:3:12,koem_3sig.eachst_clim_gy_nh4_b_ft_3(i,2:3:12)./N_MW,'b'); 
        grid on;
            errorbar(1:12,koem_3sig.eachst_clim_gy_nh4_b_ft_1(i,:)./N_MW,koem_3sig.eachst_std_gy_nh4_b_ft_1(i,:)./N_MW,'r'); 
            errorbar(1:12,koem_3sig.eachst_clim_gy_nh4_b_ft_2(i,:)./N_MW,koem_3sig.eachst_std_gy_nh4_b_ft_2(i,:)./N_MW,'g');
            errorbar(1:12,koem_3sig.eachst_clim_gy_nh4_b_ft_3(i,:)./N_MW,koem_3sig.eachst_std_gy_nh4_b_ft_3(i,:)./N_MW,'b'); 
        alpha(0.3)
        xlim([1 12]); ylim([0 12])
        title(['NH4- 3sig - ',nn(i,:),' st. seasonal climate bot']);
        xlabel('month')
        ylabel('mmol/m3')
        print(fig,['NH4_3sig_',nn(i,:),'_st._seasonal_climate_bot.png'],'-dpng')
    end
    
        %% PO4   
    if i==1
        fig=figure;hold on;
        plot(2:6:12,koem_3sig.eachst_clim_gy_po4_b_ft_1(i,2:6:12)./P_MW,'ro'); 
        plot(2:6:12,koem_3sig.eachst_clim_gy_po4_b_ft_2(i,2:6:12)./P_MW,'go');
        plot(2:6:12,koem_3sig.eachst_clim_gy_po4_b_ft_3(i,2:6:12)./P_MW,'bo'); 
        plot(2:6:12,koem_3sig.eachst_clim_gy_po4_b_ft_1(i,2:6:12)./P_MW,'r'); 
        plot(2:6:12,koem_3sig.eachst_clim_gy_po4_b_ft_2(i,2:6:12)./P_MW,'g');
        plot(2:6:12,koem_3sig.eachst_clim_gy_po4_b_ft_3(i,2:6:12)./P_MW,'b'); 
        grid on;
            errorbar(1:12,koem_3sig.eachst_clim_gy_po4_b_ft_1(i,:)./P_MW,koem_3sig.eachst_std_gy_po4_b_ft_1(i,:)./P_MW,'r'); 
            errorbar(1:12,koem_3sig.eachst_clim_gy_po4_b_ft_2(i,:)./P_MW,koem_3sig.eachst_std_gy_po4_b_ft_2(i,:)./P_MW,'g');
            errorbar(1:12,koem_3sig.eachst_clim_gy_po4_b_ft_3(i,:)./P_MW,koem_3sig.eachst_std_gy_po4_b_ft_3(i,:)./P_MW,'b'); 
        alpha(0.3); 
        xlim([1 12]); ylim([0 180./P_MW])
        title(['PO4- 3sig - ',nn(i,:),' st. seasonal climate bot']);
        xlabel('month')
        ylabel('mmol/m3')
        print(fig,['PO4_3sig_',nn(i,:),'_st._seasonal_climate_bot.png'],'-dpng')
    else
        fig=figure;hold on;
        plot(2:3:12,koem_3sig.eachst_clim_gy_po4_b_ft_1(i,2:3:12)./P_MW,'ro'); 
        plot(2:3:12,koem_3sig.eachst_clim_gy_po4_b_ft_2(i,2:3:12)./P_MW,'go');
        plot(2:3:12,koem_3sig.eachst_clim_gy_po4_b_ft_3(i,2:3:12)./P_MW,'bo'); 
        plot(2:3:12,koem_3sig.eachst_clim_gy_po4_b_ft_1(i,2:3:12)./P_MW,'r'); 
        plot(2:3:12,koem_3sig.eachst_clim_gy_po4_b_ft_2(i,2:3:12)./P_MW,'g');
        plot(2:3:12,koem_3sig.eachst_clim_gy_po4_b_ft_3(i,2:3:12)./P_MW,'b'); 
        grid on;
            errorbar(1:12,koem_3sig.eachst_clim_gy_po4_b_ft_1(i,:)./P_MW,koem_3sig.eachst_std_gy_po4_b_ft_1(i,:)./P_MW,'r'); 
            errorbar(1:12,koem_3sig.eachst_clim_gy_po4_b_ft_2(i,:)./P_MW,koem_3sig.eachst_std_gy_po4_b_ft_2(i,:)./P_MW,'g');
            errorbar(1:12,koem_3sig.eachst_clim_gy_po4_b_ft_3(i,:)./P_MW,koem_3sig.eachst_std_gy_po4_b_ft_3(i,:)./P_MW,'b'); 
        alpha(0.3)
        xlim([1 12]); ylim([0 180./P_MW])
        title(['PO4- 3sig - ',nn(i,:),' st. seasonal climate bot']);
        xlabel('month')
        ylabel('mmol/m3')
        print(fig,['PO4_3sig_',nn(i,:),'_st._seasonal_climate_bot.png'],'-dpng')
    end
    
        %% Chl    
    if i==1
        fig=figure;hold on;
        plot(2:6:12,koem_3sig.eachst_clim_gy_chl_b_ft_1(i,2:6:12),'ro'); 
        plot(2:6:12,koem_3sig.eachst_clim_gy_chl_b_ft_2(i,2:6:12),'go');
        plot(2:6:12,koem_3sig.eachst_clim_gy_chl_b_ft_3(i,2:6:12),'bo'); 
        plot(2:6:12,koem_3sig.eachst_clim_gy_chl_b_ft_1(i,2:6:12),'r'); 
        plot(2:6:12,koem_3sig.eachst_clim_gy_chl_b_ft_2(i,2:6:12),'g');
        plot(2:6:12,koem_3sig.eachst_clim_gy_chl_b_ft_3(i,2:6:12),'b'); 
        grid on;
            errorbar(1:12,koem_3sig.eachst_clim_gy_chl_b_ft_1(i,:),koem_3sig.eachst_std_gy_chl_b_ft_1(i,:),'r'); 
            errorbar(1:12,koem_3sig.eachst_clim_gy_chl_b_ft_2(i,:),koem_3sig.eachst_std_gy_chl_b_ft_2(i,:),'g');
            errorbar(1:12,koem_3sig.eachst_clim_gy_chl_b_ft_3(i,:),koem_3sig.eachst_std_gy_chl_b_ft_3(i,:),'b'); 
        alpha(0.3); 
        xlim([1 12]); ylim([0 14])
        title(['chl- 3sig - ',nn(i,:),' st. seasonal climate bot']);
        xlabel('month')
        ylabel('ug/L')
        print(fig,['chl_3sig_',nn(i,:),'_st._seasonal_climate_bot.png'],'-dpng')
    else
        fig=figure;hold on;
        plot(2:3:12,koem_3sig.eachst_clim_gy_chl_b_ft_1(i,2:3:12),'ro'); 
        plot(2:3:12,koem_3sig.eachst_clim_gy_chl_b_ft_2(i,2:3:12),'go');
        plot(2:3:12,koem_3sig.eachst_clim_gy_chl_b_ft_3(i,2:3:12),'bo'); 
        plot(2:3:12,koem_3sig.eachst_clim_gy_chl_b_ft_1(i,2:3:12),'r'); 
        plot(2:3:12,koem_3sig.eachst_clim_gy_chl_b_ft_2(i,2:3:12),'g');
        plot(2:3:12,koem_3sig.eachst_clim_gy_chl_b_ft_3(i,2:3:12),'b'); 
        grid on;
            errorbar(1:12,koem_3sig.eachst_clim_gy_chl_b_ft_1(i,:),koem_3sig.eachst_std_gy_chl_b_ft_1(i,:),'r'); 
            errorbar(1:12,koem_3sig.eachst_clim_gy_chl_b_ft_2(i,:),koem_3sig.eachst_std_gy_chl_b_ft_2(i,:),'g');
            errorbar(1:12,koem_3sig.eachst_clim_gy_chl_b_ft_3(i,:),koem_3sig.eachst_std_gy_chl_b_ft_3(i,:),'b'); 
        alpha(0.3)
        xlim([1 12]); ylim([0 14])
        title(['chl- 3sig - ',nn(i,:),' st. seasonal climate bot']);
        xlabel('month')
        ylabel('ug/L')
        print(fig,['chl_3sig_',nn(i,:),'_st._seasonal_climate_bot.png'],'-dpng')
    end
end
close all

cd D:\Àå±â»ýÅÂ\Dynamic\KOEM\compare_sewer\eachst_climate\2sig
for i = 1:size(koem_2sig.eachst_clim_gy_no3_ft_1,1)
    
%% surface
    %% NO3
    if i==1
        fig=figure;hold on;
        plot(2:6:12,koem_2sig.eachst_clim_gy_no3_ft_1(i,2:6:12)./N_MW,'ro'); 
        plot(2:6:12,koem_2sig.eachst_clim_gy_no3_ft_2(i,2:6:12)./N_MW,'go');
        plot(2:6:12,koem_2sig.eachst_clim_gy_no3_ft_3(i,2:6:12)./N_MW,'bo'); 
        plot(2:6:12,koem_2sig.eachst_clim_gy_no3_ft_1(i,2:6:12)./N_MW,'r'); 
        plot(2:6:12,koem_2sig.eachst_clim_gy_no3_ft_2(i,2:6:12)./N_MW,'g');
        plot(2:6:12,koem_2sig.eachst_clim_gy_no3_ft_3(i,2:6:12)./N_MW,'b'); 
        grid on;
            errorbar(1:12,koem_2sig.eachst_clim_gy_no3_ft_1(i,:)./N_MW,koem_2sig.eachst_std_gy_no3_ft_1(i,:)./N_MW,'r'); 
            errorbar(1:12,koem_2sig.eachst_clim_gy_no3_ft_2(i,:)./N_MW,koem_2sig.eachst_std_gy_no3_ft_2(i,:)./N_MW,'g');
            errorbar(1:12,koem_2sig.eachst_clim_gy_no3_ft_3(i,:)./N_MW,koem_2sig.eachst_std_gy_no3_ft_3(i,:)./N_MW,'b'); 
        alpha(0.3); 
        xlim([1 12]); ylim([0 24])
        title(['NO3- 2sig - ',nn(i,:),' st. seasonal climate']);
        xlabel('month')
        ylabel('mmol/m3')
        print(fig,['NO3_2sig_',nn(i,:),'_st._seasonal_climate.png'],'-dpng')
        
    else
        fig=figure;hold on;
        plot(2:3:12,koem_2sig.eachst_clim_gy_no3_ft_1(i,2:3:12)./N_MW,'ro'); 
        plot(2:3:12,koem_2sig.eachst_clim_gy_no3_ft_2(i,2:3:12)./N_MW,'go');
        plot(2:3:12,koem_2sig.eachst_clim_gy_no3_ft_3(i,2:3:12)./N_MW,'bo'); 
        plot(2:3:12,koem_2sig.eachst_clim_gy_no3_ft_1(i,2:3:12)./N_MW,'r'); 
        plot(2:3:12,koem_2sig.eachst_clim_gy_no3_ft_2(i,2:3:12)./N_MW,'g');
        plot(2:3:12,koem_2sig.eachst_clim_gy_no3_ft_3(i,2:3:12)./N_MW,'b'); 
        grid on;
            errorbar(1:12,koem_2sig.eachst_clim_gy_no3_ft_1(i,:)./N_MW, koem_2sig.eachst_std_gy_no3_ft_1(i,:)./N_MW,'r'); 
            errorbar(1:12,koem_2sig.eachst_clim_gy_no3_ft_2(i,:)./N_MW, koem_2sig.eachst_std_gy_no3_ft_2(i,:)./N_MW,'g');
            errorbar(1:12,koem_2sig.eachst_clim_gy_no3_ft_3(i,:)./N_MW, koem_2sig.eachst_std_gy_no3_ft_3(i,:)./N_MW,'b'); 
        alpha(0.3)
        xlim([1 12]); ylim([0 24])
        title(['NO3- 2sig - ',nn(i,:),' st. seasonal climate']);
        xlabel('month')
        ylabel('mmol/m3')
        print(fig,['NO3_2sig_',nn(i,:),'_st._seasonal_climate.png'],'-dpng')
    end
    
        %% NH4
    if i==1
        fig=figure;hold on;
        plot(2:6:12,koem_2sig.eachst_clim_gy_nh4_ft_1(i,2:6:12)./N_MW,'ro'); 
        plot(2:6:12,koem_2sig.eachst_clim_gy_nh4_ft_2(i,2:6:12)./N_MW,'go');
        plot(2:6:12,koem_2sig.eachst_clim_gy_nh4_ft_3(i,2:6:12)./N_MW,'bo'); 
        plot(2:6:12,koem_2sig.eachst_clim_gy_nh4_ft_1(i,2:6:12)./N_MW,'r'); 
        plot(2:6:12,koem_2sig.eachst_clim_gy_nh4_ft_2(i,2:6:12)./N_MW,'g');
        plot(2:6:12,koem_2sig.eachst_clim_gy_nh4_ft_3(i,2:6:12)./N_MW,'b'); 
        grid on;
            errorbar(1:12,koem_2sig.eachst_clim_gy_nh4_ft_1(i,:)./N_MW,koem_2sig.eachst_std_gy_nh4_ft_1(i,:)./N_MW,'r'); 
            errorbar(1:12,koem_2sig.eachst_clim_gy_nh4_ft_2(i,:)./N_MW,koem_2sig.eachst_std_gy_nh4_ft_2(i,:)./N_MW,'g');
            errorbar(1:12,koem_2sig.eachst_clim_gy_nh4_ft_3(i,:)./N_MW,koem_2sig.eachst_std_gy_nh4_ft_3(i,:)./N_MW,'b'); 
        alpha(0.3); 
        xlim([1 12]); ylim([0 12])
        title(['NH4- 2sig - ',nn(i,:),' st. seasonal climate']);
        xlabel('month')
        ylabel('mmol/m3')
        print(fig,['NH4_2sig_',nn(i,:),'_st._seasonal_climate.png'],'-dpng')
    else
        fig=figure;hold on;
        plot(2:3:12,koem_2sig.eachst_clim_gy_nh4_ft_1(i,2:3:12)./N_MW,'ro'); 
        plot(2:3:12,koem_2sig.eachst_clim_gy_nh4_ft_2(i,2:3:12)./N_MW,'go');
        plot(2:3:12,koem_2sig.eachst_clim_gy_nh4_ft_3(i,2:3:12)./N_MW,'bo'); 
        plot(2:3:12,koem_2sig.eachst_clim_gy_nh4_ft_1(i,2:3:12)./N_MW,'r'); 
        plot(2:3:12,koem_2sig.eachst_clim_gy_nh4_ft_2(i,2:3:12)./N_MW,'g');
        plot(2:3:12,koem_2sig.eachst_clim_gy_nh4_ft_3(i,2:3:12)./N_MW,'b'); 
        grid on;
            errorbar(1:12,koem_2sig.eachst_clim_gy_nh4_ft_1(i,:)./N_MW,koem_2sig.eachst_std_gy_nh4_ft_1(i,:)./N_MW,'r'); 
            errorbar(1:12,koem_2sig.eachst_clim_gy_nh4_ft_2(i,:)./N_MW,koem_2sig.eachst_std_gy_nh4_ft_2(i,:)./N_MW,'g');
            errorbar(1:12,koem_2sig.eachst_clim_gy_nh4_ft_3(i,:)./N_MW,koem_2sig.eachst_std_gy_nh4_ft_3(i,:)./N_MW,'b'); 
        alpha(0.3)
        xlim([1 12]); ylim([0 12])
        title(['NH4- 2sig - ',nn(i,:),' st. seasonal climate']);
        xlabel('month')
        ylabel('mmol/m3')
        print(fig,['NH4_2sig_',nn(i,:),'_st._seasonal_climate.png'],'-dpng')
    end
    
        %% PO4   
    if i==1
        fig=figure;hold on;
        plot(2:6:12,koem_2sig.eachst_clim_gy_po4_ft_1(i,2:6:12)./P_MW,'ro'); 
        plot(2:6:12,koem_2sig.eachst_clim_gy_po4_ft_2(i,2:6:12)./P_MW,'go');
        plot(2:6:12,koem_2sig.eachst_clim_gy_po4_ft_3(i,2:6:12)./P_MW,'bo'); 
        plot(2:6:12,koem_2sig.eachst_clim_gy_po4_ft_1(i,2:6:12)./P_MW,'r'); 
        plot(2:6:12,koem_2sig.eachst_clim_gy_po4_ft_2(i,2:6:12)./P_MW,'g');
        plot(2:6:12,koem_2sig.eachst_clim_gy_po4_ft_3(i,2:6:12)./P_MW,'b'); 
        grid on;
            errorbar(1:12,koem_2sig.eachst_clim_gy_po4_ft_1(i,:)./P_MW,koem_2sig.eachst_std_gy_po4_ft_1(i,:)./P_MW,'r'); 
            errorbar(1:12,koem_2sig.eachst_clim_gy_po4_ft_2(i,:)./P_MW,koem_2sig.eachst_std_gy_po4_ft_2(i,:)./P_MW,'g');
            errorbar(1:12,koem_2sig.eachst_clim_gy_po4_ft_3(i,:)./P_MW,koem_2sig.eachst_std_gy_po4_ft_3(i,:)./P_MW,'b'); 
        alpha(0.3); 
        xlim([1 12]); ylim([0 180./P_MW])
        title(['PO4- 2sig - ',nn(i,:),' st. seasonal climate']);
        xlabel('month')
        ylabel('mmol/m3')
        print(fig,['PO4_2sig_',nn(i,:),'_st._seasonal_climate.png'],'-dpng')
    else
        fig=figure;hold on;
        plot(2:3:12,koem_2sig.eachst_clim_gy_po4_ft_1(i,2:3:12)./P_MW,'ro'); 
        plot(2:3:12,koem_2sig.eachst_clim_gy_po4_ft_2(i,2:3:12)./P_MW,'go');
        plot(2:3:12,koem_2sig.eachst_clim_gy_po4_ft_3(i,2:3:12)./P_MW,'bo'); 
        plot(2:3:12,koem_2sig.eachst_clim_gy_po4_ft_1(i,2:3:12)./P_MW,'r'); 
        plot(2:3:12,koem_2sig.eachst_clim_gy_po4_ft_2(i,2:3:12)./P_MW,'g');
        plot(2:3:12,koem_2sig.eachst_clim_gy_po4_ft_3(i,2:3:12)./P_MW,'b'); 
        grid on;
            errorbar(1:12,koem_2sig.eachst_clim_gy_po4_ft_1(i,:)./P_MW,koem_2sig.eachst_std_gy_po4_ft_1(i,:)./P_MW,'r'); 
            errorbar(1:12,koem_2sig.eachst_clim_gy_po4_ft_2(i,:)./P_MW,koem_2sig.eachst_std_gy_po4_ft_2(i,:)./P_MW,'g');
            errorbar(1:12,koem_2sig.eachst_clim_gy_po4_ft_3(i,:)./P_MW,koem_2sig.eachst_std_gy_po4_ft_3(i,:)./P_MW,'b'); 
        alpha(0.3)
        xlim([1 12]); ylim([0 180./P_MW])
        title(['PO4- 2sig - ',nn(i,:),' st. seasonal climate']);
        xlabel('month')
        ylabel('mmol/m3')
        print(fig,['PO4_2sig_',nn(i,:),'_st._seasonal_climate.png'],'-dpng')
    end
    
        %% Chl    
    if i==1
        fig=figure;hold on;
        plot(2:6:12,koem_2sig.eachst_clim_gy_chl_ft_1(i,2:6:12),'ro'); 
        plot(2:6:12,koem_2sig.eachst_clim_gy_chl_ft_2(i,2:6:12),'go');
        plot(2:6:12,koem_2sig.eachst_clim_gy_chl_ft_3(i,2:6:12),'bo'); 
        plot(2:6:12,koem_2sig.eachst_clim_gy_chl_ft_1(i,2:6:12),'r'); 
        plot(2:6:12,koem_2sig.eachst_clim_gy_chl_ft_2(i,2:6:12),'g');
        plot(2:6:12,koem_2sig.eachst_clim_gy_chl_ft_3(i,2:6:12),'b'); 
        grid on;
            errorbar(1:12,koem_2sig.eachst_clim_gy_chl_ft_1(i,:),koem_2sig.eachst_std_gy_chl_ft_1(i,:),'r'); 
            errorbar(1:12,koem_2sig.eachst_clim_gy_chl_ft_2(i,:),koem_2sig.eachst_std_gy_chl_ft_2(i,:),'g');
            errorbar(1:12,koem_2sig.eachst_clim_gy_chl_ft_3(i,:),koem_2sig.eachst_std_gy_chl_ft_3(i,:),'b'); 
        alpha(0.3); 
        xlim([1 12]); ylim([0 14])
        title(['chl- 2sig - ',nn(i,:),' st. seasonal climate']);
        xlabel('month')
        ylabel('ug/L')
        print(fig,['chl_2sig_',nn(i,:),'_st._seasonal_climate.png'],'-dpng')
    else
        fig=figure;hold on;
        plot(2:3:12,koem_2sig.eachst_clim_gy_chl_ft_1(i,2:3:12),'ro'); 
        plot(2:3:12,koem_2sig.eachst_clim_gy_chl_ft_2(i,2:3:12),'go');
        plot(2:3:12,koem_2sig.eachst_clim_gy_chl_ft_3(i,2:3:12),'bo'); 
        plot(2:3:12,koem_2sig.eachst_clim_gy_chl_ft_1(i,2:3:12),'r'); 
        plot(2:3:12,koem_2sig.eachst_clim_gy_chl_ft_2(i,2:3:12),'g');
        plot(2:3:12,koem_2sig.eachst_clim_gy_chl_ft_3(i,2:3:12),'b'); 
        grid on;
            errorbar(1:12,koem_2sig.eachst_clim_gy_chl_ft_1(i,:),koem_2sig.eachst_std_gy_chl_ft_1(i,:),'r'); 
            errorbar(1:12,koem_2sig.eachst_clim_gy_chl_ft_2(i,:),koem_2sig.eachst_std_gy_chl_ft_2(i,:),'g');
            errorbar(1:12,koem_2sig.eachst_clim_gy_chl_ft_3(i,:),koem_2sig.eachst_std_gy_chl_ft_3(i,:),'b'); 
        alpha(0.3)
        xlim([1 12]); ylim([0 14])
        title(['chl- 2sig - ',nn(i,:),' st. seasonal climate']);
        xlabel('month')
        ylabel('ug/L')
        print(fig,['chl_2sig_',nn(i,:),'_st._seasonal_climate.png'],'-dpng')
    end
    
%% Bottom
    %% NO3
    if i==1
        fig=figure;hold on;
        plot(2:6:12,koem_2sig.eachst_clim_gy_no3_b_ft_1(i,2:6:12)./N_MW,'ro'); 
        plot(2:6:12,koem_2sig.eachst_clim_gy_no3_b_ft_2(i,2:6:12)./N_MW,'go');
        plot(2:6:12,koem_2sig.eachst_clim_gy_no3_b_ft_3(i,2:6:12)./N_MW,'bo'); 
        plot(2:6:12,koem_2sig.eachst_clim_gy_no3_b_ft_1(i,2:6:12)./N_MW,'r'); 
        plot(2:6:12,koem_2sig.eachst_clim_gy_no3_b_ft_2(i,2:6:12)./N_MW,'g');
        plot(2:6:12,koem_2sig.eachst_clim_gy_no3_b_ft_3(i,2:6:12)./N_MW,'b'); 
        grid on;
            errorbar(1:12,koem_2sig.eachst_clim_gy_no3_b_ft_1(i,:)./N_MW,koem_2sig.eachst_std_gy_no3_b_ft_1(i,:)./N_MW,'r'); 
            errorbar(1:12,koem_2sig.eachst_clim_gy_no3_b_ft_2(i,:)./N_MW,koem_2sig.eachst_std_gy_no3_b_ft_2(i,:)./N_MW,'g');
            errorbar(1:12,koem_2sig.eachst_clim_gy_no3_b_ft_3(i,:)./N_MW,koem_2sig.eachst_std_gy_no3_b_ft_3(i,:)./N_MW,'b'); 
        alpha(0.3); 
        xlim([1 12]); ylim([0 24])
        title(['NO3- 2sig - ',nn(i,:),' st. seasonal climate bot']);
        xlabel('month')
        ylabel('mmol/m3')
        print(fig,['NO3_2sig_',nn(i,:),'_st._seasonal_climate_bot.png'],'-dpng')
        
    else
        fig=figure;hold on;
        plot(2:3:12,koem_2sig.eachst_clim_gy_no3_b_ft_1(i,2:3:12)./N_MW,'ro'); 
        plot(2:3:12,koem_2sig.eachst_clim_gy_no3_b_ft_2(i,2:3:12)./N_MW,'go');
        plot(2:3:12,koem_2sig.eachst_clim_gy_no3_b_ft_3(i,2:3:12)./N_MW,'bo'); 
        plot(2:3:12,koem_2sig.eachst_clim_gy_no3_b_ft_1(i,2:3:12)./N_MW,'r'); 
        plot(2:3:12,koem_2sig.eachst_clim_gy_no3_b_ft_2(i,2:3:12)./N_MW,'g');
        plot(2:3:12,koem_2sig.eachst_clim_gy_no3_b_ft_3(i,2:3:12)./N_MW,'b'); 
        grid on;
            errorbar(1:12,koem_2sig.eachst_clim_gy_no3_b_ft_1(i,:)./N_MW, koem_2sig.eachst_std_gy_no3_b_ft_1(i,:)./N_MW,'r'); 
            errorbar(1:12,koem_2sig.eachst_clim_gy_no3_b_ft_2(i,:)./N_MW, koem_2sig.eachst_std_gy_no3_b_ft_2(i,:)./N_MW,'g');
            errorbar(1:12,koem_2sig.eachst_clim_gy_no3_b_ft_3(i,:)./N_MW, koem_2sig.eachst_std_gy_no3_b_ft_3(i,:)./N_MW,'b'); 
        alpha(0.3)
        xlim([1 12]); ylim([0 24])
        title(['NO3- 2sig - ',nn(i,:),' st. seasonal climate bot']);
        xlabel('month')
        ylabel('mmol/m3')
        print(fig,['NO3_2sig_',nn(i,:),'_st._seasonal_climate_bot.png'],'-dpng')
    end
    
        %% NH4
    if i==1
        fig=figure;hold on;
        plot(2:6:12,koem_2sig.eachst_clim_gy_nh4_b_ft_1(i,2:6:12)./N_MW,'ro'); 
        plot(2:6:12,koem_2sig.eachst_clim_gy_nh4_b_ft_2(i,2:6:12)./N_MW,'go');
        plot(2:6:12,koem_2sig.eachst_clim_gy_nh4_b_ft_3(i,2:6:12)./N_MW,'bo'); 
        plot(2:6:12,koem_2sig.eachst_clim_gy_nh4_b_ft_1(i,2:6:12)./N_MW,'r'); 
        plot(2:6:12,koem_2sig.eachst_clim_gy_nh4_b_ft_2(i,2:6:12)./N_MW,'g');
        plot(2:6:12,koem_2sig.eachst_clim_gy_nh4_b_ft_3(i,2:6:12)./N_MW,'b'); 
        grid on;
            errorbar(1:12,koem_2sig.eachst_clim_gy_nh4_b_ft_1(i,:)./N_MW,koem_2sig.eachst_std_gy_nh4_b_ft_1(i,:)./N_MW,'r'); 
            errorbar(1:12,koem_2sig.eachst_clim_gy_nh4_b_ft_2(i,:)./N_MW,koem_2sig.eachst_std_gy_nh4_b_ft_2(i,:)./N_MW,'g');
            errorbar(1:12,koem_2sig.eachst_clim_gy_nh4_b_ft_3(i,:)./N_MW,koem_2sig.eachst_std_gy_nh4_b_ft_3(i,:)./N_MW,'b'); 
        alpha(0.3); 
        xlim([1 12]); ylim([0 12])
        title(['NH4- 2sig - ',nn(i,:),' st. seasonal climate bot']);
        xlabel('month')
        ylabel('mmol/m3')
        print(fig,['NH4_2sig_',nn(i,:),'_st._seasonal_climate_bot.png'],'-dpng')
    else
        fig=figure;hold on;
        plot(2:3:12,koem_2sig.eachst_clim_gy_nh4_b_ft_1(i,2:3:12)./N_MW,'ro'); 
        plot(2:3:12,koem_2sig.eachst_clim_gy_nh4_b_ft_2(i,2:3:12)./N_MW,'go');
        plot(2:3:12,koem_2sig.eachst_clim_gy_nh4_b_ft_3(i,2:3:12)./N_MW,'bo'); 
        plot(2:3:12,koem_2sig.eachst_clim_gy_nh4_b_ft_1(i,2:3:12)./N_MW,'r'); 
        plot(2:3:12,koem_2sig.eachst_clim_gy_nh4_b_ft_2(i,2:3:12)./N_MW,'g');
        plot(2:3:12,koem_2sig.eachst_clim_gy_nh4_b_ft_3(i,2:3:12)./N_MW,'b'); 
        grid on;
            errorbar(1:12,koem_2sig.eachst_clim_gy_nh4_b_ft_1(i,:)./N_MW,koem_2sig.eachst_std_gy_nh4_b_ft_1(i,:)./N_MW,'r'); 
            errorbar(1:12,koem_2sig.eachst_clim_gy_nh4_b_ft_2(i,:)./N_MW,koem_2sig.eachst_std_gy_nh4_b_ft_2(i,:)./N_MW,'g');
            errorbar(1:12,koem_2sig.eachst_clim_gy_nh4_b_ft_3(i,:)./N_MW,koem_2sig.eachst_std_gy_nh4_b_ft_3(i,:)./N_MW,'b'); 
        alpha(0.3)
        xlim([1 12]); ylim([0 12])
        title(['NH4- 2sig - ',nn(i,:),' st. seasonal climate bot']);
        xlabel('month')
        ylabel('mmol/m3')
        print(fig,['NH4_2sig_',nn(i,:),'_st._seasonal_climate_bot.png'],'-dpng')
    end
    
        %% PO4   
    if i==1
        fig=figure;hold on;
        plot(2:6:12,koem_2sig.eachst_clim_gy_po4_b_ft_1(i,2:6:12)./P_MW,'ro'); 
        plot(2:6:12,koem_2sig.eachst_clim_gy_po4_b_ft_2(i,2:6:12)./P_MW,'go');
        plot(2:6:12,koem_2sig.eachst_clim_gy_po4_b_ft_3(i,2:6:12)./P_MW,'bo'); 
        plot(2:6:12,koem_2sig.eachst_clim_gy_po4_b_ft_1(i,2:6:12)./P_MW,'r'); 
        plot(2:6:12,koem_2sig.eachst_clim_gy_po4_b_ft_2(i,2:6:12)./P_MW,'g');
        plot(2:6:12,koem_2sig.eachst_clim_gy_po4_b_ft_3(i,2:6:12)./P_MW,'b'); 
        grid on;
            errorbar(1:12,koem_2sig.eachst_clim_gy_po4_b_ft_1(i,:)./P_MW,koem_2sig.eachst_std_gy_po4_b_ft_1(i,:)./P_MW,'r'); 
            errorbar(1:12,koem_2sig.eachst_clim_gy_po4_b_ft_2(i,:)./P_MW,koem_2sig.eachst_std_gy_po4_b_ft_2(i,:)./P_MW,'g');
            errorbar(1:12,koem_2sig.eachst_clim_gy_po4_b_ft_3(i,:)./P_MW,koem_2sig.eachst_std_gy_po4_b_ft_3(i,:)./P_MW,'b'); 
        alpha(0.3); 
        xlim([1 12]); ylim([0 180./P_MW])
        title(['PO4- 2sig - ',nn(i,:),' st. seasonal climate bot']);
        xlabel('month')
        ylabel('mmol/m3')
        print(fig,['PO4_2sig_',nn(i,:),'_st._seasonal_climate_bot.png'],'-dpng')
    else
        fig=figure;hold on;
        plot(2:3:12,koem_2sig.eachst_clim_gy_po4_b_ft_1(i,2:3:12)./P_MW,'ro'); 
        plot(2:3:12,koem_2sig.eachst_clim_gy_po4_b_ft_2(i,2:3:12)./P_MW,'go');
        plot(2:3:12,koem_2sig.eachst_clim_gy_po4_b_ft_3(i,2:3:12)./P_MW,'bo'); 
        plot(2:3:12,koem_2sig.eachst_clim_gy_po4_b_ft_1(i,2:3:12)./P_MW,'r'); 
        plot(2:3:12,koem_2sig.eachst_clim_gy_po4_b_ft_2(i,2:3:12)./P_MW,'g');
        plot(2:3:12,koem_2sig.eachst_clim_gy_po4_b_ft_3(i,2:3:12)./P_MW,'b'); 
        grid on;
            errorbar(1:12,koem_2sig.eachst_clim_gy_po4_b_ft_1(i,:)./P_MW,koem_2sig.eachst_std_gy_po4_b_ft_1(i,:)./P_MW,'r'); 
            errorbar(1:12,koem_2sig.eachst_clim_gy_po4_b_ft_2(i,:)./P_MW,koem_2sig.eachst_std_gy_po4_b_ft_2(i,:)./P_MW,'g');
            errorbar(1:12,koem_2sig.eachst_clim_gy_po4_b_ft_3(i,:)./P_MW,koem_2sig.eachst_std_gy_po4_b_ft_3(i,:)./P_MW,'b'); 
        alpha(0.3)
        xlim([1 12]); ylim([0 180./P_MW])
        title(['PO4- 2sig - ',nn(i,:),' st. seasonal climate bot']);
        xlabel('month')
        ylabel('mmol/m3')
        print(fig,['PO4_2sig_',nn(i,:),'_st._seasonal_climate_bot.png'],'-dpng')
    end
    
        %% Chl    
    if i==1
        fig=figure;hold on;
        plot(2:6:12,koem_2sig.eachst_clim_gy_chl_b_ft_1(i,2:6:12),'ro'); 
        plot(2:6:12,koem_2sig.eachst_clim_gy_chl_b_ft_2(i,2:6:12),'go');
        plot(2:6:12,koem_2sig.eachst_clim_gy_chl_b_ft_3(i,2:6:12),'bo'); 
        plot(2:6:12,koem_2sig.eachst_clim_gy_chl_b_ft_1(i,2:6:12),'r'); 
        plot(2:6:12,koem_2sig.eachst_clim_gy_chl_b_ft_2(i,2:6:12),'g');
        plot(2:6:12,koem_2sig.eachst_clim_gy_chl_b_ft_3(i,2:6:12),'b'); 
        grid on;
            errorbar(1:12,koem_2sig.eachst_clim_gy_chl_b_ft_1(i,:),koem_2sig.eachst_std_gy_chl_b_ft_1(i,:),'r'); 
            errorbar(1:12,koem_2sig.eachst_clim_gy_chl_b_ft_2(i,:),koem_2sig.eachst_std_gy_chl_b_ft_2(i,:),'g');
            errorbar(1:12,koem_2sig.eachst_clim_gy_chl_b_ft_3(i,:),koem_2sig.eachst_std_gy_chl_b_ft_3(i,:),'b'); 
        alpha(0.3); 
        xlim([1 12]); ylim([0 14])
        title(['chl- 2sig - ',nn(i,:),' st. seasonal climate bot']);
        xlabel('month')
        ylabel('ug/L')
        print(fig,['chl_2sig_',nn(i,:),'_st._seasonal_climate_bot.png'],'-dpng')
    else
        fig=figure;hold on;
        plot(2:3:12,koem_2sig.eachst_clim_gy_chl_b_ft_1(i,2:3:12),'ro'); 
        plot(2:3:12,koem_2sig.eachst_clim_gy_chl_b_ft_2(i,2:3:12),'go');
        plot(2:3:12,koem_2sig.eachst_clim_gy_chl_b_ft_3(i,2:3:12),'bo'); 
        plot(2:3:12,koem_2sig.eachst_clim_gy_chl_b_ft_1(i,2:3:12),'r'); 
        plot(2:3:12,koem_2sig.eachst_clim_gy_chl_b_ft_2(i,2:3:12),'g');
        plot(2:3:12,koem_2sig.eachst_clim_gy_chl_b_ft_3(i,2:3:12),'b'); 
        grid on;
            errorbar(1:12,koem_2sig.eachst_clim_gy_chl_b_ft_1(i,:),koem_2sig.eachst_std_gy_chl_b_ft_1(i,:),'r'); 
            errorbar(1:12,koem_2sig.eachst_clim_gy_chl_b_ft_2(i,:),koem_2sig.eachst_std_gy_chl_b_ft_2(i,:),'g');
            errorbar(1:12,koem_2sig.eachst_clim_gy_chl_b_ft_3(i,:),koem_2sig.eachst_std_gy_chl_b_ft_3(i,:),'b'); 
        alpha(0.3)
        xlim([1 12]); ylim([0 14])
        title(['chl- 2sig - ',nn(i,:),' st. seasonal climate bot']);
        xlabel('month')
        ylabel('ug/L')
        print(fig,['chl_2sig_',nn(i,:),'_st._seasonal_climate_bot.png'],'-dpng')
    end
end
close all