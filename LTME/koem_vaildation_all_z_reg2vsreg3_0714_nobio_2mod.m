close all; clc; clear;
% clim1=load('D:\장기생태\Dynamic\06_river\koem_climate_to06to15_v5_sig_gy_9points.mat'); %default
clim1=load('D:\장기생태\Dynamic\06_river\koem_climate_to06to15_v5_sig_gy_9points_no_regime_cut.mat'); %no cut regime
% clim1=load('D:\장기생태\Dynamic\06_river\koem_climate_to06to15_v5_sig_gy_8points_no_regime_cut.mat'); %no cut regime 8 points

yoonja=load('yoonjakangs_koem_data_monthly_v3_16points_2021.mat','mon_d','time');
sp_9p_st = [1:6, 14:16];
mon_d=yoonja.mon_d(sp_9p_st,:);
time_hr=yoonja.time(sp_9p_st,:);

plt_time = zeros(1,12);
plt_time(1,[2,5,8,11]) = [31, 120, 212, 304]; 

k=0
for i = 1997:2021
    for j = 1:12
        k=k+1;
        for l = 1:size(mon_d,1)
            if isnan(mon_d(l,k))==0
               ref_date{l,k} = plt_time(1,j) + mon_d(l,k);
               if leapyear(i)==1 & j > 2 & ref_date{l,k} >= 59
                   disp(['leapyear ',num2str(i),'-',num2str(j)])
                   ref_date{l,k} = ref_date{l,k} + leapyear(i);
               end
            else
                ref_date{l,k}=NaN;
            end
        end
    end
end

k=0
for i = 1997:2021
    for j = 1:12
        k=k+1;
        ref_mon{k,1} = [num2str(i) '-' num2str(j,'%02d')];
    end
end

cd J:\장기생태_2021\Dynamic\result\3regime
% regime1_kdc
% regime1_kdc_zerosw=load('kd_double_2reg.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime1_kdc_zerosw=load('merg_part_t4_2reg_mix.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime1_kdc_zerosw=load('2reg_kd_t4_mix_2.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime1_kdc_zerosw=load('2reg_kd_t4_3.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime1_kdc_zerosw=load('2reg_t4_1_mix_bio_part.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime1_kdc_zerosw=load('2reg_kd_t4_1_mix_part.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime1_kdc_zerosw=load('2reg_kd_t4_cla.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime1_kdc_zerosw=load('2reg_kd_t4_cir.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime1_kdc_zerosw=load('2reg_cir_2yr.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime1_kdc_zerosw=load('kd_t4w_2reg.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime1_kdc_zerosw=load('part_kd_t6_2reg.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime1_riv_chl=load('1regime_part_res_riv_chl.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_kdc=
% regime2_kdc_zerosw=load('kd_double_3reg.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_kdc_zerosw=load('kd_t4_2reg.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_kdc_zerosw=load('2reg_kd_t4_1_mix_part.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_kdc_zerosw=load('3reg_cir_2yr.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_kdc_zerosw=load('3reg_t4_1_mix_bio_part.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% [num2str(refyear),'reg_t4_1_mix_bio_part'];
% regime2_kdc_zerosw=load('3reg_kd_t4_1_mix_part.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_kdc_zerosw=load('2reg_kd_t4_1_mix_part.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_kdc_zerosw=load('kd_t4w_3reg.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_kdc_zerosw=load('part_kd_t6_3reg.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_noriv_chl=load('2regime_part_res_2regime_noriv_chl.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_riv_chl=load('2regime_part_res_2regime_riv_chl.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');

regime1_kdc_zerosw=load(['J:\장기생태_2021\Dynamic\result\3regime\','kd_t9_2reg.mat'],'obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_kdc_zerosw=load(['J:\장기생태_2021\Dynamic\result\3regime\','kd_t9_3reg.mat'],'obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
regime2_kdc_zerosw=load(['J:\장기생태_2021\Dynamic\result\3regime\','2reg_kd_t4_1_mix_part.mat'],'obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');

% plot_dir = 'J:\장기생태_2021\Dynamic\result\3regime\kdc_nosew\kd_double\2&3regime\4sigma'
% plot_dir = 'J:\장기생태_2021\Dynamic\result\3regime\kdc_nosew\t4_case\2&3regime\4sigma'
% plot_dir = 'J:\장기생태_2021\Dynamic\result\3regime\kdc_nosew\t4w_case\2&3regime\4sigma'
% plot_dir = 'J:\장기생태_2021\Dynamic\result\3regime\kdc_nosew\t4_case\2&3regime\'

sigma = 3;

% plot_dir = 'J:\장기생태_2021\Dynamic\result\3regime\kdc_nosew\kd_double\2&3regime\nocutregime\3sig'
% plot_dir = 'J:\장기생태_2021\Dynamic\result\3regime\kdc_nosew\t4_case\2&3regime\nocutregime\0815_mix_1'
% plot_dir = 'J:\장기생태_2021\Dynamic\result\3regime\kdc_nosew\t4_case\2&3regime\nocutregime\0815_mix_2'
% plot_dir = 'J:\장기생태_2021\Dynamic\result\3regime\kdc_nosew\t4_case\2&3regime\nocutregime\1027_cir\2yr'
% plot_dir ='J:\장기생태_2021\Dynamic\result\3regime\kdc_nosew\t4w_case\2&3regime\nocutregime\1sig'
% plot_dir ='J:\장기생태_2021\Dynamic\result\3regime\kdc_nosew\t6_case\2&3regime\nocutregime\3sig'
plot_dir = 'J:\장기생태_2021\Dynamic\result\3regime\kdc_nosew\t4_case\2&3regime\nocutregime\220702\in_t4vs30'


P_MW = 30.973762;
PO4_MW =94.971482;
N_MW = 14.006720;
NO3_MW = 62.005010;
NH4_MW = 18.038508;
DO_CONV=0.7.*44.661; %%% mg/L to mM/m^3;

xlim_set = [1 365];

%% plot obs 1st and 2nd at the same time which is time axis
plt_filter_1{1} = 32:45;
plt_filter_2{1} = 46:59;
plt_filter_1{2} = 121:135.5;
plt_filter_2{2} = 136:151;
plt_filter_1{3} = 213:227.5;
plt_filter_2{3} = 228:243;
plt_filter_1{4} = 305:319;
plt_filter_2{4} = 320:334;


if size(regime1_kdc_zerosw.sp_gy_temp_mid,3) <= size(regime2_kdc_zerosw.sp_gy_temp_mid,3)
   size_tt = size(regime1_kdc_zerosw.sp_gy_temp_mid,3);
else
   size_tt = size(regime2_kdc_zerosw.sp_gy_temp_mid,3);
end


clearvars temp_d sp_gy_*_mm
for i = 1:9 %st
    for j = 1: size_tt %time    
       
        clearvars temp_d
        temp_d=interp1(squeeze(regime1_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime1_kdc_zerosw.sp_gy_temp_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_temp_mm(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime1_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime1_kdc_zerosw.sp_gy_salt_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_salt_mm(i,:,j)=temp_d;
        
        
        % 2 regime 
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw.sp_gy_temp_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_temp_2reg_mm(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw.sp_gy_salt_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_salt_2reg_mm(i,:,j)=temp_d;
       
        
    end
end


                
% %% make same time (dd) with OBS
clearvars mod_*_pick mod_*_pick_2reg
               
                mod_temp_pick=NaN(9,32,size(ref_date,2));
                mod_salt_pick=NaN(9,32,size(ref_date,2));
                mod_temp_pick_2reg=NaN(9,32,size(ref_date,2));
                mod_salt_pick_2reg=NaN(9,32,size(ref_date,2));
                

for i = 1:9 %st
    for j = 1:size(ref_date,2) 
        if isnan(ref_date{i,j}) == 0 
            if j <= 228 %time ~ 2015.12             
                mod_temp_pick(i,:,j)=sp_gy_temp_mm(i,:,ref_date{i,j});
                mod_salt_pick(i,:,j)=sp_gy_salt_mm(i,:,ref_date{i,j});
            else % 2016.01 ~ 2019.12(2020.12) : ~276 (~288)
                % 2 regime                
                mod_temp_pick_2reg(i,:,j)=sp_gy_temp_2reg_mm(i,:,ref_date{i,j});
                mod_salt_pick_2reg(i,:,j)=sp_gy_salt_2reg_mm(i,:,ref_date{i,j});

            end
        else
             if j <= 228 %time ~ 2015.12                
                mod_temp_pick(i,:,j)=NaN;
                mod_salt_pick(i,:,j)=NaN;

            else % 2016.01 ~ 2021.12
                % 2 regime              
                mod_temp_pick_2reg(i,:,j)=NaN;
                mod_salt_pick_2reg(i,:,j)=NaN;

            end
        end
    end
end
cd(plot_dir)



     clearvars clim2*
clim2=clim1.obm_gy_temp_ft_2(:,sigma);
clim2_std = clim1.obs_std_gy_temp_ft_2(:,sigma);
clim2_2=clim1.obm_gy_temp_ft_3(:,sigma);
clim2_2_std = clim1.obs_std_gy_temp_ft_3(:,sigma);

    fig = figure; hold on;
            plot(squeeze(nanmean(sp_gy_temp_mm(:,3,:),1)),'g','linew',2);
            plot(squeeze(nanmean(sp_gy_temp_2reg_mm(:,3,:),1)),'b','linew',2);
%             clearvars *_bound_plt nonan_data_plt discon_p temp_d
%             temp_d=nanmean(mer_temp_1,1);
%             nonan_data_plt = find(isnan(temp_d)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             upper_bound_plt = temp_d + regm_std_temp_1;
%             lower_bound_plt = temp_d - regm_std_temp_1;
%             patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_temp_1; c1998.gy_temp_1; c1999.gy_temp_1; c2000.gy_temp_1; c2001.gy_temp_1; c2002.gy_temp_1; c2003.gy_temp_1;],1),'g','linew',2);
                         % 1regime
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2 + clim2_std;
            lower_bound_plt = clim2 - clim2_std;
            for i = 1:4
            patch([plt_filter_1{i} fliplr(plt_filter_1{i})], [upper_bound_plt(plt_filter_1{i})' fliplr(lower_bound_plt(plt_filter_1{i}))'],'g');
            plot(plt_filter_1{i},clim2(plt_filter_1{i}),'g-','linew',2);
            plot(plt_filter_1{i},upper_bound_plt(plt_filter_1{i}),'g-','linew',2);
            plot(plt_filter_1{i},lower_bound_plt(plt_filter_1{i}),'g-','linew',2);
            end             
 % 2regime           
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2_2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2_2 + clim2_2_std;
            lower_bound_plt = clim2_2 - clim2_2_std;
            for i = 1:4
            patch([plt_filter_2{i} fliplr(plt_filter_2{i})], [upper_bound_plt(plt_filter_2{i})' fliplr(lower_bound_plt(plt_filter_2{i}))'],'b');
            plot(plt_filter_2{i},clim2_2(plt_filter_2{i}),'b-','linew',2);
            plot(plt_filter_2{i},upper_bound_plt(plt_filter_2{i}),'b-','linew',2);
            plot(plt_filter_2{i},lower_bound_plt(plt_filter_2{i}),'b-','linew',2);
            end 
            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL temp']);
            xlabel('time','fontsize',13)
            ylabel('temp (^oC)','fontsize',13)
            grid on
            ylim([0 35])
            xlim(xlim_set)
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_temp'),'-dpng')
    
     
             clearvars clim2*
clim2=clim1.obm_gy_salt_ft_2(:,sigma);
clim2_std = clim1.obs_std_gy_salt_ft_1(:,sigma);
clim2_2=clim1.obm_gy_salt_ft_3(:,sigma);
clim2_2_std = clim1.obs_std_gy_salt_ft_2(:,sigma);   

    fig = figure; hold on; 
            plot(squeeze(nanmean(sp_gy_salt_mm(:,3,:),1)),'g','linew',2);
            plot(squeeze(nanmean(sp_gy_salt_2reg_mm(:,3,:),1)),'b','linew',2);
%             clearvars *_bound_plt nonan_data_plt discon_p salt_d
%             salt_d=nanmean(mer_salt_1,1);
%             nonan_data_plt = find(isnan(salt_d)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             upper_bound_plt = salt_d + regm_std_salt_1;
%             lower_bound_plt = salt_d - regm_std_salt_1;
%             patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_salt_1; c1998.gy_salt_1; c1999.gy_salt_1; c2000.gy_salt_1; c2001.gy_salt_1; c2002.gy_salt_1; c2003.gy_salt_1;],1),'r','linew',2);
                 % 1regime
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2 + clim2_std;
            lower_bound_plt = clim2 - clim2_std;
            for i = 1:4
            patch([plt_filter_1{i} fliplr(plt_filter_1{i})], [upper_bound_plt(plt_filter_1{i})' fliplr(lower_bound_plt(plt_filter_1{i}))'],'g');
            plot(plt_filter_1{i},clim2(plt_filter_1{i}),'g-','linew',2);
            plot(plt_filter_1{i},upper_bound_plt(plt_filter_1{i}),'g-','linew',2);
            plot(plt_filter_1{i},lower_bound_plt(plt_filter_1{i}),'g-','linew',2);
            end             
 % 2regime           
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2_2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2_2 + clim2_2_std;
            lower_bound_plt = clim2_2 - clim2_2_std;
            for i = 1:4
            patch([plt_filter_2{i} fliplr(plt_filter_2{i})], [upper_bound_plt(plt_filter_2{i})' fliplr(lower_bound_plt(plt_filter_2{i}))'],'b');
            plot(plt_filter_2{i},clim2_2(plt_filter_2{i}),'b-','linew',2);
            plot(plt_filter_2{i},upper_bound_plt(plt_filter_2{i}),'b-','linew',2);
            plot(plt_filter_2{i},lower_bound_plt(plt_filter_2{i}),'b-','linew',2);
            end 
            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL salt']);
            xlabel('time','fontsize',13)
            ylabel('salt (g/kg)','fontsize',13)
            grid on
            ylim([0 35])
            xlim(xlim_set)
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_salt'),'-dpng')
    
    
   
      


            clearvars clim2*
clim2=clim1.obm_gy_temp_b_ft_2(:,sigma);
clim2_std = clim1.obs_std_gy_temp_b_ft_2(:,sigma);
clim2_2=clim1.obm_gy_temp_b_ft_3(:,sigma);
clim2_2_std = clim1.obs_std_gy_temp_b_ft_3(:,sigma);  
                
    fig = figure; hold on;
            plot(nanmean(regime1_kdc_zerosw.gy_temp_b,1),'g','linew',2);
            plot(nanmean(regime2_kdc_zerosw.gy_temp_b,1),'b','linew',2);  
%             clearvars *_bound_plt nonan_data_plt discon_p temp_b_d
%             temp_b_d=nanmean(mer_temp_b_1,1);
%             nonan_data_plt = find(isnan(temp_b_d)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             upper_bound_plt = temp_b_d + regm_std_temp_b_1;
%             lower_bound_plt = temp_b_d - regm_std_temp_b_1;
%             patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    % plot(nanmean([c1997.gy_temp_b; c1998.gy_temp_b; c1999.gy_temp_b; c2000.gy_temp_b; c2001.gy_temp_b; c2002.gy_temp_b; c2003.gy_temp_b;],1),'g','linew',2);
                     % 1regime
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2 + clim2_std;
            lower_bound_plt = clim2 - clim2_std;
            for i = 1:4
            patch([plt_filter_1{i} fliplr(plt_filter_1{i})], [upper_bound_plt(plt_filter_1{i})' fliplr(lower_bound_plt(plt_filter_1{i}))'],'g');
            plot(plt_filter_1{i},clim2(plt_filter_1{i}),'g-','linew',2);
            plot(plt_filter_1{i},upper_bound_plt(plt_filter_1{i}),'g-','linew',2);
            plot(plt_filter_1{i},lower_bound_plt(plt_filter_1{i}),'g-','linew',2);
            end             
 % 2regime           
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2_2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2_2 + clim2_2_std;
            lower_bound_plt = clim2_2 - clim2_2_std;
            for i = 1:4
            patch([plt_filter_2{i} fliplr(plt_filter_2{i})], [upper_bound_plt(plt_filter_2{i})' fliplr(lower_bound_plt(plt_filter_2{i}))'],'b');
            plot(plt_filter_2{i},clim2_2(plt_filter_2{i}),'b-','linew',2);
            plot(plt_filter_2{i},upper_bound_plt(plt_filter_2{i}),'b-','linew',2);
            plot(plt_filter_2{i},lower_bound_plt(plt_filter_2{i}),'b-','linew',2);
            end 
            alpha(0.3) %transparency      
            title(['GY KOEM OBS vs. daily MODEL temp bot.']);
            xlabel('time','fontsize',13)
            ylabel('temp (^oC)','fontsize',13)
            grid on
            ylim([0 35])
            xlim(xlim_set)
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)

            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_temp_bot'),'-dpng')

            clearvars clim2*
clim2=clim1.obm_gy_salt_b_ft_2(:,sigma);
clim2_std = clim1.obs_std_gy_salt_b_ft_2(:,sigma);
clim2_2=clim1.obm_gy_salt_b_ft_3(:,sigma);
clim2_2_std = clim1.obs_std_gy_salt_b_ft_3(:,sigma);   
                
    fig = figure; hold on;
            plot(nanmean(regime1_kdc_zerosw.gy_salt_b,1),'g','linew',2);  
            plot(nanmean(regime2_kdc_zerosw.gy_salt_b,1),'b','linew',2);
%             clearvars *_bound_plt nonan_data_plt discon_p salt_b_d
%             salt_b_d=nanmean(mer_salt_b_1,1);
%             nonan_data_plt = find(isnan(salt_b_d)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             upper_bound_plt = salt_b_d + regm_std_salt_b_1;
%             lower_bound_plt = salt_b_d - regm_std_salt_b_1;
%             patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_salt_b; c1998.gy_salt_b; c1999.gy_salt_b; c2000.gy_salt_b; c2001.gy_salt_b; c2002.gy_salt_b; c2003.gy_salt_b;],1),'g','linew',2);
              % 1regime
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2 + clim2_std;
            lower_bound_plt = clim2 - clim2_std;
            for i = 1:4
            patch([plt_filter_1{i} fliplr(plt_filter_1{i})], [upper_bound_plt(plt_filter_1{i})' fliplr(lower_bound_plt(plt_filter_1{i}))'],'g');
            plot(plt_filter_1{i},clim2(plt_filter_1{i}),'g-','linew',2);
            plot(plt_filter_1{i},upper_bound_plt(plt_filter_1{i}),'g-','linew',2);
            plot(plt_filter_1{i},lower_bound_plt(plt_filter_1{i}),'g-','linew',2);
            end             
 % 2regime           
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2_2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2_2 + clim2_2_std;
            lower_bound_plt = clim2_2 - clim2_2_std;
            for i = 1:4
            patch([plt_filter_2{i} fliplr(plt_filter_2{i})], [upper_bound_plt(plt_filter_2{i})' fliplr(lower_bound_plt(plt_filter_2{i}))'],'b');
            plot(plt_filter_2{i},clim2_2(plt_filter_2{i}),'b-','linew',2);
            plot(plt_filter_2{i},upper_bound_plt(plt_filter_2{i}),'b-','linew',2);
            plot(plt_filter_2{i},lower_bound_plt(plt_filter_2{i}),'b-','linew',2);
            end 
            alpha(0.3) %transparency      
            title(['GY KOEM OBS vs. daily MODEL salt bot.']);
            xlabel('time','fontsize',13)
            ylabel('salt (g/kg)','fontsize',13)
            grid on
            ylim([0 35])
            xlim(xlim_set)
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_salt_bot'),'-dpng') 
            
            close all;

