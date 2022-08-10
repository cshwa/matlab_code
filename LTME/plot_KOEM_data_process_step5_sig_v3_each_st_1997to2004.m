close all; clc; clear;

cd D:\Àå±â»ıÅÂ\Dynamic\06_river
sig_extract = 1 % 1 =sigma extracted data, else not extracted data

% name_tag{sp_gy}
% ±¤¾çÇ×, ±¤¾ç4, ±¤¾ç3, ±¤¾ç2, ±¤¾ç1, ±¤¾ç5, ¿©¼ö2, ¿©¼ö3, ¿©¼ö1
% [res I]=sort([4,3,2,1,5]);
P_MW = 30.973762;
PO4_MW =94.971482;
N_MW = 14.006720;
NO3_MW = 62.005010;
NH4_MW = 18.038508;

% yoonja=load('yoonjakangs_koem_data_monthly.mat'); %
% yoonjakangs_koem_data_processing.m  - 9 points 
% yoonja=load('yoonjakangs_koem_data_monthly_v2_16points.mat'); % yoonjakangs_koem_data_processing_v2_13points.m
yoonja=load('D:\Àå±â»ıÅÂ\Dynamic\KOEM\yoonjakangs_koem_data_monthly_v3_16points_2021.mat'); % yoonjakangs_koem_data_processing_v2_13points.m

% input the sigma
sig = 3; %% sigma
% sig = 2; %% sigma

% % station name tag
% % (for preventing the miss picking using tag from 1 col's name on excel) 
[res I]=sort([1,5,4,3,2,6,8,9,7]);
% sorting confirm
nn=['±¤¾çÇ×'; '±¤¾ç4'; '±¤¾ç3'; '±¤¾ç2'; '±¤¾ç1'; '±¤¾ç5'; '¿©¼ö2'; '¿©¼ö3'; '¿©¼ö1'];

% combining the tag and outter point excluding
name_tag = nn(I,:); 
size_tag = length(name_tag);

%% pick matched name with tag
% %port
% for i = 1:length(name_tag)
%    if  sum(strcmp(name_tag{i}, txt_matc_p)) ~= 0
%        indx{i} = find([strcmp(name_tag{i}, txt_matc_p)] == 1)     
%    end
% end


%% make 1997 to 2019 'yymm' form
k=0
for i = 1997:2021
    for j = 1:12
        k=k+1;
        ref_date{k,1} = [num2str(i) '-' num2str(j,'%02d')];
    end
end

%% make 1997 to 2019 'yymm' form
for j = 1:12
 ref_date_mm{j,1} = [num2str(j,'%02d')];
end

% spatial region index
do_sur_clim = yoonja.do_sur;
no3_sur_clim = yoonja.no3_sur;
din_sur_clim = yoonja.DIN_sur;
temp_sur_clim = yoonja.temp_sur;
salt_sur_clim = yoonja.salt_sur;
chl_sur_clim = yoonja.chl_sur;
nh4_sur_clim = yoonja.nh4_sur;
ss_sur_clim = yoonja.ss_sur;
si_sur_clim = yoonja.si_sur;
secchi_sur_clim = yoonja.secchi_sur;
po4_sur_clim = yoonja.po4_sur;
% mon_d_obs = yoonja.mon_d;
% time_obs = yoonja.time;

po4_bot_clim = yoonja.po4_bot;        
do_bot_clim = yoonja.do_bot;
no3_bot_clim = yoonja.no3_bot;
din_bot_clim = yoonja.DIN_bot;
temp_bot_clim = yoonja.temp_bot;
salt_bot_clim = yoonja.salt_bot;
chl_bot_clim = yoonja.chl_bot;
nh4_bot_clim = yoonja.nh4_bot;
ss_bot_clim = yoonja.ss_bot;
si_bot_clim = yoonja.si_bot;

%% cut 3 regime
cut_1 = 96 % 2004-12
% cut_2 = 228 % 2015-12

              
for varlist = { 'do','chl','ss','si','po4','no3','nh4','din','temp','salt','secchi'}
    clearvars varname 
    varname = char(varlist);
    eval([varname,'_sur_cut_1',' = ',varname,'_sur_clim(:,1:cut_1);']);
%     eval([varname,'_sur_cut_2',' = ',varname,'_sur_clim(:,cut_1+1:cut_2);']);
%     eval([varname,'_sur_cut_3',' = ',varname,'_sur_clim(:,cut_2+1:end);']);

    if strcmp(varname,'secchi') == 0 % there are no secchi_bot
        eval([varname,'_bot_cut_1',' = ',varname,'_bot_clim(:,1:cut_1);']);
%         eval([varname,'_bot_cut_2',' = ',varname,'_bot_clim(:,cut_1+1:cut_2);']);
%         eval([varname,'_bot_cut_3',' = ',varname,'_bot_clim(:,cut_2+1:end);']);
    end
end

if sig_extract == 1
%% extract over 3sig
sp_gy=size(yoonja.do_bot,1); %% num of spatial point

for varlist = { 'do','chl','ss','si','po4','no3','nh4','din','temp','salt','secchi'}
        clearvars varname 
        varname = char(varlist);
    for i = 1:1 % num of regime
        for j=1:sp_gy %% num of spatial point
            for sig=1:6 %% sigma
                clearvars data data_ref
                eval(['data_ref = ', varname,'_sur_cut_',num2str(i),'(j,:);']);
                eval(['data = ', varname,'_sur_cut_',num2str(i),'(j,:);']);
                data(data_ref > nanmean(data_ref) + sig*nanstd(data_ref)) =NaN;
                data(data_ref < nanmean(data_ref) - sig*nanstd(data_ref)) =NaN;
                eval(['mtx_regime_',varname,'_',num2str(i),'_s(j,:,sig) = data;']);  % extracted data OUT [sp_gy, T, sig]
                disp(['mtx_regime_',varname,'_',num2str(i),'_s'])
           
                if strcmp(varname,'secchi') == 0 % there are no secchi_bot
                    clearvars data data_ref
                    eval(['data_ref = ', varname,'_bot_cut_',num2str(i),'(j,:);']);
                    eval(['data = ', varname,'_bot_cut_',num2str(i),'(j,:);']);
                    data(data_ref > nanmean(data_ref) + sig*nanstd(data_ref)) =NaN;
                    data(data_ref < nanmean(data_ref) - sig*nanstd(data_ref)) =NaN;
                    eval(['mtx_regime_',varname,'_',num2str(i),'_b_s(j,:,sig) = data;']);  % extracted data OUT [sp_gy, T, sig]
                   disp(['mtx_regime_',varname,'_',num2str(i),'_b_s'])
                end
            end
        end
    end
end




% save('koem_timeseires_monthly_gy_only_2sig_v2.mat');
% save(['koem_timeseires_monthly_gy_only_1to3sig_v2_each_st.mat']); % 9 points 
% save(['koem_monthly_gyonly_1to6sig_v2_each_st_16p_to06to15.mat']); % 16 points 
save(['koem_monthly_gyonly_1to6sig_v3_each_st_16p_1997to2004.mat']); % 16 points 

else
    %% no sigma extract
    sp_gy=size(yoonja.do_bot,1); %% num of spatial point

    for varlist = { 'do','chl','ss','si','po4','no3','nh4','din','temp','salt','secchi'}
            clearvars varname 
            varname = char(varlist);
        for i = 1:1 % num of regime
            for j=1:sp_gy %% num of spatial point
                    clearvars data data_ref
                    eval(['data_ref = ', varname,'_sur_cut_',num2str(i),'(j,:);']);
                    eval(['data = ', varname,'_sur_cut_',num2str(i),'(j,:);']);
                    eval(['mtx_regime_',varname,'_',num2str(i),'_s(j,:) = data;']);  % no extracted data OUT [sp_gy, T, sig]
                    disp(['mtx_regime_',varname,'_',num2str(i),'_s'])

                    if strcmp(varname,'secchi') == 0 % there are no secchi_bot
                        clearvars data data_ref
                        eval(['data_ref = ', varname,'_bot_cut_',num2str(i),'(j,:);']);
                        eval(['data = ', varname,'_bot_cut_',num2str(i),'(j,:);']);
                        eval(['mtx_regime_',varname,'_',num2str(i),'_b_s(j,:) = data;']);  % no extracted data OUT [sp_gy, T, sig]
                       disp(['mtx_regime_',varname,'_',num2str(i),'_b_s'])
                    end

            end
        end
    end

save(['koem_monthly_gyonly_nosig_v3_each_st_16p_1997to2004.mat']); % 16 points 
end
return

%% extraction confirm

c_c = {'r','g','b','k','m','c'};
figure;  hold on;
for i = 6:-1:1
    plot(squeeze(mtx_regime_nh4_1_s(3,:,i))./14,'o','color',c_c{i});
end

figure;  hold on;
for i = 6:-1:1
    plot(squeeze(mtx_regime_nh4_2_s(3,:,i))./14,'o','color',c_c{i});
end

figure;  hold on;
for i = 6:-1:1
    plot(squeeze(mtx_regime_po4_1_s(3,:,i))./30.9,'o','color',c_c{i});
end

figure;  hold on;
for i = 6:-1:1
    plot(squeeze(mtx_regime_po4_2_s(3,:,i))./30.9,'o','color',c_c{i});
end

figure;  hold on;
for i = 6:-1:1
    plot(squeeze(mtx_regime_po4_3_s(3,:,i))./30.9,'o','color',c_c{i});
end

%% make vaildation data for to06to15
% make eom_d
k=0;
for i = 2001:2001
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end
t_tick_pre=sum(eom_d,2);

for i = 1:size(eom_d,1)
    for j = 1:size(eom_d,2)
        eom_d_each(i,j) = sum(eom_d(i,1:j));
    end
end

if sig_extract == 1
sp_9p_st = [1:6, 14:16];    
for varlist = { 'do','chl','ss','si','po4','no3','nh4','din','temp','salt','secchi'}
    clearvars varname 
    varname = char(varlist);
    for itt = 1:3 %regime
        for sig = 1:6 %sigma
        for i = 1:12
             if i==1
                 clearvars t_eom_interval
                 t_eom_interval = 1:eom_d_each(1,i);
             else    
                 clearvars t_eom_interval
                 t_eom_interval = eom_d_each(1,i-1)+1:eom_d_each(1,i);
             end
            %% time & spatial std & mean
                clearvars len_temp data_temp
                eval(['data_temp = squeeze(mtx_regime_',varname,'_',num2str(itt),'_s(sp_9p_st,:,sig));']);  % extracted data OUT [sp_gy, T, sig]
                len_temp=size(data_temp(:,i:12:end),1)*size(data_temp(:,i:12:end),2);
                % time & spatial std 
                eval(['obs_std_gy_',varname,'_ft_',num2str(itt),'(t_eom_interval,sig)= nanstd(reshape(data_temp(:,i:12:end),1,len_temp));']);  % extracted data OUT [sp_gy, T, sig]
                % time std 
                eval(['spm_obs_std_gy_',varname,'_ft_',num2str(itt),'(t_eom_interval,sig)= nanstd(nanmean(data_temp(:,i:12:end),1));']);  % extracted data OUT [sp_gy, T, sig]
                % spatial std
                eval(['tm_obs_std_gy_',varname,'_ft_',num2str(itt),'(t_eom_interval,sig)= nanstd(nanmean(data_temp(:,i:12:end),2));']);  % extracted data OUT [sp_gy, T, sig]
                % mean
                eval(['obm_gy_',varname,'_ft_',num2str(itt),'(t_eom_interval,sig)= nanmean(nanmean(data_temp(:,i:12:end)));']);  % extracted data OUT [sp_gy, T, sig]
                if strcmp(varname,'secchi') == 0 % there are no secchi_bot
                    clearvars len_temp data_temp
                   eval(['data_temp = squeeze(mtx_regime_',varname,'_',num2str(itt),'_b_s(sp_9p_st,:,sig));']);  % extracted data OUT [sp_gy, T, sig]
                   len_temp=size(data_temp(:,i:12:end),1)*size(data_temp(:,i:12:end),2);
                   % time & spatial std 
                   eval(['obs_std_gy_',varname,'_b_ft_',num2str(itt),'(t_eom_interval,sig)= nanstd(reshape(data_temp(:,i:12:end),1,len_temp));']);  % extracted data OUT [sp_gy, T, sig]
                   % time std 
                   eval(['spm_obs_std_gy_',varname,'_b_ft_',num2str(itt),'(t_eom_interval,sig)= nanstd(nanmean(data_temp(:,i:12:end),1));']);  % extracted data OUT [sp_gy, T, sig]
                   % spatial std
                   eval(['tm_obs_std_gy_',varname,'_b_ft_',num2str(itt),'(t_eom_interval,sig)= nanstd(nanmean(data_temp(:,i:12:end),2));']);  % extracted data OUT [sp_gy, T, sig]
                   % mean
                 eval(['obm_gy_',varname,'_b_ft_',num2str(itt),'(t_eom_interval,sig)= nanmean(nanmean(data_temp(:,i:12:end)));']);  % extracted data OUT [sp_gy, T, sig]
                end
        end
        end
    end
end

sp_9p_st = [1:6, 14:16];
%% each std and mean cal.
for varlist = { 'do','chl','ss','si','po4','no3','nh4','din','temp','salt','secchi'}
    for itt= 1:3 %regime
        for sig = 1:6 %sigma
            for j=1:length(sp_9p_st) %% num of spatial point
                for i = 1:12
                     if i==1
                         clearvars t_eom_interval
                         t_eom_interval = 1:eom_d_each(1,i);
                     else    
                         clearvars t_eom_interval
                         t_eom_interval = eom_d_each(1,i-1)+1:eom_d_each(1,i);
                     end
                        %% time & spatial std & mean
                        clearvars len_temp data_temp
                        eval(['data_temp = squeeze(mtx_regime_',varname,'_',num2str(itt),'_s(sp_9p_st(j),:,sig));']);  % extracted data OUT [sp_gy, T, sig]
                        len_temp=size(data_temp(:,i:12:end),1)*size(data_temp(:,i:12:end),2);
                        % time & spatial std 
                        eval(['eachst_std_gy_',varname,'_ft_',num2str(itt),'(t_eom_interval,sig)= nanstd(data_temp(i:12:end));']);  % extracted data OUT [sp_gy, T, sig]
                        % mean
                        eval(['eachst_clim_gy_',varname,'_ft_',num2str(itt),'(t_eom_interval,sig)= nanmean(data_temp(i:12:end));']);  % extracted data OUT [sp_gy, T, sig]
                        if strcmp(varname,'secchi') == 0 % there are no secchi_bot
                           clearvars len_temp data_temp
                           eval(['data_temp = squeeze(mtx_regime_',varname,'_',num2str(itt),'_b_s(sp_9p_st(j),:,sig));']);  % extracted data OUT [sp_gy, T, sig]
                           % time & spatial std 
                           eval(['eachst_std_gy_',varname,'_b_ft_',num2str(itt),'(t_eom_interval,sig)= nanstd(data_temp(i:12:end));']);  % extracted data OUT [sp_gy, T, sig]
                           % mean
                           eval(['eachst_clim_gy_',varname,'_b_ft_',num2str(itt),'(t_eom_interval,sig)= nanmean(data_temp(i:12:end));']);  % extracted data OUT [sp_gy, T, sig]
                        end
                end
            end
        end
    end
end                
save(['koem_climate_to06to15_v5_sig_gy_9points.mat'],'obs*','obm*','spm_obs*','tm_obs*'); % full 9points data (extract_sigma_for_each_st)

else
sp_9p_st = [1:6, 14:16];
for varlist = { 'do','chl','ss','si','po4','no3','nh4','din','temp','salt','secchi'}
    clearvars varname 
    varname = char(varlist);
    for itt = 1:3 %regime
%         for sig = 1:3 %sigma
        for i = 1:12
             if i==1
                 clearvars t_eom_interval
                 t_eom_interval = 1:eom_d_each(1,i);
             else    
                 clearvars t_eom_interval
                 t_eom_interval = eom_d_each(1,i-1)+1:eom_d_each(1,i);
             end
            %% time & spatial std & mean
                clearvars len_temp data_temp
                eval(['data_temp = squeeze(mtx_regime_',varname,'_',num2str(itt),'_s(sp_9p_st,:));']);  % extracted data OUT [sp_gy, T, sig]
                len_temp=size(data_temp(:,i:12:end),1)*size(data_temp(:,i:12:end),2);
                % time & spatial std 
                eval(['obs_std_gy_',varname,'_ft_',num2str(itt),'(t_eom_interval)= nanstd(reshape(data_temp(:,i:12:end),1,len_temp));']);  % extracted data OUT [sp_gy, T, sig]
                % time std 
                eval(['spm_obs_std_gy_',varname,'_ft_',num2str(itt),'(t_eom_interval)= nanstd(nanmean(data_temp(:,i:12:end),1));']);  % extracted data OUT [sp_gy, T, sig]
                % spatial std
                eval(['tm_obs_std_gy_',varname,'_ft_',num2str(itt),'(t_eom_interval)= nanstd(nanmean(data_temp(:,i:12:end),2));']);  % extracted data OUT [sp_gy, T, sig]
                % mean
                eval(['obm_gy_',varname,'_ft_',num2str(itt),'(t_eom_interval)= nanmean(nanmean(data_temp(:,i:12:end)));']);  % extracted data OUT [sp_gy, T, sig]
                if strcmp(varname,'secchi') == 0 % there are no secchi_bot
                    clearvars len_temp data_temp
                   eval(['data_temp = squeeze(mtx_regime_',varname,'_',num2str(itt),'_b_s(sp_9p_st,:));']);  % extracted data OUT [sp_gy, T, sig]
                   len_temp=size(data_temp(:,i:12:end),1)*size(data_temp(:,i:12:end),2);
                   % time & spatial std 
                   eval(['obs_std_gy_',varname,'_b_ft_',num2str(itt),'(t_eom_interval)= nanstd(reshape(data_temp(:,i:12:end),1,len_temp));']);  % extracted data OUT [sp_gy, T, sig]
                   % time std 
                   eval(['spm_obs_std_gy_',varname,'_b_ft_',num2str(itt),'(t_eom_interval)= nanstd(nanmean(data_temp(:,i:12:end),1));']);  % extracted data OUT [sp_gy, T, sig]
                   % spatial std
                   eval(['tm_obs_std_gy_',varname,'_b_ft_',num2str(itt),'(t_eom_interval)= nanstd(nanmean(data_temp(:,i:12:end),2));']);  % extracted data OUT [sp_gy, T, sig]
                   % mean
                 eval(['obm_gy_',varname,'_b_ft_',num2str(itt),'(t_eom_interval)= nanmean(nanmean(data_temp(:,i:12:end)));']);  % extracted data OUT [sp_gy, T, sig]
                end
        end
%         end
    end
end

sp_9p_st = [1:6, 14:16];
%% each std and mean cal.
for varlist = { 'do','chl','ss','si','po4','no3','nh4','din','temp','salt','secchi'}
    for itt= 1:3 %regime
%         for sig = 1:3 %sigma
            for j=1:length(sp_9p_st) %% num of spatial point
                for i = 1:12
                     if i==1
                         clearvars t_eom_interval
                         t_eom_interval = 1:eom_d_each(1,i);
                     else    
                         clearvars t_eom_interval
                         t_eom_interval = eom_d_each(1,i-1)+1:eom_d_each(1,i);
                     end
                        %% time & spatial std & mean
                        clearvars len_temp data_temp
                        eval(['data_temp = squeeze(mtx_regime_',varname,'_',num2str(itt),'_s(sp_9p_st(j),:));']);  % extracted data OUT [sp_gy, T, sig]
                        len_temp=size(data_temp(:,i:12:end),1)*size(data_temp(:,i:12:end),2);
                        % time & spatial std 
                        eval(['eachst_std_gy_',varname,'_ft_',num2str(itt),'(t_eom_interval)= nanstd(data_temp(i:12:end));']);  % extracted data OUT [sp_gy, T, sig]
                        % mean
                        eval(['eachst_clim_gy_',varname,'_ft_',num2str(itt),'(t_eom_interval)= nanmean(data_temp(i:12:end));']);  % extracted data OUT [sp_gy, T, sig]
                        if strcmp(varname,'secchi') == 0 % there are no secchi_bot
                           clearvars len_temp data_temp
                           eval(['data_temp = squeeze(mtx_regime_',varname,'_',num2str(itt),'_b_s(sp_9p_st(j),:));']);  % extracted data OUT [sp_gy, T, sig]
                           % time & spatial std 
                           eval(['eachst_std_gy_',varname,'_b_ft_',num2str(itt),'(t_eom_interval)= nanstd(data_temp(i:12:end));']);  % extracted data OUT [sp_gy, T, sig]
                           % mean
                           eval(['eachst_clim_gy_',varname,'_b_ft_',num2str(itt),'(t_eom_interval)= nanmean(data_temp(i:12:end));']);  % extracted data OUT [sp_gy, T, sig]
                        end
                end
            end
%         end
    end
end                
save(['koem_climate_to06to15_v5_nosig_gy_9points.mat'],'obs*','obm*','spm_obs*','tm_obs*'); % full 9points data (no extract_sigma_for_each_st)
end


% right_here
%% make linear coeff. and save file

%% salt
color_pic = lines(size(regime_salt_yr,1));
marker_sty = {'o','+','x','^','>','h','p','s','d','.','*','v','<','o','+','x','^','>','h','p','s','d','.','*','v','<',...
    'o','+','x','^','>','h','p','s','d','.','*','v','<','o','+','x','^','>','h','p','s','d','.','*','v','<','o','+','x','^','>','h','p','s','d','.','*','v','<'};
xp = 1:22;
j=0
figure; hold on;
for i = 1:size(regime_salt_yr,1)
  clearvars reg_data_salt xp_w_salt pf_w_salt
reg_data_salt = regime_salt_yr(i,:);
if isnan(reg_data_salt(1)) == 0 
    j = j+1;
    xp_w_salt = find(isnan(reg_data_salt)==0);
    pf_w_salt = polyfit(xp_w_salt,reg_data_salt(xp_w_salt),1);
    yp_w_salt(i,:) = polyval(pf_w_salt,xp);
    scatter(1:22,regime_salt_yr(i,:),marker_sty{i},'MarkerEdgeColor',color_pic(i,:));
    plot(1:22, yp_w_salt(i,:),'color',color_pic(i,:));
    coeff_salt(i,:) = pf_w_salt;
    temp_case(j) = i;
end
hold on
end
xlabel('time(year)','fontsize',13)
ylabel('salt (mg/m^3)','fontsize',13)
set(gca,'xtick',[1:2:22]);
set(gca,'xlim',[1 22]);
set(gca,'xticklabel',1997:2:2018);
title('KOEM-Ç¥Ãş¿°ºĞ ¿¬Æò±Õ','fontsize',13)
grid on
set(gca,'fontsize',13)
% ylim([32 35])


%% no3
j=0
figure; hold on;
for i = 1:size(regime_no3_yr,1)
  clearvars reg_data_no3 xp_w_no3 pf_w_no3
reg_data_no3 = regime_no3_yr(i,:);
if isnan(reg_data_no3(1)) == 0 
    j = j+1;
    xp_w_no3 = find(isnan(reg_data_no3)==0);
    pf_w_no3 = polyfit(xp_w_no3,reg_data_no3(xp_w_no3),1);
    yp_w_no3(i,:) = polyval(pf_w_no3,xp);
    scatter(1:22,regime_no3_yr(i,:),marker_sty{i},'MarkerEdgeColor',color_pic(i,:));
    plot(1:22, yp_w_no3(i,:),'color',color_pic(i,:));
    coeff_no3(i,:) = pf_w_no3;
    temp_case(j) = i;
end
hold on
end
xlabel('time(year)','fontsize',13)
ylabel('NO3 (ug/L)','fontsize',13)
set(gca,'xtick',[1:2:22]);
set(gca,'xlim',[1 22]);
set(gca,'xticklabel',1997:2:2018);
title('KOEM°üÃø-Ç¥ÃşÁú»ê¿° ¿¬Æò±Õ','fontsize',13)
grid on
set(gca,'fontsize',13)
% ylim([32 35])
% legend('205-01','205-02','205-03','205-04','205-05')

%temp
j=0
figure; hold on;
for i = 1:size(regime_temp_yr,1)
  clearvars reg_data_temp xp_w_temp pf_w_temp
reg_data_temp = regime_temp_yr(i,:);
if isnan(reg_data_temp(1)) == 0 
    j = j+1;
    xp_w_temp = find(isnan(reg_data_temp)==0);
    pf_w_temp = polyfit(xp_w_temp,reg_data_temp(xp_w_temp),1);
    yp_w_temp(i,:) = polyval(pf_w_temp,xp);
    scatter(1:22,regime_temp_yr(i,:),marker_sty{i},'MarkerEdgeColor',color_pic(i,:));
    plot(1:22, yp_w_temp(i,:),'color',color_pic(i,:));
    coeff_temp(i,:) = pf_w_temp;
    temp_case(j) = i;
end
hold on
end
xlabel('time(year)','fontsize',13)
ylabel('temperature (^oC)','fontsize',13)
set(gca,'xtick',[1:2:22]);
set(gca,'xlim',[1 22]);
set(gca,'xticklabel',1997:2:2018);
title('KOEM°üÃø-Ç¥Ãş¼ö¿Â ¿¬Æò±Õ','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([13 20])

save('kodc_just_mean_yearly_linear_trend.mat','-v7.3');

filename=['KOEM_chl_bot_' num2str(i) 'mth']; 

print('-dpng',filename); 


close all
