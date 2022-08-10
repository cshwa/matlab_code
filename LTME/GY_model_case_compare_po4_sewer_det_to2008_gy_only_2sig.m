close all; clear;clc

cd D:\장기생태\observation\관측\장기생태프로젝트\남해_DB
ltme =load('LTME_observation_data_v2.mat'); %LTME obs.

% cd D:\장기생태\Dynamic\KOEM\gy_2001\gy_2001_koem_daily\spmean
% load('koem_result_processed_std_spmean.mat'); % KOEM obs.
% % mpdata= load('mpdata_result.mat');

cd D:\장기생태\Dynamic\result\2001
c2001 = load('2001_mp_p_sewer_det_result.mat','gy_*','jj_*','sgy_*','egy_*','eom_d_each');
c2001ob = load('2001_koem_compare_data_gy_2sig.mat','obm_*','obs_std_*');
cd D:\장기생태\Dynamic\result\2002
c2002 = load('2002_mp_p_sewer_det_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each');
c2002ob = load('2002_koem_compare_data_gy_2sig.mat','obm_*','obs_std_*');
cd D:\장기생태\Dynamic\result\2003
c2003 = load('2003_mp_p_sewer_det_f_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each');
c2003ob = load('2003_koem_compare_data_gy_2sig.mat','obm_*','obs_std_*');
cd D:\장기생태\Dynamic\result\2004
c2004 = load('2004_mp_p_sewer_det_f_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each');
c2004ob = load('2004_koem_compare_data_gy_2sig.mat','obm_*','obs_std_*');
cd D:\장기생태\Dynamic\result\2005
c2005 = load('2005_mp_p_sewer_det_f_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each');
c2005ob = load('2005_koem_compare_data_gy_2sig.mat','obm_*','obs_std_*');
cd D:\장기생태\Dynamic\result\2006
c2006 = load('2006_mp_p_sewer_det_f_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each');
c2006ob = load('2006_koem_compare_data_gy_2sig.mat','obm_*','obs_std_*');
cd D:\장기생태\Dynamic\result\2007
c2007 = load('2007_mp_p_sewer_det_f_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each');
c2007ob = load('2007_koem_compare_data_gy_2sig.mat','obm_*','obs_std_*');
cd D:\장기생태\Dynamic\result\2008
c2008 = load('2008_mp_p_sewer_det_f_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each');
c2008ob = load('2008_koem_compare_data_gy_2sig.mat','obm_*','obs_std_*');

cd D:\장기생태\Dynamic\koem
clim1=load('koem_climate_3regime_v2_to_2003_12.mat');


% cd D:\장기생태\Dynamic\KOEM\gy_2001\phymin_half_saturation
% phymin= load('phymin_result.mat');
% 
% cd D:\장기생태\Dynamic\KOEM\gy_2001\mp_p_ws_half_saturation
% mp_p_ws= load('mp_p_ws_result.mat');

% % mpdata.gy_temp
cd D:\장기생태\Dynamic\KOEM\compare_sewer\det_2001~8_gy_2sig


%make time tick
for i = 1:12
t_tIc{i} = {num2str(i,'%02d')};
if i == 1
    t_tIc{i} = {['2001.' num2str(i,'%02d')]};
end
end

for i = 1:12
t_tIc2{i} = {num2str(i,'%02d')};
if i == 1
    t_tIc2{i} = {['2002.' num2str(i,'%02d')]};
end
end

% for i = 1:12
% t_tIc2{i} = {num2str(i,'%02d')};
% if i == 1
%     t_tIc3{i} = {['2003.' num2str(i,'%02d')]};
% end
% end


c2001_eom = [1 (c2001.eom_d_each(1:end-1) + 1)];
c2002_eom = [(c2001.eom_d_each(end)+1) (c2001.eom_d_each(end)+c2002.eom_d_each(1:end-1) + 1)];
c2003_eom = [(c2001.eom_d_each(end)+c2002.eom_d_each(end)+1) (c2001.eom_d_each(end)+c2002.eom_d_each(end)+c2003.eom_d_each(1:end-1) + 1)];
c2004_eom = [(c2001.eom_d_each(end)+c2002.eom_d_each(end)+c2003.eom_d_each(end)+1) (c2001.eom_d_each(end)+c2002.eom_d_each(end)+c2003.eom_d_each(end) + c2004.eom_d_each(1:end-1)+ 1)];
c2005_eom = [(c2001.eom_d_each(end)+c2002.eom_d_each(end)+c2003.eom_d_each(end)+c2004.eom_d_each(end)+1) (c2001.eom_d_each(end)+c2002.eom_d_each(end)+c2003.eom_d_each(end) + c2004.eom_d_each(end) + c2005.eom_d_each(1:end-1)+ 1)];
c2006_eom = [(c2001.eom_d_each(end)+c2002.eom_d_each(end)+c2003.eom_d_each(end)+c2004.eom_d_each(end)+c2005.eom_d_each(end)+1) (c2001.eom_d_each(end)+c2002.eom_d_each(end)+c2003.eom_d_each(end) + c2004.eom_d_each(end) + c2005.eom_d_each(end)+ c2006.eom_d_each(1:end-1)+ 1)];
c2007_eom = [(c2001.eom_d_each(end)+c2002.eom_d_each(end)+c2003.eom_d_each(end)+c2004.eom_d_each(end)+c2005.eom_d_each(end)+c2006.eom_d_each(end)+1) (c2001.eom_d_each(end)+c2002.eom_d_each(end)+c2003.eom_d_each(end) + c2004.eom_d_each(end) + c2005.eom_d_each(end)+ c2006.eom_d_each(end)+c2007.eom_d_each(1:end-1)+ 1)];
c2008_eom = [(c2001.eom_d_each(end)+c2002.eom_d_each(end)+c2003.eom_d_each(end)+c2004.eom_d_each(end)+c2005.eom_d_each(end)+c2006.eom_d_each(end)+c2007.eom_d_each(end)+1) (c2001.eom_d_each(end)+c2002.eom_d_each(end)+c2003.eom_d_each(end) + c2004.eom_d_each(end) + c2005.eom_d_each(end)+ c2006.eom_d_each(end)+c2007.eom_d_each(end)+c2008.eom_d_each(1:end-1)+ 1)];

com_eom = [c2001_eom c2002_eom c2003_eom c2004_eom c2005_eom c2006_eom c2007_eom c2008_eom];

% % for i = 1:length(t_tIc)+length(t_tIc2)
% %     if i < 13
% %         com_t_tic{i}=t_tIc{i}
% %     else
% %         com_t_tic{i}=t_tIc2{i-12}
% %     end
% % end

clearvars com_t_tic
% com_t_tic = {'2001.01','02','03','04','05','06','07','08','09','10','11','12', ...
%             '2002.01','02','03','04','05','06','07','08','09','10','11','12', ...
%             '2003.01','02','03','04','05','06','07','08','09','10','11','12',...}
%             '2004.01','02','03','04','05','06','07','08','09','10','11','12'}

        
com_t_tic = {'2001', ...
            '2002', ...
            '2003',...
            '2004',...
            '2005',...
            '2006',...
            '2007',...
            '2008'}
        
        
       
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%% 2001
        plot(c2001.gy_no3,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2001ob.obm_gy_no3)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2001ob.obm_gy_no3./14 + c2001ob.obs_std_gy_no3./14;
        lower_bound_plt = c2001ob.obm_gy_no3./14 - c2001ob.obs_std_gy_no3./14;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2001ob.obm_gy_no3./14,'r-','linew',2); 
        plot(1:length(c2001ob.obm_gy_no3),upper_bound_plt,'r-','linew',2);
        plot(1:length(c2001ob.obm_gy_no3),lower_bound_plt,'r-','linew',2);
        
%% 2002
        yr1_t=length(c2001.gy_no3)
        plot(yr1_t+1:yr1_t+length(c2002.gy_no3),c2002.gy_no3,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2002ob.obm_gy_no3)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2002ob.obm_gy_no3./14 + c2002ob.obs_std_gy_no3./14;
        lower_bound_plt = c2002ob.obm_gy_no3./14 - c2002ob.obs_std_gy_no3./14;
        for i = 1:4
        if i == 1
        patch([yr1_t+nonan_data_plt(1:discon_p(1)) fliplr(yr1_t+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([yr1_t+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(yr1_t+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(yr1_t+1:yr1_t+length(c2002.gy_no3),c2002ob.obm_gy_no3./14,'r-','linew',2);
        plot(yr1_t+1:yr1_t+length(c2002.gy_no3),upper_bound_plt,'r-','linew',2);
        plot(yr1_t+1:yr1_t+length(c2002.gy_no3),lower_bound_plt,'r-','linew',2);
        
 %% 2003
        yr1_t=length(c2001.gy_no3)
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003.gy_no3),c2003.gy_no3,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2003ob.obm_gy_no3)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2003ob.obm_gy_no3./14 + c2003ob.obs_std_gy_no3./14;
        lower_bound_plt = c2003ob.obm_gy_no3./14 - c2003ob.obs_std_gy_no3./14;
        for i = 1:4
        if i == 1
        patch([c2003_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2003_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2003_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2003_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_no3),c2003ob.obm_gy_no3./14,'r-','linew',2);
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_no3),upper_bound_plt,'r-','linew',2);
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_no3),lower_bound_plt,'r-','linew',2);
 %% 2004
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004.gy_no3),c2004.gy_no3,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2004ob.obm_gy_no3)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2004ob.obm_gy_no3./14 + c2004ob.obs_std_gy_no3./14;
        lower_bound_plt = c2004ob.obm_gy_no3./14 - c2004ob.obs_std_gy_no3./14;
        for i = 1:4
        if i == 1
        patch([c2004_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2004_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2004_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2004_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_no3),c2004ob.obm_gy_no3./14,'r-','linew',2);
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_no3),upper_bound_plt,'r-','linew',2);
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_no3),lower_bound_plt,'r-','linew',2);       
 %% 2005
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005.gy_no3),c2005.gy_no3,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2005ob.obm_gy_no3)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2005ob.obm_gy_no3./14 + c2005ob.obs_std_gy_no3./14;
        lower_bound_plt = c2005ob.obm_gy_no3./14 - c2005ob.obs_std_gy_no3./14;
        for i = 1:4
        if i == 1
        patch([c2005_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2005_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2005_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2005_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_no3),c2005ob.obm_gy_no3./14,'r-','linew',2);
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_no3),upper_bound_plt,'r-','linew',2);
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_no3),lower_bound_plt,'r-','linew',2); 
%% 2006
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006.gy_no3),c2006.gy_no3,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2006ob.obm_gy_no3)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2006ob.obm_gy_no3./14 + c2006ob.obs_std_gy_no3./14;
        lower_bound_plt = c2006ob.obm_gy_no3./14 - c2006ob.obs_std_gy_no3./14;
        for i = 1:4
        if i == 1
        patch([c2006_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2006_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2006_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2006_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_no3),c2006ob.obm_gy_no3./14,'r-','linew',2);
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_no3),upper_bound_plt,'r-','linew',2);
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_no3),lower_bound_plt,'r-','linew',2); 
%% 2007
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007.gy_no3),c2007.gy_no3,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2007ob.obm_gy_no3)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2007ob.obm_gy_no3./14 + c2007ob.obs_std_gy_no3./14;
        lower_bound_plt = c2007ob.obm_gy_no3./14 - c2007ob.obs_std_gy_no3./14;
        for i = 1:4
        if i == 1
        patch([c2007_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2007_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2007_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2007_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_no3),c2007ob.obm_gy_no3./14,'r-','linew',2);
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_no3),upper_bound_plt,'r-','linew',2);
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_no3),lower_bound_plt,'r-','linew',2); 
%% 2008
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008.gy_no3),c2008.gy_no3,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2008ob.obm_gy_no3)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2008ob.obm_gy_no3./14 + c2008ob.obs_std_gy_no3./14;
        lower_bound_plt = c2008ob.obm_gy_no3./14 - c2008ob.obs_std_gy_no3./14;
        for i = 1:4
        if i == 1
        patch([c2008_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2008_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2008_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2008_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_no3),c2008ob.obm_gy_no3./14,'r-','linew',2);
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_no3),upper_bound_plt,'r-','linew',2);
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_no3),lower_bound_plt,'r-','linew',2);         
        
        alpha(0.3) %transparency      
        
        title(['GY KOEM OBS vs. daily MODEL NO3']);
        xlabel('time','fontsize',13)
        ylabel('NO3 (umol/m^3)','fontsize',13)
        grid on
        ylim([0 75])
        xlim([1 c2008_eom(1)-1+length(c2008.gy_no3)])
        set(gca,'fontsize',13)
        xticks(com_eom(1:12:end))
        xticklabels(com_t_tic)
        %xtickangle(90)
        
        print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_no3'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%% 2001
        plot(c2001.gy_nh4,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2001ob.obm_gy_nh4)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2001ob.obm_gy_nh4./14 + c2001ob.obs_std_gy_nh4./14;
        lower_bound_plt = c2001ob.obm_gy_nh4./14 - c2001ob.obs_std_gy_nh4./14;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2001ob.obm_gy_nh4./14,'r-','linew',2); 
        plot(1:length(c2001ob.obm_gy_nh4),upper_bound_plt,'r-','linew',2);
        plot(1:length(c2001ob.obm_gy_nh4),lower_bound_plt,'r-','linew',2);
        
%% 2002
        yr1_t=length(c2001.gy_nh4)
        plot(yr1_t+1:yr1_t+length(c2002.gy_nh4),c2002.gy_nh4,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2002ob.obm_gy_nh4)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2002ob.obm_gy_nh4./14 + c2002ob.obs_std_gy_nh4./14;
        lower_bound_plt = c2002ob.obm_gy_nh4./14 - c2002ob.obs_std_gy_nh4./14;
        for i = 1:4
        if i == 1
        patch([yr1_t+nonan_data_plt(1:discon_p(1)) fliplr(yr1_t+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([yr1_t+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(yr1_t+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(yr1_t+1:yr1_t+length(c2002.gy_nh4),c2002ob.obm_gy_nh4./14,'r-','linew',2);
        plot(yr1_t+1:yr1_t+length(c2002.gy_nh4),upper_bound_plt,'r-','linew',2);
        plot(yr1_t+1:yr1_t+length(c2002.gy_nh4),lower_bound_plt,'r-','linew',2);
        
 %% 2003
        yr1_t=length(c2001.gy_nh4)
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003.gy_nh4),c2003.gy_nh4,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2003ob.obm_gy_nh4)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2003ob.obm_gy_nh4./14 + c2003ob.obs_std_gy_nh4./14;
        lower_bound_plt = c2003ob.obm_gy_nh4./14 - c2003ob.obs_std_gy_nh4./14;
        for i = 1:4
        if i == 1
        patch([c2003_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2003_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2003_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2003_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_nh4),c2003ob.obm_gy_nh4./14,'r-','linew',2);
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_nh4),upper_bound_plt,'r-','linew',2);
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_nh4),lower_bound_plt,'r-','linew',2);
        alpha(0.3) %transparency
%% 2004
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004.gy_nh4),c2004.gy_nh4,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2004ob.obm_gy_nh4)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2004ob.obm_gy_nh4./14 + c2004ob.obs_std_gy_nh4./14;
        lower_bound_plt = c2004ob.obm_gy_nh4./14 - c2004ob.obs_std_gy_nh4./14;
        for i = 1:4
        if i == 1
        patch([c2004_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2004_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2004_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2004_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_nh4),c2004ob.obm_gy_nh4./14,'r-','linew',2);
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_nh4),upper_bound_plt,'r-','linew',2);
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_nh4),lower_bound_plt,'r-','linew',2); 
 %% 2005
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005.gy_nh4),c2005.gy_nh4,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2005ob.obm_gy_nh4)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2005ob.obm_gy_nh4./14 + c2005ob.obs_std_gy_nh4./14;
        lower_bound_plt = c2005ob.obm_gy_nh4./14 - c2005ob.obs_std_gy_nh4./14;
        for i = 1:4
        if i == 1
        patch([c2005_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2005_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2005_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2005_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_nh4),c2005ob.obm_gy_nh4./14,'r-','linew',2);
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_nh4),upper_bound_plt,'r-','linew',2);
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_nh4),lower_bound_plt,'r-','linew',2); 
 %% 2006
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006.gy_nh4),c2006.gy_nh4,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2006ob.obm_gy_nh4)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2006ob.obm_gy_nh4./14 + c2006ob.obs_std_gy_nh4./14;
        lower_bound_plt = c2006ob.obm_gy_nh4./14 - c2006ob.obs_std_gy_nh4./14;
        for i = 1:4
        if i == 1
        patch([c2006_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2006_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2006_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2006_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_nh4),c2006ob.obm_gy_nh4./14,'r-','linew',2);
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_nh4),upper_bound_plt,'r-','linew',2);
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_nh4),lower_bound_plt,'r-','linew',2); 
 %% 2007
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007.gy_nh4),c2007.gy_nh4,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2007ob.obm_gy_nh4)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2007ob.obm_gy_nh4./14 + c2007ob.obs_std_gy_nh4./14;
        lower_bound_plt = c2007ob.obm_gy_nh4./14 - c2007ob.obs_std_gy_nh4./14;
        for i = 1:4
        if i == 1
        patch([c2007_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2007_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2007_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2007_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_nh4),c2007ob.obm_gy_nh4./14,'r-','linew',2);
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_nh4),upper_bound_plt,'r-','linew',2);
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_nh4),lower_bound_plt,'r-','linew',2); 
        
 %% 2008
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008.gy_nh4),c2008.gy_nh4,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2008ob.obm_gy_nh4)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2008ob.obm_gy_nh4./14 + c2008ob.obs_std_gy_nh4./14;
        lower_bound_plt = c2008ob.obm_gy_nh4./14 - c2008ob.obs_std_gy_nh4./14;
        for i = 1:4
        if i == 1
        patch([c2008_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2008_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2008_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2008_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_nh4),c2008ob.obm_gy_nh4./14,'r-','linew',2);
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_nh4),upper_bound_plt,'r-','linew',2);
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_nh4),lower_bound_plt,'r-','linew',2);         
        
         alpha(0.3) %transparency  
              
        title(['GY KOEM OBS vs. daily MODEL nh4']);
        xlabel('time','fontsize',13)
        ylabel('nh4 (umol/m^3)','fontsize',13)
        grid on
        ylim([0 14])
        yticks(1:1:14)
        xlim([1 c2008_eom(1)-1+length(c2008.gy_no3)])
        xticks(com_eom(1:12:end))
        xticklabels(com_t_tic)
        %xtickangle(90)
        set(gca,'fontsize',13)
        print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_nh4'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%% 2001
        plot(c2001.gy_chl,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2001ob.obm_gy_chl)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2001ob.obm_gy_chl + c2001ob.obs_std_gy_chl;
        lower_bound_plt = c2001ob.obm_gy_chl - c2001ob.obs_std_gy_chl;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2001ob.obm_gy_chl,'r-','linew',2); 
        plot(1:length(c2001ob.obm_gy_chl),upper_bound_plt,'r-','linew',2);
        plot(1:length(c2001ob.obm_gy_chl),lower_bound_plt,'r-','linew',2);
        
%% 2002
        yr1_t=length(c2001.gy_chl)
        plot(yr1_t+1:yr1_t+length(c2002.gy_chl),c2002.gy_chl,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2002ob.obm_gy_chl)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2002ob.obm_gy_chl + c2002ob.obs_std_gy_chl;
        lower_bound_plt = c2002ob.obm_gy_chl - c2002ob.obs_std_gy_chl;
        for i = 1:4
        if i == 1
        patch([yr1_t+nonan_data_plt(1:discon_p(1)) fliplr(yr1_t+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([yr1_t+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(yr1_t+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(yr1_t+1:yr1_t+length(c2002.gy_chl),c2002ob.obm_gy_chl,'r-','linew',2);
        plot(yr1_t+1:yr1_t+length(c2002.gy_chl),upper_bound_plt,'r-','linew',2);
        plot(yr1_t+1:yr1_t+length(c2002.gy_chl),lower_bound_plt,'r-','linew',2);
        
 %% 2003
        yr1_t=length(c2001.gy_chl)
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003.gy_chl),c2003.gy_chl,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2003ob.obm_gy_chl)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2003ob.obm_gy_chl + c2003ob.obs_std_gy_chl;
        lower_bound_plt = c2003ob.obm_gy_chl - c2003ob.obs_std_gy_chl;
        for i = 1:4
        if i == 1
        patch([c2003_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2003_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2003_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2003_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_chl),c2003ob.obm_gy_chl,'r-','linew',2);
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_chl),upper_bound_plt,'r-','linew',2);
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_chl),lower_bound_plt,'r-','linew',2);

%% 2004
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004.gy_chl),c2004.gy_chl,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2004ob.obm_gy_chl)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2004ob.obm_gy_chl + c2004ob.obs_std_gy_chl;
        lower_bound_plt = c2004ob.obm_gy_chl - c2004ob.obs_std_gy_chl;
        for i = 1:4
        if i == 1
        patch([c2004_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2004_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2004_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2004_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_chl),c2004ob.obm_gy_chl,'r-','linew',2);
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_chl),upper_bound_plt,'r-','linew',2);
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_chl),lower_bound_plt,'r-','linew',2);        
%% 2005
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005.gy_chl),c2005.gy_chl,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2005ob.obm_gy_chl)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2005ob.obm_gy_chl + c2005ob.obs_std_gy_chl;
        lower_bound_plt = c2005ob.obm_gy_chl - c2005ob.obs_std_gy_chl;
        for i = 1:4
        if i == 1
        patch([c2005_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2005_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2005_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2005_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_chl),c2005ob.obm_gy_chl,'r-','linew',2);
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_chl),upper_bound_plt,'r-','linew',2);
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_chl),lower_bound_plt,'r-','linew',2);       
%% 2006
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006.gy_chl),c2006.gy_chl,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2006ob.obm_gy_chl)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2006ob.obm_gy_chl + c2006ob.obs_std_gy_chl;
        lower_bound_plt = c2006ob.obm_gy_chl - c2006ob.obs_std_gy_chl;
        for i = 1:4
        if i == 1
        patch([c2006_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2006_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2006_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2006_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_chl),c2006ob.obm_gy_chl,'r-','linew',2);
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_chl),upper_bound_plt,'r-','linew',2);
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_chl),lower_bound_plt,'r-','linew',2);       
        
%% 2007
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007.gy_chl),c2007.gy_chl,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2007ob.obm_gy_chl)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2007ob.obm_gy_chl + c2007ob.obs_std_gy_chl;
        lower_bound_plt = c2007ob.obm_gy_chl - c2007ob.obs_std_gy_chl;
        for i = 1:4
        if i == 1
        patch([c2007_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2007_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2007_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2007_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_chl),c2007ob.obm_gy_chl,'r-','linew',2);
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_chl),upper_bound_plt,'r-','linew',2);
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_chl),lower_bound_plt,'r-','linew',2);       

%% 2008
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008.gy_chl),c2008.gy_chl,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2008ob.obm_gy_chl)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2008ob.obm_gy_chl + c2008ob.obs_std_gy_chl;
        lower_bound_plt = c2008ob.obm_gy_chl - c2008ob.obs_std_gy_chl;
        for i = 1:4
        if i == 1
        patch([c2008_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2008_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2008_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2008_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_chl),c2008ob.obm_gy_chl,'r-','linew',2);
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_chl),upper_bound_plt,'r-','linew',2);
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_chl),lower_bound_plt,'r-','linew',2);       
        
        alpha(0.3) %transparency      
        
        title(['GY KOEM OBS vs. daily MODEL chl']);
        xlabel('time','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        ylim([0 15])
        yticks(1:1:15)
        xlim([1 c2008_eom(1)-1+length(c2008.gy_no3)])
        xticks(com_eom(1:12:end))
        xticklabels(com_t_tic)
        %xtickangle(90)
        set(gca,'fontsize',13)
        print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_chl'),'-dpng')
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%% 2001
        plot(c2001.gy_temp,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2001ob.obm_gy_temp)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2001ob.obm_gy_temp + c2001ob.obs_std_gy_temp;
        lower_bound_plt = c2001ob.obm_gy_temp - c2001ob.obs_std_gy_temp;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2001ob.obm_gy_temp,'r-','linew',2); 
        plot(1:length(c2001ob.obm_gy_temp),upper_bound_plt,'r-','linew',2);
        plot(1:length(c2001ob.obm_gy_temp),lower_bound_plt,'r-','linew',2);
        
%% 2002
        yr1_t=length(c2001.gy_temp)
        plot(yr1_t+1:yr1_t+length(c2002.gy_temp),c2002.gy_temp,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2002ob.obm_gy_temp)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2002ob.obm_gy_temp + c2002ob.obs_std_gy_temp;
        lower_bound_plt = c2002ob.obm_gy_temp - c2002ob.obs_std_gy_temp;
        for i = 1:4
        if i == 1
        patch([yr1_t+nonan_data_plt(1:discon_p(1)) fliplr(yr1_t+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([yr1_t+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(yr1_t+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(yr1_t+1:yr1_t+length(c2002.gy_temp),c2002ob.obm_gy_temp,'r-','linew',2);
        plot(yr1_t+1:yr1_t+length(c2002.gy_temp),upper_bound_plt,'r-','linew',2);
        plot(yr1_t+1:yr1_t+length(c2002.gy_temp),lower_bound_plt,'r-','linew',2);
        
%% 2003
        yr1_t=length(c2001.gy_temp)
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003.gy_temp),c2003.gy_temp,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2003ob.obm_gy_temp)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2003ob.obm_gy_temp + c2003ob.obs_std_gy_temp;
        lower_bound_plt = c2003ob.obm_gy_temp - c2003ob.obs_std_gy_temp;
        for i = 1:4
        if i == 1
        patch([c2003_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2003_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2003_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2003_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_temp),c2003ob.obm_gy_temp,'r-','linew',2);
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_temp),upper_bound_plt,'r-','linew',2);
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_temp),lower_bound_plt,'r-','linew',2);
        
%% 2004
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004.gy_temp),c2004.gy_temp,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2004ob.obm_gy_temp)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2004ob.obm_gy_temp + c2004ob.obs_std_gy_temp;
        lower_bound_plt = c2004ob.obm_gy_temp - c2004ob.obs_std_gy_temp;
        for i = 1:4
        if i == 1
        patch([c2004_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2004_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2004_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2004_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_temp),c2004ob.obm_gy_temp,'r-','linew',2);
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_temp),upper_bound_plt,'r-','linew',2);
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_temp),lower_bound_plt,'r-','linew',2);
 %% 2005
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005.gy_temp),c2005.gy_temp,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2005ob.obm_gy_temp)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2005ob.obm_gy_temp + c2005ob.obs_std_gy_temp;
        lower_bound_plt = c2005ob.obm_gy_temp - c2005ob.obs_std_gy_temp;
        for i = 1:4
        if i == 1
        patch([c2005_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2005_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2005_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2005_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_temp),c2005ob.obm_gy_temp,'r-','linew',2);
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_temp),upper_bound_plt,'r-','linew',2);
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_temp),lower_bound_plt,'r-','linew',2);

 %% 2006
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006.gy_temp),c2006.gy_temp,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2006ob.obm_gy_temp)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2006ob.obm_gy_temp + c2006ob.obs_std_gy_temp;
        lower_bound_plt = c2006ob.obm_gy_temp - c2006ob.obs_std_gy_temp;
        for i = 1:4
        if i == 1
        patch([c2006_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2006_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2006_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2006_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_temp),c2006ob.obm_gy_temp,'r-','linew',2);
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_temp),upper_bound_plt,'r-','linew',2);
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_temp),lower_bound_plt,'r-','linew',2);
 
 %% 2007
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007.gy_temp),c2007.gy_temp,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2007ob.obm_gy_temp)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2007ob.obm_gy_temp + c2007ob.obs_std_gy_temp;
        lower_bound_plt = c2007ob.obm_gy_temp - c2007ob.obs_std_gy_temp;
        for i = 1:4
        if i == 1
        patch([c2007_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2007_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2007_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2007_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_temp),c2007ob.obm_gy_temp,'r-','linew',2);
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_temp),upper_bound_plt,'r-','linew',2);
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_temp),lower_bound_plt,'r-','linew',2);
        
%% 2008
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008.gy_temp),c2008.gy_temp,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2008ob.obm_gy_temp)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2008ob.obm_gy_temp + c2008ob.obs_std_gy_temp;
        lower_bound_plt = c2008ob.obm_gy_temp - c2008ob.obs_std_gy_temp;
        for i = 1:4
        if i == 1
        patch([c2008_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2008_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2008_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2008_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_temp),c2008ob.obm_gy_temp,'r-','linew',2);
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_temp),upper_bound_plt,'r-','linew',2);
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_temp),lower_bound_plt,'r-','linew',2);        
        
        alpha(0.3) %transparency      
        
        title(['GY KOEM OBS vs. daily MODEL temp']);
        xlabel('time','fontsize',13)
        ylabel('temp (o^C)','fontsize',13)
        grid on
        ylim([0 35])
        xlim([1 c2008_eom(1)-1+length(c2008.gy_no3)])
        xticks(com_eom(1:12:end))
        xticklabels(com_t_tic)
        %xtickangle(90)
        set(gca,'fontsize',13)
        print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_temp'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%% 2001
        plot(c2001.gy_salt,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2001ob.obm_gy_salt)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2001ob.obm_gy_salt + c2001ob.obs_std_gy_salt;
        lower_bound_plt = c2001ob.obm_gy_salt - c2001ob.obs_std_gy_salt;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2001ob.obm_gy_salt,'r-','linew',2); 
        plot(1:length(c2001ob.obm_gy_salt),upper_bound_plt,'r-','linew',2);
        plot(1:length(c2001ob.obm_gy_salt),lower_bound_plt,'r-','linew',2);
        
%% 2002
        yr1_t=length(c2001.gy_salt)
        plot(yr1_t+1:yr1_t+length(c2002.gy_salt),c2002.gy_salt,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2002ob.obm_gy_salt)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2002ob.obm_gy_salt + c2002ob.obs_std_gy_salt;
        lower_bound_plt = c2002ob.obm_gy_salt - c2002ob.obs_std_gy_salt;
        for i = 1:4
        if i == 1
        patch([yr1_t+nonan_data_plt(1:discon_p(1)) fliplr(yr1_t+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([yr1_t+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(yr1_t+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(yr1_t+1:yr1_t+length(c2002.gy_salt),c2002ob.obm_gy_salt,'r-','linew',2);
        plot(yr1_t+1:yr1_t+length(c2002.gy_salt),upper_bound_plt,'r-','linew',2);
        plot(yr1_t+1:yr1_t+length(c2002.gy_salt),lower_bound_plt,'r-','linew',2); 
%% 2003
        yr1_t=length(c2001.gy_salt)
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003.gy_salt),c2003.gy_salt,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2003ob.obm_gy_salt)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2003ob.obm_gy_salt + c2003ob.obs_std_gy_salt;
        lower_bound_plt = c2003ob.obm_gy_salt - c2003ob.obs_std_gy_salt;
        for i = 1:4
        if i == 1
        patch([c2003_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2003_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2003_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2003_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_salt),c2003ob.obm_gy_salt,'r-','linew',2);
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_salt),upper_bound_plt,'r-','linew',2);
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_salt),lower_bound_plt,'r-','linew',2);
 %% 2004
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004.gy_salt),c2004.gy_salt,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2004ob.obm_gy_salt)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2004ob.obm_gy_salt + c2004ob.obs_std_gy_salt;
        lower_bound_plt = c2004ob.obm_gy_salt - c2004ob.obs_std_gy_salt;
        for i = 1:4
        if i == 1
        patch([c2004_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2004_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2004_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2004_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_salt),c2004ob.obm_gy_salt,'r-','linew',2);
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_salt),upper_bound_plt,'r-','linew',2);
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_salt),lower_bound_plt,'r-','linew',2);
%% 2005
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005.gy_salt),c2005.gy_salt,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2005ob.obm_gy_salt)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2005ob.obm_gy_salt + c2005ob.obs_std_gy_salt;
        lower_bound_plt = c2005ob.obm_gy_salt - c2005ob.obs_std_gy_salt;
        for i = 1:4
        if i == 1
        patch([c2005_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2005_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2005_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2005_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_salt),c2005ob.obm_gy_salt,'r-','linew',2);
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_salt),upper_bound_plt,'r-','linew',2);
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_salt),lower_bound_plt,'r-','linew',2);

%% 2006
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006.gy_salt),c2006.gy_salt,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2006ob.obm_gy_salt)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2006ob.obm_gy_salt + c2006ob.obs_std_gy_salt;
        lower_bound_plt = c2006ob.obm_gy_salt - c2006ob.obs_std_gy_salt;
        for i = 1:4
        if i == 1
        patch([c2006_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2006_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2006_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2006_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_salt),c2006ob.obm_gy_salt,'r-','linew',2);
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_salt),upper_bound_plt,'r-','linew',2);
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_salt),lower_bound_plt,'r-','linew',2);

%% 2007
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007.gy_salt),c2007.gy_salt,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2007ob.obm_gy_salt)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2007ob.obm_gy_salt + c2007ob.obs_std_gy_salt;
        lower_bound_plt = c2007ob.obm_gy_salt - c2007ob.obs_std_gy_salt;
        for i = 1:4
        if i == 1
        patch([c2007_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2007_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2007_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2007_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_salt),c2007ob.obm_gy_salt,'r-','linew',2);
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_salt),upper_bound_plt,'r-','linew',2);
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_salt),lower_bound_plt,'r-','linew',2);

%% 2008
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008.gy_salt),c2008.gy_salt,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2008ob.obm_gy_salt)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2008ob.obm_gy_salt + c2008ob.obs_std_gy_salt;
        lower_bound_plt = c2008ob.obm_gy_salt - c2008ob.obs_std_gy_salt;
        for i = 1:4
        if i == 1
        patch([c2008_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2008_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2008_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2008_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_salt),c2008ob.obm_gy_salt,'r-','linew',2);
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_salt),upper_bound_plt,'r-','linew',2);
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_salt),lower_bound_plt,'r-','linew',2);
        
        alpha(0.3) %transparency      
        
        title(['GY KOEM OBS vs. daily MODEL salt']);
        xlabel('time','fontsize',13)
        ylabel('salt (psu)','fontsize',13)
        grid on
        ylim([0 35])
        xlim([1 c2008_eom(1)-1+length(c2008.gy_no3)])
        xticks(com_eom(1:12:end))
        xticklabels(com_t_tic)
        %xtickangle(90)
        set(gca,'fontsize',13)
        print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_salt'),'-dpng')            
 %% bot
 
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%% 2001
        plot(c2001.gy_no3_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2001ob.obm_gy_no3_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2001ob.obm_gy_no3_b./14 + c2001ob.obs_std_gy_no3_b./14;
        lower_bound_plt = c2001ob.obm_gy_no3_b./14 - c2001ob.obs_std_gy_no3_b./14;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2001ob.obm_gy_no3_b./14,'r-','linew',2); 
        plot(1:length(c2001ob.obm_gy_no3_b),upper_bound_plt,'r-','linew',2);
        plot(1:length(c2001ob.obm_gy_no3_b),lower_bound_plt,'r-','linew',2);
        
%% 2002
        yr1_t=length(c2001.gy_no3_b)
        plot(yr1_t+1:yr1_t+length(c2002.gy_no3_b),c2002.gy_no3_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2002ob.obm_gy_no3_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2002ob.obm_gy_no3_b./14 + c2002ob.obs_std_gy_no3_b./14;
        lower_bound_plt = c2002ob.obm_gy_no3_b./14 - c2002ob.obs_std_gy_no3_b./14;
        for i = 1:4
        if i == 1
        patch([yr1_t+nonan_data_plt(1:discon_p(1)) fliplr(yr1_t+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([yr1_t+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(yr1_t+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(yr1_t+1:yr1_t+length(c2002.gy_no3_b),c2002ob.obm_gy_no3_b./14,'r-','linew',2);
        plot(yr1_t+1:yr1_t+length(c2002.gy_no3_b),upper_bound_plt,'r-','linew',2);
        plot(yr1_t+1:yr1_t+length(c2002.gy_no3_b),lower_bound_plt,'r-','linew',2);
        
%% 2003
        yr1_t=length(c2001.gy_no3_b)
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003.gy_no3_b),c2003.gy_no3_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2003ob.obm_gy_no3_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2003ob.obm_gy_no3_b./14 + c2003ob.obs_std_gy_no3_b./14;
        lower_bound_plt = c2003ob.obm_gy_no3_b./14 - c2003ob.obs_std_gy_no3_b./14;
        for i = 1:4
        if i == 1
        patch([c2003_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2003_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2003_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2003_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_no3_b),c2003ob.obm_gy_no3_b./14,'r-','linew',2);
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_no3_b),upper_bound_plt,'r-','linew',2);
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_no3_b),lower_bound_plt,'r-','linew',2);
        
%% 2004
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004.gy_no3_b),c2004.gy_no3_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2004ob.obm_gy_no3_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2004ob.obm_gy_no3_b./14 + c2004ob.obs_std_gy_no3_b./14;
        lower_bound_plt = c2004ob.obm_gy_no3_b./14 - c2004ob.obs_std_gy_no3_b./14;
        for i = 1:4
        if i == 1
        patch([c2004_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2004_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2004_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2004_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_no3_b),c2004ob.obm_gy_no3_b./14,'r-','linew',2);
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_no3_b),upper_bound_plt,'r-','linew',2);
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_no3_b),lower_bound_plt,'r-','linew',2); 
%% 2005
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005.gy_no3_b),c2005.gy_no3_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2005ob.obm_gy_no3_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2005ob.obm_gy_no3_b./14 + c2005ob.obs_std_gy_no3_b./14;
        lower_bound_plt = c2005ob.obm_gy_no3_b./14 - c2005ob.obs_std_gy_no3_b./14;
        for i = 1:4
        if i == 1
        patch([c2005_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2005_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2005_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2005_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_no3_b),c2005ob.obm_gy_no3_b./14,'r-','linew',2);
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_no3_b),upper_bound_plt,'r-','linew',2);
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_no3_b),lower_bound_plt,'r-','linew',2); 

%% 2006
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006.gy_no3_b),c2006.gy_no3_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2006ob.obm_gy_no3_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2006ob.obm_gy_no3_b./14 + c2006ob.obs_std_gy_no3_b./14;
        lower_bound_plt = c2006ob.obm_gy_no3_b./14 - c2006ob.obs_std_gy_no3_b./14;
        for i = 1:4
        if i == 1
        patch([c2006_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2006_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2006_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2006_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_no3_b),c2006ob.obm_gy_no3_b./14,'r-','linew',2);
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_no3_b),upper_bound_plt,'r-','linew',2);
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_no3_b),lower_bound_plt,'r-','linew',2); 
        
%% 2007
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007.gy_no3_b),c2007.gy_no3_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2007ob.obm_gy_no3_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2007ob.obm_gy_no3_b./14 + c2007ob.obs_std_gy_no3_b./14;
        lower_bound_plt = c2007ob.obm_gy_no3_b./14 - c2007ob.obs_std_gy_no3_b./14;
        for i = 1:4
        if i == 1
        patch([c2007_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2007_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2007_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2007_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_no3_b),c2007ob.obm_gy_no3_b./14,'r-','linew',2);
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_no3_b),upper_bound_plt,'r-','linew',2);
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_no3_b),lower_bound_plt,'r-','linew',2); 
        
%% 2008
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008.gy_no3_b),c2008.gy_no3_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2008ob.obm_gy_no3_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2008ob.obm_gy_no3_b./14 + c2008ob.obs_std_gy_no3_b./14;
        lower_bound_plt = c2008ob.obm_gy_no3_b./14 - c2008ob.obs_std_gy_no3_b./14;
        for i = 1:4
        if i == 1
        patch([c2008_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2008_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2008_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2008_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_no3_b),c2008ob.obm_gy_no3_b./14,'r-','linew',2);
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_no3_b),upper_bound_plt,'r-','linew',2);
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_no3_b),lower_bound_plt,'r-','linew',2); 
      
        alpha(0.3) %transparency      
        
        title(['GY KOEM OBS vs. daily MODEL no3 bot.']);
        xlabel('time','fontsize',13)
        ylabel('no3 bot (umol/m^3)','fontsize',13)
        grid on
        ylim([0 75])
        xlim([1 c2008_eom(1)-1+length(c2008.gy_no3)])
        xticks(com_eom(1:12:end))
        xticklabels(com_t_tic)
        %xtickangle(90)
        set(gca,'fontsize',13)
        print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_no3_bot'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%% 2001
        plot(c2001.gy_nh4_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2001ob.obm_gy_nh4_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2001ob.obm_gy_nh4_b./14 + c2001ob.obs_std_gy_nh4_b./14;
        lower_bound_plt = c2001ob.obm_gy_nh4_b./14 - c2001ob.obs_std_gy_nh4_b./14;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2001ob.obm_gy_nh4_b./14,'r-','linew',2); 
        plot(1:length(c2001ob.obm_gy_nh4_b),upper_bound_plt,'r-','linew',2);
        plot(1:length(c2001ob.obm_gy_nh4_b),lower_bound_plt,'r-','linew',2);
        
%% 2002
        yr1_t=length(c2001.gy_nh4_b)
        plot(yr1_t+1:yr1_t+length(c2002.gy_nh4_b),c2002.gy_nh4_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2002ob.obm_gy_nh4_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2002ob.obm_gy_nh4_b./14 + c2002ob.obs_std_gy_nh4_b./14;
        lower_bound_plt = c2002ob.obm_gy_nh4_b./14 - c2002ob.obs_std_gy_nh4_b./14;
        for i = 1:4
        if i == 1
        patch([yr1_t+nonan_data_plt(1:discon_p(1)) fliplr(yr1_t+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([yr1_t+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(yr1_t+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(yr1_t+1:yr1_t+length(c2002.gy_nh4_b),c2002ob.obm_gy_nh4_b./14,'r-','linew',2);
        plot(yr1_t+1:yr1_t+length(c2002.gy_nh4_b),upper_bound_plt,'r-','linew',2);
        plot(yr1_t+1:yr1_t+length(c2002.gy_nh4_b),lower_bound_plt,'r-','linew',2);
        
%% 2003
        yr1_t=length(c2001.gy_nh4_b)
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003.gy_nh4_b),c2003.gy_nh4_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2003ob.obm_gy_nh4_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2003ob.obm_gy_nh4_b./14 + c2003ob.obs_std_gy_nh4_b./14;
        lower_bound_plt = c2003ob.obm_gy_nh4_b./14 - c2003ob.obs_std_gy_nh4_b./14;
        for i = 1:4
        if i == 1
        patch([c2003_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2003_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2003_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2003_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_nh4_b),c2003ob.obm_gy_nh4_b./14,'r-','linew',2);
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_nh4_b),upper_bound_plt,'r-','linew',2);
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_nh4_b),lower_bound_plt,'r-','linew',2);
        
%% 2004
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004.gy_nh4_b),c2004.gy_nh4_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2004ob.obm_gy_nh4_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2004ob.obm_gy_nh4_b./14 + c2004ob.obs_std_gy_nh4_b./14;
        lower_bound_plt = c2004ob.obm_gy_nh4_b./14 - c2004ob.obs_std_gy_nh4_b./14;
        for i = 1:4
        if i == 1
        patch([c2004_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2004_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2004_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2004_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_nh4_b),c2004ob.obm_gy_nh4_b./14,'r-','linew',2);
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_nh4_b),upper_bound_plt,'r-','linew',2);
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_nh4_b),lower_bound_plt,'r-','linew',2); 
%% 2005
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005.gy_nh4_b),c2005.gy_nh4_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2005ob.obm_gy_nh4_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2005ob.obm_gy_nh4_b./14 + c2005ob.obs_std_gy_nh4_b./14;
        lower_bound_plt = c2005ob.obm_gy_nh4_b./14 - c2005ob.obs_std_gy_nh4_b./14;
        for i = 1:4
        if i == 1
        patch([c2005_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2005_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2005_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2005_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_nh4_b),c2005ob.obm_gy_nh4_b./14,'r-','linew',2);
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_nh4_b),upper_bound_plt,'r-','linew',2);
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_nh4_b),lower_bound_plt,'r-','linew',2); 

%% 2006
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006.gy_nh4_b),c2006.gy_nh4_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2006ob.obm_gy_nh4_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2006ob.obm_gy_nh4_b./14 + c2006ob.obs_std_gy_nh4_b./14;
        lower_bound_plt = c2006ob.obm_gy_nh4_b./14 - c2006ob.obs_std_gy_nh4_b./14;
        for i = 1:4
        if i == 1
        patch([c2006_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2006_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2006_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2006_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_nh4_b),c2006ob.obm_gy_nh4_b./14,'r-','linew',2);
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_nh4_b),upper_bound_plt,'r-','linew',2);
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_nh4_b),lower_bound_plt,'r-','linew',2); 

%% 2007
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007.gy_nh4_b),c2007.gy_nh4_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2007ob.obm_gy_nh4_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2007ob.obm_gy_nh4_b./14 + c2007ob.obs_std_gy_nh4_b./14;
        lower_bound_plt = c2007ob.obm_gy_nh4_b./14 - c2007ob.obs_std_gy_nh4_b./14;
        for i = 1:4
        if i == 1
        patch([c2007_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2007_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2007_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2007_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_nh4_b),c2007ob.obm_gy_nh4_b./14,'r-','linew',2);
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_nh4_b),upper_bound_plt,'r-','linew',2);
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_nh4_b),lower_bound_plt,'r-','linew',2); 

%% 2008
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008.gy_nh4_b),c2008.gy_nh4_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2008ob.obm_gy_nh4_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2008ob.obm_gy_nh4_b./14 + c2008ob.obs_std_gy_nh4_b./14;
        lower_bound_plt = c2008ob.obm_gy_nh4_b./14 - c2008ob.obs_std_gy_nh4_b./14;
        for i = 1:4
        if i == 1
        patch([c2008_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2008_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2008_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2008_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_nh4_b),c2008ob.obm_gy_nh4_b./14,'r-','linew',2);
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_nh4_b),upper_bound_plt,'r-','linew',2);
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_nh4_b),lower_bound_plt,'r-','linew',2); 

        alpha(0.3) %transparency      
        
        title(['GY KOEM OBS vs. daily MODEL nh4 bot.']);
        xlabel('time','fontsize',13)
        ylabel('nh4 bot (umol/m^3)','fontsize',13)
        grid on
        ylim([0 14])
        yticks(1:1:14)
        xlim([1 c2008_eom(1)-1+length(c2008.gy_no3)])
        xticks(com_eom(1:12:end))
        xticklabels(com_t_tic)
        %xtickangle(90)
        set(gca,'fontsize',13)
        print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_nh4_bot'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%% 2001
        plot(c2001.gy_chl_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2001ob.obm_gy_chl_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2001ob.obm_gy_chl_b + c2001ob.obs_std_gy_chl_b;
        lower_bound_plt = c2001ob.obm_gy_chl_b - c2001ob.obs_std_gy_chl_b;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2001ob.obm_gy_chl_b,'r-','linew',2); 
        plot(1:length(c2001ob.obm_gy_chl_b),upper_bound_plt,'r-','linew',2);
        plot(1:length(c2001ob.obm_gy_chl_b),lower_bound_plt,'r-','linew',2);
        
%% 2002
        yr1_t=length(c2001.gy_chl_b)
        plot(yr1_t+1:yr1_t+length(c2002.gy_chl_b),c2002.gy_chl_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2002ob.obm_gy_chl_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2002ob.obm_gy_chl_b + c2002ob.obs_std_gy_chl_b;
        lower_bound_plt = c2002ob.obm_gy_chl_b - c2002ob.obs_std_gy_chl_b;
        for i = 1:4
        if i == 1
        patch([yr1_t+nonan_data_plt(1:discon_p(1)) fliplr(yr1_t+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([yr1_t+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(yr1_t+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(yr1_t+1:yr1_t+length(c2002.gy_chl_b),c2002ob.obm_gy_chl_b,'r-','linew',2);
        plot(yr1_t+1:yr1_t+length(c2002.gy_chl_b),upper_bound_plt,'r-','linew',2);
        plot(yr1_t+1:yr1_t+length(c2002.gy_chl_b),lower_bound_plt,'r-','linew',2);
        
         %% 2003
        yr1_t=length(c2001.gy_chl_b)
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003.gy_chl_b),c2003.gy_chl_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2003ob.obm_gy_chl_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2003ob.obm_gy_chl_b + c2003ob.obs_std_gy_chl_b;
        lower_bound_plt = c2003ob.obm_gy_chl_b - c2003ob.obs_std_gy_chl_b;
        for i = 1:4
        if i == 1
        patch([c2003_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2003_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2003_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2003_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_chl_b),c2003ob.obm_gy_chl_b,'r-','linew',2);
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_chl_b),upper_bound_plt,'r-','linew',2);
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_chl_b),lower_bound_plt,'r-','linew',2);
        
%% 2004
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004.gy_chl_b),c2004.gy_chl_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2004ob.obm_gy_chl_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2004ob.obm_gy_chl_b + c2004ob.obs_std_gy_chl_b;
        lower_bound_plt = c2004ob.obm_gy_chl_b - c2004ob.obs_std_gy_chl_b;
        for i = 1:4
        if i == 1
        patch([c2004_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2004_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2004_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2004_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_chl_b),c2004ob.obm_gy_chl_b,'r-','linew',2);
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_chl_b),upper_bound_plt,'r-','linew',2);
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_chl_b),lower_bound_plt,'r-','linew',2); 
%% 2005
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005.gy_chl_b),c2005.gy_chl_b,'b','linew',2);
%         clearvars *_bound_plt nonan_data_plt discon_p
%         nonan_data_plt = find(isnan(c2005ob.obm_gy_chl_b)==0);
%         nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%         discon_p = find(nonan_diff ~= 1); % find discontinuous point
%         discon_p(4) = length(nonan_data_plt); %end point
%         upper_bound_plt = c2005ob.obm_gy_chl_b + c2005ob.obs_std_gy_chl_b;
%         lower_bound_plt = c2005ob.obm_gy_chl_b - c2005ob.obs_std_gy_chl_b;
%         for i = 1:4
%         if i == 1
%         patch([c2005_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2005_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
%         else
%         patch([c2005_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2005_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
%         end
%         end 
%         plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_chl_b),c2005ob.obm_gy_chl_b,'r-','linew',2);
%         plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_chl_b),upper_bound_plt,'r-','linew',2);
%         plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_chl_b),lower_bound_plt,'r-','linew',2); 

%% 2006
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006.gy_chl_b),c2006.gy_chl_b,'b','linew',2);
%         clearvars *_bound_plt nonan_data_plt discon_p
%         nonan_data_plt = find(isnan(c2006ob.obm_gy_chl_b)==0);
%         nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%         discon_p = find(nonan_diff ~= 1); % find discontinuous point
%         discon_p(4) = length(nonan_data_plt); %end point
%         upper_bound_plt = c2006ob.obm_gy_chl_b + c2006ob.obs_std_gy_chl_b;
%         lower_bound_plt = c2006ob.obm_gy_chl_b - c2006ob.obs_std_gy_chl_b;
%         for i = 1:4
%         if i == 1
%         patch([c2006_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2006_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
%         else
%         patch([c2006_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2006_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
%         end
%         end 
%         plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_chl_b),c2006ob.obm_gy_chl_b,'r-','linew',2);
%         plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_chl_b),upper_bound_plt,'r-','linew',2);
%         plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_chl_b),lower_bound_plt,'r-','linew',2); 

%% 2007
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007.gy_chl_b),c2007.gy_chl_b,'b','linew',2);
%         clearvars *_bound_plt nonan_data_plt discon_p
%         nonan_data_plt = find(isnan(c2007ob.obm_gy_chl_b)==0);
%         nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%         discon_p = find(nonan_diff ~= 1); % find discontinuous point
%         discon_p(4) = length(nonan_data_plt); %end point
%         upper_bound_plt = c2007ob.obm_gy_chl_b + c2007ob.obs_std_gy_chl_b;
%         lower_bound_plt = c2007ob.obm_gy_chl_b - c2007ob.obs_std_gy_chl_b;
%         for i = 1:4
%         if i == 1
%         patch([c2007_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2007_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
%         else
%         patch([c2007_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2007_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
%         end
%         end 
%         plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_chl_b),c2007ob.obm_gy_chl_b,'r-','linew',2);
%         plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_chl_b),upper_bound_plt,'r-','linew',2);
%         plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_chl_b),lower_bound_plt,'r-','linew',2); 

%% 2008
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008.gy_chl_b),c2008.gy_chl_b,'b','linew',2);
%         clearvars *_bound_plt nonan_data_plt discon_p
%         nonan_data_plt = find(isnan(c2008ob.obm_gy_chl_b)==0);
%         nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%         discon_p = find(nonan_diff ~= 1); % find discontinuous point
%         discon_p(4) = length(nonan_data_plt); %end point
%         upper_bound_plt = c2008ob.obm_gy_chl_b + c2008ob.obs_std_gy_chl_b;
%         lower_bound_plt = c2008ob.obm_gy_chl_b - c2008ob.obs_std_gy_chl_b;
%         for i = 1:4
%         if i == 1
%         patch([c2008_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2008_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
%         else
%         patch([c2008_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2008_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
%         end
%         end 
%         plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_chl_b),c2008ob.obm_gy_chl_b,'r-','linew',2);
%         plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_chl_b),upper_bound_plt,'r-','linew',2);
%         plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_chl_b),lower_bound_plt,'r-','linew',2); 
       
        alpha(0.3) %transparency      
        
        title(['GY KOEM OBS vs. daily MODEL chl bot.']);
        xlabel('time','fontsize',13)
        ylabel('chl bot (ug/L)','fontsize',13)
        grid on
        ylim([0 15])
        xlim([1 c2008_eom(1)-1+length(c2008.gy_no3)])
        yticks(1:1:15)
        xticks(com_eom(1:12:end))
        xticklabels(com_t_tic)
        %xtickangle(90)
        set(gca,'fontsize',13)
        print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_chl_bot'),'-dpng')
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%% 2001
        plot(c2001.gy_temp_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2001ob.obm_gy_temp_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2001ob.obm_gy_temp_b + c2001ob.obs_std_gy_temp_b;
        lower_bound_plt = c2001ob.obm_gy_temp_b - c2001ob.obs_std_gy_temp_b;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2001ob.obm_gy_temp_b,'r-','linew',2); 
        plot(1:length(c2001ob.obm_gy_temp_b),upper_bound_plt,'r-','linew',2);
        plot(1:length(c2001ob.obm_gy_temp_b),lower_bound_plt,'r-','linew',2);
        
%% 2002
        yr1_t=length(c2001.gy_temp_b)
        plot(yr1_t+1:yr1_t+length(c2002.gy_temp_b),c2002.gy_temp_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2002ob.obm_gy_temp_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2002ob.obm_gy_temp_b + c2002ob.obs_std_gy_temp_b;
        lower_bound_plt = c2002ob.obm_gy_temp_b - c2002ob.obs_std_gy_temp_b;
        for i = 1:4
        if i == 1
        patch([yr1_t+nonan_data_plt(1:discon_p(1)) fliplr(yr1_t+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([yr1_t+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(yr1_t+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(yr1_t+1:yr1_t+length(c2002.gy_temp_b),c2002ob.obm_gy_temp_b,'r-','linew',2);
        plot(yr1_t+1:yr1_t+length(c2002.gy_temp_b),upper_bound_plt,'r-','linew',2);
        plot(yr1_t+1:yr1_t+length(c2002.gy_temp_b),lower_bound_plt,'r-','linew',2);
        
%% 2003
        yr1_t=length(c2001.gy_temp_b)
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003.gy_temp_b),c2003.gy_temp_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2003ob.obm_gy_temp_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2003ob.obm_gy_temp_b + c2003ob.obs_std_gy_temp_b;
        lower_bound_plt = c2003ob.obm_gy_temp_b - c2003ob.obs_std_gy_temp_b;
        for i = 1:4
        if i == 1
        patch([c2003_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2003_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2003_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2003_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_temp_b),c2003ob.obm_gy_temp_b,'r-','linew',2);
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_temp_b),upper_bound_plt,'r-','linew',2);
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_temp_b),lower_bound_plt,'r-','linew',2);
        
%% 2004
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004.gy_temp_b),c2004.gy_temp_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2004ob.obm_gy_temp_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2004ob.obm_gy_temp_b + c2004ob.obs_std_gy_temp_b;
        lower_bound_plt = c2004ob.obm_gy_temp_b - c2004ob.obs_std_gy_temp_b;
        for i = 1:4
        if i == 1
        patch([c2004_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2004_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2004_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2004_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_temp_b),c2004ob.obm_gy_temp_b,'r-','linew',2);
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_temp_b),upper_bound_plt,'r-','linew',2);
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_temp_b),lower_bound_plt,'r-','linew',2); 
%% 2005
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005.gy_temp_b),c2005.gy_temp_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2005ob.obm_gy_temp_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2005ob.obm_gy_temp_b + c2005ob.obs_std_gy_temp_b;
        lower_bound_plt = c2005ob.obm_gy_temp_b - c2005ob.obs_std_gy_temp_b;
        for i = 1:4
        if i == 1
        patch([c2005_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2005_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2005_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2005_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_temp_b),c2005ob.obm_gy_temp_b,'r-','linew',2);
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_temp_b),upper_bound_plt,'r-','linew',2);
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_temp_b),lower_bound_plt,'r-','linew',2); 
 
%% 2006
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006.gy_temp_b),c2006.gy_temp_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2006ob.obm_gy_temp_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2006ob.obm_gy_temp_b + c2006ob.obs_std_gy_temp_b;
        lower_bound_plt = c2006ob.obm_gy_temp_b - c2006ob.obs_std_gy_temp_b;
        for i = 1:4
        if i == 1
        patch([c2006_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2006_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2006_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2006_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_temp_b),c2006ob.obm_gy_temp_b,'r-','linew',2);
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_temp_b),upper_bound_plt,'r-','linew',2);
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_temp_b),lower_bound_plt,'r-','linew',2); 

%% 2007
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007.gy_temp_b),c2007.gy_temp_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2007ob.obm_gy_temp_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2007ob.obm_gy_temp_b + c2007ob.obs_std_gy_temp_b;
        lower_bound_plt = c2007ob.obm_gy_temp_b - c2007ob.obs_std_gy_temp_b;
        for i = 1:4
        if i == 1
        patch([c2007_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2007_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2007_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2007_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_temp_b),c2007ob.obm_gy_temp_b,'r-','linew',2);
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_temp_b),upper_bound_plt,'r-','linew',2);
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_temp_b),lower_bound_plt,'r-','linew',2); 

%% 2008
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008.gy_temp_b),c2008.gy_temp_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2008ob.obm_gy_temp_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2008ob.obm_gy_temp_b + c2008ob.obs_std_gy_temp_b;
        lower_bound_plt = c2008ob.obm_gy_temp_b - c2008ob.obs_std_gy_temp_b;
        for i = 1:4
        if i == 1
        patch([c2008_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2008_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2008_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2008_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_temp_b),c2008ob.obm_gy_temp_b,'r-','linew',2);
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_temp_b),upper_bound_plt,'r-','linew',2);
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_temp_b),lower_bound_plt,'r-','linew',2); 
       
        alpha(0.3) %transparency      
        
        title(['GY KOEM OBS vs. daily MODEL temp bot.']);
        xlabel('time','fontsize',13)
        ylabel('temp bot (o^C)','fontsize',13)
        grid on
        ylim([0 35])
        xlim([1 c2008_eom(1)-1+length(c2008.gy_no3)])
         xticks(com_eom(1:12:end))
        xticklabels(com_t_tic)
        %xtickangle(90)
        set(gca,'fontsize',13)
        print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_temp_bot'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%% 2001
        plot(c2001.gy_salt_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2001ob.obm_gy_salt_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2001ob.obm_gy_salt_b + c2001ob.obs_std_gy_salt_b;
        lower_bound_plt = c2001ob.obm_gy_salt_b - c2001ob.obs_std_gy_salt_b;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2001ob.obm_gy_salt_b,'r-','linew',2); 
        plot(1:length(c2001ob.obm_gy_salt_b),upper_bound_plt,'r-','linew',2);
        plot(1:length(c2001ob.obm_gy_salt_b),lower_bound_plt,'r-','linew',2);
        
%% 2002
        yr1_t=length(c2001.gy_salt_b)
        plot(yr1_t+1:yr1_t+length(c2002.gy_salt_b),c2002.gy_salt_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2002ob.obm_gy_salt_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2002ob.obm_gy_salt_b + c2002ob.obs_std_gy_salt_b;
        lower_bound_plt = c2002ob.obm_gy_salt_b - c2002ob.obs_std_gy_salt_b;
        for i = 1:4
        if i == 1
        patch([yr1_t+nonan_data_plt(1:discon_p(1)) fliplr(yr1_t+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([yr1_t+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(yr1_t+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(yr1_t+1:yr1_t+length(c2002.gy_salt_b),c2002ob.obm_gy_salt_b,'r-','linew',2);
        plot(yr1_t+1:yr1_t+length(c2002.gy_salt_b),upper_bound_plt,'r-','linew',2);
        plot(yr1_t+1:yr1_t+length(c2002.gy_salt_b),lower_bound_plt,'r-','linew',2);
        
%% 2003
        yr1_t=length(c2001.gy_salt_b)
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003.gy_salt_b),c2003.gy_salt_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2003ob.obm_gy_salt_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2003ob.obm_gy_salt_b + c2003ob.obs_std_gy_salt_b;
        lower_bound_plt = c2003ob.obm_gy_salt_b - c2003ob.obs_std_gy_salt_b;
        for i = 1:4
        if i == 1
        patch([c2003_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2003_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2003_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2003_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_salt_b),c2003ob.obm_gy_salt_b,'r-','linew',2);
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_salt_b),upper_bound_plt,'r-','linew',2);
        plot(c2003_eom(1):c2003_eom(1)-1+length(c2003ob.obm_gy_salt_b),lower_bound_plt,'r-','linew',2);
        
%% 2004
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004.gy_salt_b),c2004.gy_salt_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2004ob.obm_gy_salt_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2004ob.obm_gy_salt_b + c2004ob.obs_std_gy_salt_b;
        lower_bound_plt = c2004ob.obm_gy_salt_b - c2004ob.obs_std_gy_salt_b;
        for i = 1:4
        if i == 1
        patch([c2004_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2004_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2004_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2004_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_salt_b),c2004ob.obm_gy_salt_b,'r-','linew',2);
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_salt_b),upper_bound_plt,'r-','linew',2);
        plot(c2004_eom(1):c2004_eom(1)-1+length(c2004ob.obm_gy_salt_b),lower_bound_plt,'r-','linew',2); 
%% 2005
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005.gy_salt_b),c2005.gy_salt_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2005ob.obm_gy_salt_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2005ob.obm_gy_salt_b + c2005ob.obs_std_gy_salt_b;
        lower_bound_plt = c2005ob.obm_gy_salt_b - c2005ob.obs_std_gy_salt_b;
        for i = 1:4
        if i == 1
        patch([c2005_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2005_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2005_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2005_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_salt_b),c2005ob.obm_gy_salt_b,'r-','linew',2);
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_salt_b),upper_bound_plt,'r-','linew',2);
        plot(c2005_eom(1):c2005_eom(1)-1+length(c2005ob.obm_gy_salt_b),lower_bound_plt,'r-','linew',2); 
        
%% 2006
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006.gy_salt_b),c2006.gy_salt_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2006ob.obm_gy_salt_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2006ob.obm_gy_salt_b + c2006ob.obs_std_gy_salt_b;
        lower_bound_plt = c2006ob.obm_gy_salt_b - c2006ob.obs_std_gy_salt_b;
        for i = 1:4
        if i == 1
        patch([c2006_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2006_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2006_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2006_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_salt_b),c2006ob.obm_gy_salt_b,'r-','linew',2);
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_salt_b),upper_bound_plt,'r-','linew',2);
        plot(c2006_eom(1):c2006_eom(1)-1+length(c2006ob.obm_gy_salt_b),lower_bound_plt,'r-','linew',2); 
        
%% 2007
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007.gy_salt_b),c2007.gy_salt_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2007ob.obm_gy_salt_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2007ob.obm_gy_salt_b + c2007ob.obs_std_gy_salt_b;
        lower_bound_plt = c2007ob.obm_gy_salt_b - c2007ob.obs_std_gy_salt_b;
        for i = 1:4
        if i == 1
        patch([c2007_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2007_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2007_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2007_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_salt_b),c2007ob.obm_gy_salt_b,'r-','linew',2);
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_salt_b),upper_bound_plt,'r-','linew',2);
        plot(c2007_eom(1):c2007_eom(1)-1+length(c2007ob.obm_gy_salt_b),lower_bound_plt,'r-','linew',2); 
        
%% 2008
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008.gy_salt_b),c2008.gy_salt_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(c2008ob.obm_gy_salt_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = c2008ob.obm_gy_salt_b + c2008ob.obs_std_gy_salt_b;
        lower_bound_plt = c2008ob.obm_gy_salt_b - c2008ob.obs_std_gy_salt_b;
        for i = 1:4
        if i == 1
        patch([c2008_eom(1)-1+nonan_data_plt(1:discon_p(1)) fliplr(c2008_eom(1)-1+nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([c2008_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(c2008_eom(1)-1+nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_salt_b),c2008ob.obm_gy_salt_b,'r-','linew',2);
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_salt_b),upper_bound_plt,'r-','linew',2);
        plot(c2008_eom(1):c2008_eom(1)-1+length(c2008ob.obm_gy_salt_b),lower_bound_plt,'r-','linew',2); 
        
        alpha(0.3) %transparency      
        
        title(['GY KOEM OBS vs. daily MODEL salt bot.']);
        xlabel('time','fontsize',13)
        ylabel('salt bot (psu)','fontsize',13)
        grid on
        ylim([0 35])
        xlim([1 c2008_eom(1)-1+length(c2008.gy_no3)])
                xticks(com_eom(1:12:end))
        xticklabels(com_t_tic)
        %xtickangle(90)
        set(gca,'fontsize',13)
        print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_salt_bot'),'-dpng')   
        
return        
        
%%south gy plot

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_no3,'b','linew',2);
        plot(pcase.sgy_no3,'k','linew',2);
        plot(mp_p_ws.sgy_no3,'g','linew',2); 
        plot(obm_sgy_no3./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_no3),obm_sgy_no3./14 + obs_std_sgy_no3./14,'m-','linew',2);
        plot(1:length(obm_sgy_no3),obm_sgy_no3./14 - obs_std_sgy_no3./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NO3 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_south_sgy'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_nh4,'b','linew',2);
        plot(pcase.sgy_nh4,'k','linew',2);
        plot(mp_p_ws.sgy_nh4,'g','linew',2);
        plot(obm_sgy_nh4./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_nh4),obm_sgy_nh4./14 + obs_std_sgy_nh4./14,'m-','linew',2);
        plot(1:length(obm_sgy_nh4),obm_sgy_nh4./14 - obs_std_sgy_nh4./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NH4 (umol/m^3)','fontsize',16)
        xticks([0:50:366])
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_south_sgy'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_chl,'b','linew',2);
        plot(pcase.sgy_chl,'k','linew',2);
        plot(mp_p_ws.sgy_chl,'g','linew',2);
        plot(obm_sgy_chl,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_chl),obm_sgy_chl + obs_std_sgy_chl,'m-','linew',2);
        plot(1:length(obm_sgy_chl),obm_sgy_chl - obs_std_sgy_chl,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL chl']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('chl (ug/L)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_south_sgy'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_temp,'b','linew',2);
        plot(pcase.sgy_temp,'k','linew',2);
        plot(mp_p_ws.sgy_temp,'g','linew',2);
        plot(obm_sgy_temp,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_temp),obm_sgy_temp + obs_std_sgy_temp,'m-','linew',2);
        plot(1:length(obm_sgy_temp),obm_sgy_temp - obs_std_sgy_temp,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL temp']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('temp (o^C)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_south_sgy'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_salt,'b','linew',2);
        plot(pcase.sgy_salt,'k','linew',2);
        plot(mp_p_ws.sgy_salt,'g','linew',2);   
        plot(obm_sgy_salt,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_salt),obm_sgy_salt + obs_std_sgy_salt,'m-','linew',2);
        plot(1:length(obm_sgy_salt),obm_sgy_salt - obs_std_sgy_salt,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL salt']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('salt (psu)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_south_sgy'),'-dpng')            
 %% bot
 
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_no3_b,'b','linew',2);
        plot(pcase.sgy_no3_b,'k','linew',2);
        plot(mp_p_ws.sgy_no3_b,'g','linew',2);
        plot(obm_sgy_no3_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_no3_b),obm_sgy_no3_b./14 + obs_std_sgy_no3_b./14,'m-','linew',2);
        plot(1:length(obm_sgy_no3_b),obm_sgy_no3_b./14 - obs_std_sgy_no3_b./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3 bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NO3 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_bot_south_sgy'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_nh4_b,'b','linew',2);
        plot(pcase.sgy_nh4_b,'k','linew',2);
        plot(mp_p_ws.sgy_nh4_b,'g','linew',2); 
        plot(obm_sgy_nh4_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_nh4_b),obm_sgy_nh4_b./14 + obs_std_sgy_nh4_b./14,'m-','linew',2);
        plot(1:length(obm_sgy_nh4_b),obm_sgy_nh4_b./14 - obs_std_sgy_nh4_b./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4 bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NH4 (umol/m^3)','fontsize',16)
        xticks([0:50:366])
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_bot_south_sgy'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_chl_b,'b','linew',2);
        plot(pcase.sgy_chl_b,'k','linew',2);
        plot(mp_p_ws.sgy_chl_b,'g','linew',2);  
        plot(obm_sgy_chl_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_chl_b),obm_sgy_chl_b + obs_std_sgy_chl_b,'m-','linew',2);
        plot(1:length(obm_sgy_chl_b),obm_sgy_chl_b - obs_std_sgy_chl_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL chl bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('chl (ug/L)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_bot_south_sgy'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_temp_b,'b','linew',2);
        plot(pcase.sgy_temp_b,'k','linew',2);
        plot(mp_p_ws.sgy_temp_b,'g','linew',2);        
        plot(obm_sgy_temp_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_temp_b),obm_sgy_temp_b + obs_std_sgy_temp_b,'m-','linew',2);
        plot(1:length(obm_sgy_temp_b),obm_sgy_temp_b - obs_std_sgy_temp_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL temp bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('temp (o^C)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_bot_south_sgy'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_salt_b,'b','linew',2);
        plot(pcase.sgy_salt_b,'k','linew',2);
        plot(mp_p_ws.sgy_salt_b,'g','linew',2);    
        plot(obm_sgy_salt_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_salt_b),obm_sgy_salt_b + obs_std_sgy_salt_b,'m-','linew',2);
        plot(1:length(obm_sgy_salt_b),obm_sgy_salt_b - obs_std_sgy_salt_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL salt bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('salt (psu)','fontsize',16)
        xticks([0:50:366])
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_bot_south_sgy'),'-dpng')    
        
%%east gy plot

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_no3,'b','linew',2);
        plot(pcase.egy_no3,'k','linew',2);
        plot(mp_p_ws.egy_no3,'g','linew',2);  
        plot(obm_egy_no3./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_no3),obm_egy_no3./14 + obs_std_egy_no3./14,'m-','linew',2);
        plot(1:length(obm_egy_no3),obm_egy_no3./14 - obs_std_egy_no3./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NO3 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_egy'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines       
        plot(mpdata.egy_nh4,'b','linew',2);
        plot(pcase.egy_nh4,'k','linew',2);
        plot(mp_p_ws.egy_nh4,'g','linew',2); 
        plot(obm_egy_nh4./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_nh4),obm_egy_nh4./14 + obs_std_egy_nh4./14,'m-','linew',2);
        plot(1:length(obm_egy_nh4),obm_egy_nh4./14 - obs_std_egy_nh4./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NH4 (umol/m^3)','fontsize',16)
        xticks([0:50:366])
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_egy'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_chl,'b','linew',2);
        plot(pcase.egy_chl,'k','linew',2);
        plot(mp_p_ws.egy_chl,'g','linew',2); 
        plot(obm_egy_chl,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_chl),obm_egy_chl + obs_std_egy_chl,'m-','linew',2);
        plot(1:length(obm_egy_chl),obm_egy_chl - obs_std_egy_chl,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL chl']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('chl (ug/L)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_egy'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_temp,'b','linew',2);
        plot(pcase.egy_temp,'k','linew',2);
        plot(mp_p_ws.egy_temp,'g','linew',2); 
        plot(obm_egy_temp,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_temp),obm_egy_temp + obs_std_egy_temp,'m-','linew',2);
        plot(1:length(obm_egy_temp),obm_egy_temp - obs_std_egy_temp,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL temp']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('temp (o^C)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_egy'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_salt,'b','linew',2);
        plot(pcase.egy_salt,'k','linew',2);
        plot(mp_p_ws.egy_salt,'g','linew',2); 
        plot(obm_egy_salt,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_salt),obm_egy_salt + obs_std_egy_salt,'m-','linew',2);
        plot(1:length(obm_egy_salt),obm_egy_salt - obs_std_egy_salt,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL salt']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('salt (psu)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_egy'),'-dpng')            
 %% bot
 
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_no3_b,'b','linew',2);
        plot(pcase.egy_no3_b,'k','linew',2);
        plot(mp_p_ws.egy_no3_b,'g','linew',2); 
        plot(obm_egy_no3_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_no3_b),obm_egy_no3_b./14 + obs_std_egy_no3_b./14,'m-','linew',2);
        plot(1:length(obm_egy_no3_b),obm_egy_no3_b./14 - obs_std_egy_no3_b./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3 bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NO3 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_bot_egy'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_nh4_b,'b','linew',2);
        plot(pcase.egy_nh4_b,'k','linew',2);
        plot(mp_p_ws.egy_nh4_b,'g','linew',2); 
        plot(obm_egy_nh4_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_nh4_b),obm_egy_nh4_b./14 + obs_std_egy_nh4_b./14,'m-','linew',2);
        plot(1:length(obm_egy_nh4_b),obm_egy_nh4_b./14 - obs_std_egy_nh4_b./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4 bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NH4 (umol/m^3)','fontsize',16)
        xticks([0:50:366])
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_bot_egy'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_chl_b,'b','linew',2);
        plot(pcase.egy_chl_b,'k','linew',2);
        plot(mp_p_ws.egy_chl_b,'g','linew',2); 
        plot(obm_egy_chl_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_chl_b),obm_egy_chl_b + obs_std_egy_chl_b,'m-','linew',2);
        plot(1:length(obm_egy_chl_b),obm_egy_chl_b - obs_std_egy_chl_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL chl bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('chl (ug/L)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_bot_egy'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_temp_b,'b','linew',2);
        plot(pcase.egy_temp_b,'k','linew',2);
        plot(mp_p_ws.egy_temp_b,'g','linew',2); 
        plot(obm_egy_temp_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_temp_b),obm_egy_temp_b + obs_std_egy_temp_b,'m-','linew',2);
        plot(1:length(obm_egy_temp_b),obm_egy_temp_b - obs_std_egy_temp_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL temp bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('temp (o^C)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_bot_egy'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_salt_b,'b','linew',2);
        plot(pcase.egy_salt_b,'k','linew',2);
        plot(mp_p_ws.egy_salt_b,'g','linew',2); 
        plot(obm_egy_salt_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_salt_b),obm_egy_salt_b + obs_std_egy_salt_b,'m-','linew',2);
        plot(1:length(obm_egy_salt_b),obm_egy_salt_b - obs_std_egy_salt_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL salt bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('salt (psu)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_bot_egy'),'-dpng')    
        
%% jinju

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_no3,'b','linew',2);
        plot(pcase.jj_no3,'k','linew',2);
        plot(mp_p_ws.jj_no3,'g','linew',2); 
        plot(obm_jj_no3./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_no3),obm_jj_no3./14 + obs_std_jj_no3./14,'m-','linew',2);
        plot(1:length(obm_jj_no3),obm_jj_no3./14 - obs_std_jj_no3./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NO3 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_jj'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_nh4,'b','linew',2);
        plot(pcase.jj_nh4,'k','linew',2);
        plot(mp_p_ws.jj_nh4,'g','linew',2); 
        plot(obm_jj_nh4./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_nh4),obm_jj_nh4./14 + obs_std_jj_nh4./14,'m-','linew',2);
        plot(1:length(obm_jj_nh4),obm_jj_nh4./14 - obs_std_jj_nh4./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NH4 (umol/m^3)','fontsize',16)
        xticks([0:50:366])
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_jj'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_chl,'b','linew',2);
        plot(pcase.jj_chl,'k','linew',2);
        plot(mp_p_ws.jj_chl,'g','linew',2); 
        plot(obm_jj_chl,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_chl),obm_jj_chl + obs_std_jj_chl,'m-','linew',2);
        plot(1:length(obm_jj_chl),obm_jj_chl - obs_std_jj_chl,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL chl']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('chl (ug/L)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_jj'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_temp,'b','linew',2);
        plot(pcase.jj_temp,'k','linew',2);
        plot(mp_p_ws.jj_temp,'g','linew',2); 
        plot(obm_jj_temp,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_temp),obm_jj_temp + obs_std_jj_temp,'m-','linew',2);
        plot(1:length(obm_jj_temp),obm_jj_temp - obs_std_jj_temp,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL temp']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('temp (o^C)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_jj'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_salt,'b','linew',2);
        plot(pcase.jj_salt,'k','linew',2);
        plot(mp_p_ws.jj_salt,'g','linew',2); 
        plot(obm_jj_salt,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_salt),obm_jj_salt + obs_std_jj_salt,'m-','linew',2);
        plot(1:length(obm_jj_salt),obm_jj_salt - obs_std_jj_salt,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL salt']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('salt (psu)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_jj'),'-dpng')            
 %% bot
 
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_no3_b,'b','linew',2);
        plot(pcase.jj_no3_b,'k','linew',2);
        plot(mp_p_ws.jj_no3_b,'g','linew',2); 
        plot(obm_jj_no3_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_no3_b),obm_jj_no3_b./14 + obs_std_jj_no3_b./14,'m-','linew',2);
        plot(1:length(obm_jj_no3_b),obm_jj_no3_b./14 - obs_std_jj_no3_b./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3 bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NO3 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_bot_jj'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_nh4_b,'b','linew',2);
        plot(pcase.jj_nh4_b,'k','linew',2);
        plot(mp_p_ws.jj_nh4_b,'g','linew',2); 
        plot(obm_jj_nh4_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_nh4_b),obm_jj_nh4_b./14 + obs_std_jj_nh4_b./14,'m-','linew',2);
        plot(1:length(obm_jj_nh4_b),obm_jj_nh4_b./14 - obs_std_jj_nh4_b./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4 bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NH4 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_bot_jj'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_chl_b,'b','linew',2);
        plot(pcase.jj_chl_b,'k','linew',2);
        plot(mp_p_ws.jj_chl_b,'g','linew',2); 
        plot(obm_jj_chl_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_chl_b),obm_jj_chl_b + obs_std_jj_chl_b,'m-','linew',2);
        plot(1:length(obm_jj_chl_b),obm_jj_chl_b - obs_std_jj_chl_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL chl bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('chl (ug/L)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_bot_jj'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_temp_b,'b','linew',2);
        plot(pcase.jj_temp_b,'k','linew',2);
        plot(mp_p_ws.jj_temp_b,'g','linew',2); 
        plot(obm_jj_temp_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_temp_b),obm_jj_temp_b + obs_std_jj_temp_b,'m-','linew',2);
        plot(1:length(obm_jj_temp_b),obm_jj_temp_b - obs_std_jj_temp_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL temp bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('temp (o^C)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_bot_jj'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_salt_b,'b','linew',2);
        plot(pcase.jj_salt_b,'k','linew',2);
        plot(mp_p_ws.jj_salt_b,'g','linew',2); 
        plot(obm_jj_salt_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_salt_b),obm_jj_salt_b + obs_std_jj_salt_b,'m-','linew',2);
        plot(1:length(obm_jj_salt_b),obm_jj_salt_b - obs_std_jj_salt_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL salt bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('salt (psu)','fontsize',16)
        xticks([0:50:366])
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_bot_jj'),'-dpng')    


