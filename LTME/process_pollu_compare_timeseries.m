clc;clear all;close all

raw=xlsread('광양만 자료 송부_ed_ed.xlsx','데이터');
raw2=xlsread('광양만 자료 송부_ed_ed.xlsx','섬진강 하동 유량, 수질');
 % 1-year 3-month 4-discharge 5-BOD 6-COD 7-SS 8-TN 9-TP
for i=1:11
    data(i,:,:)=raw(1+(i-1)*108:108*i,:);  % year, data , 1~9 
end

for i=1:11
    data_sum(i,:,:)=raw2(1+(i-1)*12:12*i,:);
end

meandata=squeeze(nanmean(data));
meandatas=squeeze(nanmean(data_sum));

% 1.광양산업단지폐수종말처리장
% 2.광양중앙하수종말처리장
% 3.광양하수종말처리장
% 4.광영하수종말처리장
% 5.여수월내산업단지폐수종말처리장
% 6.여수율촌산업단지폐수종말처리장
% 7.여수중흥산업단지폐수종말처리장
% 8.여수하수종말처리장
% 9.진월하수종말처리장
clearvars data_st
for i=1:11
    for j = 1:9
        data_st(i,:,:,j)= data(i,1+(12*(j-1)):12*j,:); % year, month, variable, station
    end
end

for i = 1:9
    data_st_clim(:,:,i) = meandata(1+(12*(i-1)):12*i,:);
end
% data_st_clim = mon, variable, station


load('gy_pollution.mat')


clearvars com_t_tic
com_t_tic = {'2007.01','02','03','04','05','06','07','08','09','10','11','12', ...
            '2008.01','02','03','04','05','06','07','08','09','10','11','12', ...
            '2009.01','02','03','04','05','06','07','08','09','10','11','12', ...
            '2010.01','02','03','04','05','06','07','08','09','10','11','12', ...
            '2011.01','02','03','04','05','06','07','08','09','10','11','12', ...
            '2012.01','02','03','04','05','06','07','08','09','10','11','12', ...
            '2013.01','02','03','04','05','06','07','08','09','10','11','12', ...
            '2014.01','02','03','04','05','06','07','08','09','10','11','12', ...
            '2015.01','02','03','04','05','06','07','08','09','10','11','12', ...
            '2016.01','02','03','04','05','06','07','08','09','10','11','12', ...
            '2017.01','02','03','04','05','06','07','08','09','10','11','12'};
clearvars com_t_tic
com_t_tic = {'2007', ...
            '2008', ...
            '2009', ...
            '2010', ...
            '2011', ...
            '2012', ...
            '2013', ...
            '2014', ...
            '2015', ...
            '2016', ...
            '2017'};        
        


% TN 
for i = 1:9
clearvars re re_clim
fig = figure; hold on;
re = reshape(squeeze(data_st(:,:,8,i))',11*12,1); % TN for st. i timeseries     
re_clim=repmat(data_st_clim(:,8,i),11,1); % TN for st. i climate    
nnan_dx =find(isnan(re)==0);
re_clim(1:nnan_dx(1)-1)=NaN;

%         plot(zeros(365,1),'g--','linew',2); % zero lines
%% 2001
        plot(re .*(1000/14),'b','linew',2);
        plot(re_clim .*(1000/14),'r','linew',2);
        
        tn_sewer(i,:)=re;
%         clearvars *_bound_plt nonan_data_plt discon_p
%         nonan_data_plt = find(isnan(c2001.obm_gy_no3)==0);
%         nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%         discon_p = find(nonan_diff ~= 1); % find discontinuous point
%         discon_p(4) = length(nonan_data_plt); %end point
%         upper_bound_plt = c2001.obm_gy_no3./14 + c2001.obs_std_gy_no3./14;
%         lower_bound_plt = c2001.obm_gy_no3./14 - c2001.obs_std_gy_no3./14;
%         for i = 1:4
%         if i == 1
%         patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
%         else
%         patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
%         end
%         end 
%         plot(c2001.obm_gy_no3./14,'r-','linew',2); 
%         plot(1:length(c2001.obm_gy_no3),upper_bound_plt,'r-','linew',2);
%         plot(1:length(c2001.obm_gy_no3),lower_bound_plt,'r-','linew',2);
            
        alpha(0.3) %transparency      
        
        title(['GY Sewer TN timeseries vs. climate st.',num2str(i)]);
        xlabel('time','fontsize',13)
        ylabel('TN (umol/m^3)','fontsize',13)
        grid on
        ylim([0 2600])
        xlim([1 length(re)])
        xticks(1:12:length(re))
        xticklabels(com_t_tic)
        xtickangle(60)
        set(gca,'fontsize',13)
%         print(fig,strcat('sewer_tn_timeseriese_vs_climate_st_',num2str(i)),'-dpng')
close
end
        
% TP
for i = 1:9
clearvars re re_clim
fig = figure; hold on;
re = reshape(squeeze(data_st(:,:,9,i))',11*12,1); % TN for st. i timeseries     
re_clim=repmat(data_st_clim(:,9,i),11,1); % TN for st. i climate    
nnan_dx =find(isnan(re)==0);
re_clim(1:nnan_dx(1)-1)=NaN;

%         plot(zeros(365,1),'g--','linew',2); % zero lines
%% 2001
        plot(re .*(1000/30.973762),'b','linew',2);
        plot(re_clim .*(1000/30.973762),'r','linew',2);
        tp_sewer(i,:)=re;
%         clearvars *_bound_plt nonan_data_plt discon_p
%         nonan_data_plt = find(isnan(c2001.obm_gy_no3)==0);
%         nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%         discon_p = find(nonan_diff ~= 1); % find discontinuous point
%         discon_p(4) = length(nonan_data_plt); %end point
%         upper_bound_plt = c2001.obm_gy_no3./14 + c2001.obs_std_gy_no3./14;
%         lower_bound_plt = c2001.obm_gy_no3./14 - c2001.obs_std_gy_no3./14;
%         for i = 1:4
%         if i == 1
%         patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
%         else
%         patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
%         end
%         end 
%         plot(c2001.obm_gy_no3./14,'r-','linew',2); 
%         plot(1:length(c2001.obm_gy_no3),upper_bound_plt,'r-','linew',2);
%         plot(1:length(c2001.obm_gy_no3),lower_bound_plt,'r-','linew',2);
            
        alpha(0.3) %transparency      
        
        title(['GY Sewer TP timeseries vs. climate st.',num2str(i)]);
        xlabel('time','fontsize',13)
        ylabel('TP (umol/m^3)','fontsize',13)
        grid on
        ylim([0 140])
        xlim([1 length(re)])
        xticks(1:12:length(re))
        xticklabels(com_t_tic)
        xtickangle(60)
        set(gca,'fontsize',13)
%         print(fig,strcat('sewer_tp_timeseriese_vs_climate_st_',num2str(i)),'-dpng')
close
        end

%% total mass
% mass (discharge)
for i = 1:9
clearvars re re_clim
fig = figure; hold on;
re = reshape(squeeze(data_st(:,:,4,i))',11*12,1); % TN for st. i timeseries     
re_clim=repmat(data_st_clim(:,4,i),11,1); % TN for st. i climate    
nnan_dx =find(isnan(re)==0);
re_clim(1:nnan_dx(1)-1)=NaN;

%         plot(zeros(365,1),'g--','linew',2); % zero lines
%% 2001
        plot(re ,'b','linew',2);
        plot(re_clim ,'r','linew',2);
        dis_sewer(i,:)=re;
%         clearvars *_bound_plt nonan_data_plt discon_p
%         nonan_data_plt = find(isnan(c2001.obm_gy_no3)==0);
%         nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%         discon_p = find(nonan_diff ~= 1); % find discontinuous point
%         discon_p(4) = length(nonan_data_plt); %end point
%         upper_bound_plt = c2001.obm_gy_no3./14 + c2001.obs_std_gy_no3./14;
%         lower_bound_plt = c2001.obm_gy_no3./14 - c2001.obs_std_gy_no3./14;
%         for i = 1:4
%         if i == 1
%         patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
%         else
%         patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
%         end
%         end 
%         plot(c2001.obm_gy_no3./14,'r-','linew',2); 
%         plot(1:length(c2001.obm_gy_no3),upper_bound_plt,'r-','linew',2);
%         plot(1:length(c2001.obm_gy_no3),lower_bound_plt,'r-','linew',2);
            
        alpha(0.3) %transparency      
        
        title(['GY Sewer TN timeseries vs. climate st.',num2str(i)]);
        xlabel('time','fontsize',13)
        ylabel('TN (umol/m^3)','fontsize',13)
        grid on
        ylim([-inf inf])
        xlim([1 length(re)])
        xticks(1:12:length(re))
        xticklabels(com_t_tic)
        xtickangle(60)
        set(gca,'fontsize',13)
%         print(fig,strcat('sewer_tn_timeseriese_vs_climate_st_',num2str(i)),'-dpng')
end

save('sewer_from_LTME.mat','*_sewer');

% TN 
for i = 1:9
clearvars re re_clim
fig = figure; hold on;
re = reshape(squeeze(data_st(:,:,8,i))',11*12,1); % TN for st. i timeseries     
re_clim=repmat(data_st_clim(:,8,i),11,1); % TN for st. i climate    
nnan_dx =find(isnan(re)==0);
re_clim(1:nnan_dx(1)-1)=NaN;

%         plot(zeros(365,1),'g--','linew',2); % zero lines
%% 2001
        plot(re .*(1000/14),'b','linew',2);
        plot(re_clim .*(1000/14),'r','linew',2);
%         clearvars *_bound_plt nonan_data_plt discon_p
%         nonan_data_plt = find(isnan(c2001.obm_gy_no3)==0);
%         nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%         discon_p = find(nonan_diff ~= 1); % find discontinuous point
%         discon_p(4) = length(nonan_data_plt); %end point
%         upper_bound_plt = c2001.obm_gy_no3./14 + c2001.obs_std_gy_no3./14;
%         lower_bound_plt = c2001.obm_gy_no3./14 - c2001.obs_std_gy_no3./14;
%         for i = 1:4
%         if i == 1
%         patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
%         else
%         patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
%         end
%         end 
%         plot(c2001.obm_gy_no3./14,'r-','linew',2); 
%         plot(1:length(c2001.obm_gy_no3),upper_bound_plt,'r-','linew',2);
%         plot(1:length(c2001.obm_gy_no3),lower_bound_plt,'r-','linew',2);
            
        alpha(0.3) %transparency      
        
        title(['GY Sewer TN timeseries vs. climate st.',num2str(i)]);
        xlabel('time','fontsize',13)
        ylabel('TN (umol/m^3)','fontsize',13)
        grid on
        ylim([0 2600])
        xlim([1 length(re)])
        xticks(1:12:length(re))
        xticklabels(com_t_tic)
        xtickangle(60)
        set(gca,'fontsize',13)
        print(fig,strcat('sewer_tn_timeseriese_vs_climate_st_',num2str(i)),'-dpng')
end
        
% TP
for i = 1:9
clearvars re re_clim
fig = figure; hold on;
re = reshape(squeeze(data_st(:,:,9,i))',11*12,1); % TN for st. i timeseries     
re_dis = reshape(squeeze(data_st(:,:,4,i))',11*12,1); % TN for st. i timeseries 
re_clim=repmat(data_st_clim(:,9,i),11,1); % TN for st. i climate    
re_clim_dis=repmat(data_st_clim(:,4,i),11,1); % TN for st. i climate    

nnan_dx =find(isnan(re)==0);
re_clim(1:nnan_dx(1)-1)=NaN;

%         plot(zeros(365,1),'g--','linew',2); % zero lines
%% 2001
        plot(re .*(1000/30.973762) ,'b','linew',2);
        plot(re_clim .*(1000/30.973762),'r','linew',2);
%         clearvars *_bound_plt nonan_data_plt discon_p
%         nonan_data_plt = find(isnan(c2001.obm_gy_no3)==0);
%         nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%         discon_p = find(nonan_diff ~= 1); % find discontinuous point
%         discon_p(4) = length(nonan_data_plt); %end point
%         upper_bound_plt = c2001.obm_gy_no3./14 + c2001.obs_std_gy_no3./14;
%         lower_bound_plt = c2001.obm_gy_no3./14 - c2001.obs_std_gy_no3./14;
%         for i = 1:4
%         if i == 1
%         patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
%         else
%         patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
%         end
%         end 
%         plot(c2001.obm_gy_no3./14,'r-','linew',2); 
%         plot(1:length(c2001.obm_gy_no3),upper_bound_plt,'r-','linew',2);
%         plot(1:length(c2001.obm_gy_no3),lower_bound_plt,'r-','linew',2);
            
        alpha(0.3) %transparency      
        
        title(['GY Sewer TP timeseries vs. climate st.',num2str(i)]);
        xlabel('time','fontsize',13)
        ylabel('TP (umol/m^3)','fontsize',13)
        grid on
        ylim([0 140])
        xlim([1 length(re)])
        xticks(1:12:length(re))
        xticklabels(com_t_tic)
        xtickangle(60)
        set(gca,'fontsize',13)
        print(fig,strcat('sewer_tp_timeseriese_vs_climate_st_',num2str(i)),'-dpng')
        end



xtick




re = reshape(squeeze(data_st(:,:,8,1))',11*12,1); % TN for st. 1 timeseries
plot(re)



