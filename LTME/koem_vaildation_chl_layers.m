 clearvars clim2*
clim2=clim1.obm_gy_chl_ft_1(:,3);
clim2_std = clim1.obs_std_gy_chl_ft_1(:,3);
clim2_2=clim1.obm_gy_chl_ft_2(:,3);
clim2_2_std = clim1.obs_std_gy_chl_ft_2(:,3);
 clearvars clim2*
clim2=clim1.obm_gy_chl_ft_1(:,3);
clim2_std = clim1.obs_std_gy_chl_ft_1(:,3);
clim2_2=clim1.obm_gy_chl_ft_2(:,3);
clim2_2_std = clim1.obs_std_gy_chl_ft_2(:,3);

    fig = figure; hold on;
    for i = 1:20
            plot(regime1_kdc_zerosw.gy_chl_mid(i,:),'linew',2);
    end
%             plot(nanmean(regime2_kdc_zerosw.gy_chl,1),'g','linew',2);

             % 1regime
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2 + clim2_std;
            lower_bound_plt = clim2 - clim2_std;
            for i = 1:4
            patch([plt_filter_1{i} fliplr(plt_filter_1{i})], [upper_bound_plt(plt_filter_1{i})' fliplr(lower_bound_plt(plt_filter_1{i}))'],'r');
            plot(plt_filter_1{i},clim2(plt_filter_1{i}),'r-','linew',2);
            plot(plt_filter_1{i},upper_bound_plt(plt_filter_1{i}),'r-','linew',2);
            plot(plt_filter_1{i},lower_bound_plt(plt_filter_1{i}),'r-','linew',2);
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
            patch([plt_filter_2{i} fliplr(plt_filter_2{i})], [upper_bound_plt(plt_filter_2{i})' fliplr(lower_bound_plt(plt_filter_2{i}))'],'g');
            plot(plt_filter_2{i},clim2_2(plt_filter_2{i}),'g-','linew',2);
            plot(plt_filter_2{i},upper_bound_plt(plt_filter_2{i}),'g-','linew',2);
            plot(plt_filter_2{i},lower_bound_plt(plt_filter_2{i}),'g-','linew',2);
            end 

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL chl']);
            xlabel('time','fontsize',13)
            ylabel('chl (umol/m^3)','fontsize',13)
            grid on
            ylim([0 30])
            xlim(xlim_set)
            set(gca,'fontsize',13)
