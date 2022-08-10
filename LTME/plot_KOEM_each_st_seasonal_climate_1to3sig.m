close all; clear all; close all;

% from : plot_KOEM_st_seasonal_climate_3regime_v5_spatial_sigma_extract_10y
koem_1sig=load(['koem_climate_10yr_v5_1sig_gy_each_st.mat']); % full (9point) data (extract_sigma_for_each_st)
koem_2sig=load(['koem_climate_10yr_v5_2sig_gy_each_st.mat']); % full (9point) data (extract_sigma_for_each_st)
koem_3sig=load(['koem_climate_10yr_v5_3sig_gy_each_st.mat']); % full (9point) data (extract_sigma_for_each_st)

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