   close all; clear; clc;
   
   cd F:\장기생태_output_임시
   c1997=load('1997_sewer_re_v9grid_t1.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*');
   c2007=load('2007_sewer_re_v9grid_t1.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*');
   
   
    plt_no3 = c1997.gy_no3;
    plt_nh4 = c1997.gy_nh4;
    plt_po4 = c1997.gy_po4;
    
    plt_no3_b = c1997.gy_no3_b;
    plt_nh4_b = c1997.gy_nh4_b;
    plt_po4_b = c1997.gy_po4_b;
    
    plt_no3_2 = c2007.gy_no3;
    plt_nh4_2 = c2007.gy_nh4;
    plt_po4_2 = c2007.gy_po4;
    
    plt_no3_b_2 = c2007.gy_no3_b;
    plt_nh4_b_2 = c2007.gy_nh4_b;
    plt_po4_b_2 = c2007.gy_po4_b;
   
%     save('half_saturation_setting_10yr_1st.mat','plt_*');
    
% NH4 vs. PO4
    np_ratio = plt_nh4 ./ plt_po4;
    np_ratio_b = plt_nh4_b ./ plt_po4_b;
    
    figure; plot(np_ratio,'linew',2); hold on; plot(np_ratio_b,'linew',2,'color','r');
    plot(1:365,repmat(16,1,365),'g--','linew',2) ;xlim([1 365]); alpha(.3);  grid on; 
    legend('np-ratio','np-ratio(bot)'); set(gca,'fontsize',15);
    
    np_ratio_2 = plt_nh4_2 ./ plt_po4_2;
    np_ratio_b_2 = plt_nh4_b_2 ./ plt_po4_b_2;
    
   figure; plot(np_ratio_2,'linew',2); hold on; plot(np_ratio_b_2,'linew',2,'color','r');
    plot(1:365,repmat(16,1,365),'g--','linew',2) ;xlim([1 365]); alpha(.3);  grid on; 
    legend('np-ratio','np-ratio(bot)'); set(gca,'fontsize',15);

% NH4+NO3 vs. PO4    
    np_ratio = (plt_nh4 + plt_no3) ./ plt_po4;
    np_ratio_b = (plt_nh4_b + plt_no3_b) ./ plt_po4_b;
    
    figure; plot(np_ratio,'linew',2); hold on; plot(np_ratio_b,'linew',2,'color','r');
    plot(1:365,repmat(16,1,365),'g--','linew',2) ;xlim([1 365]); alpha(.3);  grid on; 
    legend('np-ratio','np-ratio(bot)'); set(gca,'fontsize',15);

    np_ratio_2 = (plt_nh4_2 + plt_no3_2) ./ plt_po4_2;
    np_ratio_b_2 = (plt_nh4_b_2 + plt_no3_b_2) ./ plt_po4_b_2;
    
    figure; plot(np_ratio_2,'linew',2); hold on; plot(np_ratio_b_2,'linew',2,'color','r');
    plot(1:365,repmat(16,1,365),'g--','linew',2) ;xlim([1 365]); alpha(.3);  grid on; 
    legend('np-ratio','np-ratio(bot)'); set(gca,'fontsize',15);

    
%default
    Kn = 0.6;
    Ka = 0.6;
    Kp = 0.0599;
    
    lim_n = (plt_no3 ./ (Kn + plt_no3)) .* (1./(1 + plt_nh4)./ Ka) + (plt_nh4./(Ka + plt_nh4));
    lim_p = (plt_po4 ./ (Kp + plt_po4));
    figure; plot(lim_n,'linew',2); hold on; plot(lim_p,'linew',2,'color','r');xlim([1 365]);
    alpha(.3);  grid on; 
    
    Kn = 0.30;
    Ka = 0.30;
    Kp = 0.0187;
    
    lim_n = (plt_no3 ./ (Kn + plt_no3)) .* (1./(1 + plt_nh4)./ Ka) + (plt_nh4./(Ka + plt_nh4));
    lim_p = (plt_po4 ./ (Kp + plt_po4));
    figure; plot(lim_n,'linew',2); hold on; plot(lim_p,'linew',2,'color','r');xlim([1 365]);
    alpha(.3);  grid on; 
    
    Kn = 1.2;
    Ka = 1.2;
    Kp = 0.1;
    
    lim_n = (plt_no3 ./ (Kn + plt_no3)) .* (1./(1 + plt_nh4)./ Ka) + (plt_nh4./(Ka + plt_nh4));
    lim_p = (plt_po4 ./ (Kp + plt_po4));
    figure; plot(lim_n,'linew',2); hold on; plot(lim_p,'linew',2,'color','r');xlim([1 365]);
    alpha(.3);  grid on;
    
    
    % add po4
    Kn = 0.80;
    Ka = 0.80;
    Kp = 0.0599;
    
    lim_n = (plt_no3 ./ (Kn + plt_no3)) .* (1./(1 + plt_nh4)./ Ka) + (plt_nh4./(Ka + plt_nh4));
    lim_p = ((plt_po4+1) ./ (Kp + (plt_po4+1)));
    figure; plot(lim_n,'linew',2); hold on; plot(lim_p,'linew',2,'color','r');xlim([1 365]);
    alpha(.3);  grid on; 
    
    figure; plot(lim_n - lim_p,'linew',2); hold on; 
    
    Kn = 0.80;
    Ka = 0.80;
    Kp = 0.0599;
    
    lim_n = (plt_no3 ./ (Kn + plt_no3)) .* (1./(1 + plt_nh4)./ Ka) + (plt_nh4./(Ka + plt_nh4));
    lim_p = (plt_po4 ./ (Kp + plt_po4));
    figure; plot(lim_n,'linew',2); hold on; plot(lim_p,'linew',2,'color','r');xlim([1 365]);
    alpha(.3);  grid on; 
    
    Kn = 0.80;
    Ka = 0.80;
    Kp = 0.000001;
    lim_n = (plt_no3 ./ (Kn + plt_no3)) .* (1./(1 + plt_nh4)./ Ka) + (plt_nh4./(Ka + plt_nh4));
    lim_p = (plt_po4 ./ (Kp + plt_po4));
    figure; plot(lim_n - lim_p,'linew',2); hold on; xlim([1 365]);
    alpha(.3);  grid on;
    
    figure; plot(lim_n,'linew',2); hold on; plot(lim_p,'linew',2,'color','r');xlim([1 365]);
    alpha(.3);  grid on; 
    
    
%     Kn = 0.80^-1;
%     Ka = 0.80^-1;
%     Kp = 0.60^-1;
    Kn = 0.50;
    Ka = 0.50;
    Kp = 0.28;
    
    lim_n = (plt_no3 ./ (Kn + plt_no3)) .* (1./(1 + plt_nh4)./ Ka) + (plt_nh4./(Ka + plt_nh4));
    lim_p = (plt_po4 ./ (Kp + plt_po4));
    plot(lim_n,'b--','linew',2); hold on; plot(lim_p,'r--','linew',2);
    yticks(0.1:.1:1.5); grid on; xlim([1 365]); alpha(.3); ylim([0 1.5]);
    legend('0.8 N-uptake','0.6 P-uptake','0.5 N-uptake','0.28 P-uptake');
    
    Kn = 0.80;
    Ka = 0.50;
    Kp = 0.20;
    
    lim_n = (plt_no3 ./ (Kn + plt_no3)) .* (1./(1 + plt_nh4)./ Ka) + (plt_nh4./(Ka + plt_nh4));
    lim_p = (plt_po4 ./ (Kp + plt_po4));
    plot(lim_n,'b-o'); hold on; plot(lim_p,'r-o');

% bot
    Kn = 0.80;
    Ka = 0.80;
    Kp = 0.60;
    
    lim_n_b = (plt_no3_b ./ (Kn + plt_no3_b)) .* (1./(1 + plt_nh4_b)./ Ka) + (plt_nh4_b./(Ka + plt_nh4_b));
    lim_p_b = (plt_po4_b ./ (Kp + plt_po4_b));
    figure; plot(lim_n_b,'linew',2); hold on; plot(lim_p_b,'linew',2,'color','r');
    
%     Kn = 0.80^-1;
%     Ka = 0.80^-1;
%     Kp = 0.60^-1;
    Kn = 0.50;
    Ka = 0.50;
    Kp = 0.28;
    
    lim_n_b = (plt_no3_b ./ (Kn + plt_no3_b)) .* (1./(1 + plt_nh4_b)./ Ka) + (plt_nh4_b./(Ka + plt_nh4_b));
    lim_p_b = (plt_po4_b ./ (Kp + plt_po4_b));
    plot(lim_n_b,'b--','linew',2); hold on; plot(lim_p_b,'r--','linew',2);
    yticks(0.1:.1:1.5); grid on; xlim([1 365]); alpha(.3); ylim([0 1.5]);
    legend('0.8 N-uptake','0.6 P-uptake','0.5 N-uptake','0.28 P-uptake');