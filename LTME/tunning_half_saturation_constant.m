
close all;clear;clc;
ten_1st = load('half_saturation_setting_10yr_1st.mat');
ten_2nd = load('half_saturation_setting_10yr_2nd.mat');

%% tunn1~2
%surf    
    %fix 1st
    Kn = 1;
    Ka = 1;
    Kp = 0.025;
    
    lim_n = (ten_1st.plt_no3 ./ (Kn + ten_1st.plt_no3)) .* (1./(1 + ten_1st.plt_nh4)./ Ka) + (ten_1st.plt_nh4./(Ka + ten_1st.plt_nh4));
    lim_p = (ten_1st.plt_po4 ./ (Kp + ten_1st.plt_po4));
    figure; plot(lim_n,'linew',2); hold on; plot(lim_p,'linew',2,'color','r');xlim([1 365]);
    alpha(.3);  grid on; 


    %fix 2nd
    lim_n = (ten_2nd.plt_no3 ./ (Kn + ten_2nd.plt_no3)) .* (1./(1 + ten_2nd.plt_nh4)./ Ka) + (ten_2nd.plt_nh4./(Ka + ten_2nd.plt_nh4));
    lim_p = (ten_2nd.plt_po4 ./ (Kp + ten_2nd.plt_po4));
    plot(lim_n,'linew',2,'color','g'); hold on; plot(lim_p,'linew',2,'color','c');xlim([1 365]);
    alpha(.3);  grid on; 
    
    legend('N_1_s_t','P_1_s_t','N_2_n_d','P_2_n_d');  set(gca,'fontsize',13); ylim([0.5 1.1]);

%bot
%default 1st
    Kn = 1;
    Ka = 1;
    Kp = 0.025;
    
    lim_n = (ten_1st.plt_no3_b ./ (Kn + ten_1st.plt_no3_b)) .* (1./(1 + ten_1st.plt_nh4_b)./ Ka) + (ten_1st.plt_nh4_b./(Ka + ten_1st.plt_nh4_b));
    lim_p = (ten_1st.plt_po4_b ./ (Kp + ten_1st.plt_po4_b));
    figure; plot(lim_n,'linew',2); hold on; plot(lim_p,'linew',2,'color','r');xlim([1 365]);
    alpha(.3);  grid on;
   

%default 2nd
    lim_n = (ten_2nd.plt_no3_b ./ (Kn + ten_2nd.plt_no3_b)) .* (1./(1 + ten_2nd.plt_nh4_b)./ Ka) + (ten_2nd.plt_nh4_b./(Ka + ten_2nd.plt_nh4_b));
    lim_p = (ten_2nd.plt_po4_b ./ (Kp + ten_2nd.plt_po4_b));
    plot(lim_n,'linew',2,'color','g'); hold on; plot(lim_p,'linew',2,'color','c');xlim([1 365]);
    alpha(.3);  grid on;
    legend('N_1_s_t','P_1_s_t','N_2_n_d','P_2_n_d'); set(gca,'fontsize',13);   ylim([0.5 1.1]);



%fix 1st
    Kn = 1;
    Ka = 1;
    Kp = 0.038;
    
    lim_n = (ten_1st.plt_no3 ./ (Kn + ten_1st.plt_no3)) .* (1./(1 + ten_1st.plt_nh4)./ Ka) + (ten_1st.plt_nh4./(Ka + ten_1st.plt_nh4));
    lim_p = (ten_1st.plt_po4 ./ (Kp + ten_1st.plt_po4));
    figure; plot(lim_n,'linew',2); hold on; plot(lim_p,'linew',2,'color','r');xlim([1 365]);
    alpha(.3);  grid on;
   

%fix 2nd
    lim_n = (ten_2nd.plt_no3 ./ (Kn + ten_2nd.plt_no3)) .* (1./(1 + ten_2nd.plt_nh4)./ Ka) + (ten_2nd.plt_nh4./(Ka + ten_2nd.plt_nh4));
    lim_p = (ten_2nd.plt_po4 ./ (Kp + ten_2nd.plt_po4));
    plot(lim_n,'linew',2,'color','g'); hold on; plot(lim_p,'linew',2,'color','c');xlim([1 365]);
    alpha(.3);  grid on; 
    
    legend('N_1_s_t','P_1_s_t','N_2_e_d','P_2_e_d');
    
    
    
%% default
%surf
%default 1st
    Kn = 0.8;
    Ka = 0.8;
    Kp = 16.7^-1;
    
    lim_n = (ten_1st.plt_no3 ./ (Kn + ten_1st.plt_no3)) .* (1./(1 + ten_1st.plt_nh4)./ Ka) + (ten_1st.plt_nh4./(Ka + ten_1st.plt_nh4));
    lim_p = (ten_1st.plt_po4 ./ (Kp + ten_1st.plt_po4));
    figure; plot(lim_n,'linew',2); hold on; plot(lim_p,'linew',2,'color','r');xlim([1 365]);
    alpha(.3);  grid on;
   

%default 2nd
    lim_n = (ten_2nd.plt_no3 ./ (Kn + ten_2nd.plt_no3)) .* (1./(1 + ten_2nd.plt_nh4)./ Ka) + (ten_2nd.plt_nh4./(Ka + ten_2nd.plt_nh4));
    lim_p = (ten_2nd.plt_po4 ./ (Kp + ten_2nd.plt_po4));
    plot(lim_n,'linew',2,'color','g'); hold on; plot(lim_p,'linew',2,'color','c');xlim([1 365]);
    alpha(.3);  grid on;
    legend('N_1_s_t','P_1_s_t','N_2_n_d','P_2_n_d'); set(gca,'fontsize',13);ylim([0.5 1.1]);

    
 %bot   
 %default 1st
    Kn = 0.8;
    Ka = 0.8;
    Kp = 16.7^-1;
    
    lim_n = (ten_1st.plt_no3_b ./ (Kn + ten_1st.plt_no3_b)) .* (1./(1 + ten_1st.plt_nh4_b)./ Ka) + (ten_1st.plt_nh4_b./(Ka + ten_1st.plt_nh4_b));
    lim_p = (ten_1st.plt_po4_b ./ (Kp + ten_1st.plt_po4_b));
    figure; plot(lim_n,'linew',2); hold on; plot(lim_p,'linew',2,'color','r');xlim([1 365]);
    alpha(.3);  grid on;
   

%default 2nd
    lim_n = (ten_2nd.plt_no3_b ./ (Kn + ten_2nd.plt_no3_b)) .* (1./(1 + ten_2nd.plt_nh4_b)./ Ka) + (ten_2nd.plt_nh4_b./(Ka + ten_2nd.plt_nh4_b));
    lim_p = (ten_2nd.plt_po4_b ./ (Kp + ten_2nd.plt_po4_b));
    plot(lim_n,'linew',2,'color','g'); hold on; plot(lim_p,'linew',2,'color','c');xlim([1 365]);
    alpha(.3);  grid on;
    legend('N_1_s_t','P_1_s_t','N_2_n_d','P_2_n_d'); set(gca,'fontsize',13);   ylim([0.5 1.1]);

    
    
    
    
    
    