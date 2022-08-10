close all; clear all; clc;

cd D:\장기생태\Dynamic\06_river\data\sj_1980to1996
%river discharge 1980~1989, calculated from water lv. & discharge relation
%from ykang
sj = load('songjung_climate_days_bio_to07to16.mat'); % 3regime climate (1997~2006, 2007~2015, 2016~2019);
ng = load('namgang_climate_days_bio_to07to16.mat');  % 3regime climate (1997~2006, 2007~2015, 2016~2019);

cd D:\장기생태\Dynamic\06_river\하수종말처리장\여수_중흥\
load('sewer_monthly_climate_fix_to06to15.mat');
load('D:\장기생태\Dynamic\06_river\data\sj_1980to1996\songjung_climate_days_load_to07to16.mat','sj_*_load_*');

cycle = 365;
%if it's not leap yr.. delete 2/29
if cycle == 365 
%     sumjin_tn.yp_w_tn_04(60) =[];
%     nam_tn.yp_w_tn_04(60)=[];
            %dis
        clim_dis_g_d(60)=[];
        clim_dis_gc_d(60)=[];
        clim_dis_gh_d(60)=[];
        clim_dis_gh2_d(60)=[];
        clim_dis_d(60)=[];
        clim_dis_y_d(60)=[];
        clim_dis_j_d(60)=[];
        clim_dis_ye_d(60)=[];
        clim_dis_jw_d(60)=[];
        %TN
        clim_tn_g_ns_d(60)=[];
        clim_tn_gc_ns_d(60)=[];
        clim_tn_gh_ns_d(60)=[];
        clim_tn_gh2_ns_d(60)=[];
        clim_tn_w_1st_ns_d(60)=[];
        clim_tn_w_2nd_ns_d(60)=[];
        clim_tn_w_3rd_ns_d(60)=[];
        clim_tn_y_ns_d(60)=[];
        clim_tn_j_1st_ns_d(60)=[];
        clim_tn_j_2nd_ns_d(60)=[];
        clim_tn_j_3rd_ns_d(60)=[];
        clim_tn_ye_ns_d(60)=[];
        clim_tn_jw_ns_d(60)=[];
        %TP
        clim_tp_g_ns_d(60)=[];
        clim_tp_gc_ns_d(60)=[];
        clim_tp_gh_ns_d(60)=[];
        clim_tp_gh2_ns_d(60)=[];
        clim_tp_w_1st_ns_d(60)=[];
        clim_tp_w_2nd_ns_d(60)=[];
        clim_tp_w_3rd_ns_d(60)=[];
        clim_tp_y_ns_d(60)=[];
        clim_tp_j_1st_ns_d(60)=[];
        clim_tp_j_2nd_ns_d(60)=[];
        clim_tp_j_3rd_ns_d(60)=[];
        clim_tp_ye_ns_d(60)=[];
        clim_tp_jw_ns_d(60)=[];
end


%--- 광양산업단지폐수종말처리장----------------------------------------------
% MW_N = 14.006720;
total_dis_tn=(clim_dis_g_d/86400 .* clim_tn_g_ns_d) + (clim_dis_gc_d/86400 .*  clim_tn_gc_ns_d) + ...
    (clim_dis_gh_d/86400 .* clim_tn_gh_ns_d) + (clim_dis_gh2_d/86400 .* clim_tn_gh2_ns_d) + (clim_dis_d/86400 .* clim_tn_w_1st_ns_d) +...
    (clim_dis_y_d/86400 .* clim_tn_y_ns_d) + (clim_dis_j_d/86400 .* clim_tn_j_1st_ns_d) + (clim_dis_ye_d/86400 .* clim_tn_ye_ns_d) + ...
    (clim_dis_jw_d/86400 .* clim_tn_jw_ns_d);

total_dis_tn_2nd=(clim_dis_g_d/86400 .* clim_tn_g_ns_d) + (clim_dis_gc_d/86400 .*  clim_tn_gc_ns_d) + ...
    (clim_dis_gh_d/86400 .* clim_tn_gh_ns_d) + (clim_dis_gh2_d/86400 .* clim_tn_gh2_ns_d) + (clim_dis_d/86400 .* clim_tn_w_2nd_ns_d) +...
    (clim_dis_y_d/86400 .* clim_tn_y_ns_d) + (clim_dis_j_d/86400 .* clim_tn_j_2nd_ns_d) + (clim_dis_ye_d/86400 .* clim_tn_ye_ns_d) + ...
    (clim_dis_jw_d/86400 .* clim_tn_jw_ns_d);

total_dis_tn_3rd=(clim_dis_g_d/86400 .* clim_tn_g_ns_d) + (clim_dis_gc_d/86400 .*  clim_tn_gc_ns_d) + ...
    (clim_dis_gh_d/86400 .* clim_tn_gh_ns_d) + (clim_dis_gh2_d/86400 .* clim_tn_gh2_ns_d) + (clim_dis_d/86400 .* clim_tn_w_3rd_ns_d) +...
    (clim_dis_y_d/86400 .* clim_tn_y_ns_d) + (clim_dis_j_d/86400 .* clim_tn_j_3rd_ns_d) + (clim_dis_ye_d/86400 .* clim_tn_ye_ns_d) + ...
    (clim_dis_jw_d/86400 .* clim_tn_jw_ns_d);

total_dis_tp=(clim_dis_g_d/86400 .* clim_tp_g_ns_d) + (clim_dis_gc_d/86400 .*  clim_tp_gc_ns_d) + ...
    (clim_dis_gh_d/86400 .* clim_tp_gh_ns_d) + (clim_dis_gh2_d/86400 .* clim_tp_gh2_ns_d) + (clim_dis_d/86400 .* clim_tp_w_1st_ns_d) +...
    (clim_dis_y_d/86400 .* clim_tp_y_ns_d) + (clim_dis_j_d/86400 .* clim_tp_j_1st_ns_d) + (clim_dis_ye_d/86400 .* clim_tp_ye_ns_d) + ...
    (clim_dis_jw_d/86400 .* clim_tp_jw_ns_d);

total_dis_tp_2nd=(clim_dis_g_d/86400 .* clim_tp_g_ns_d) + (clim_dis_gc_d/86400 .*  clim_tp_gc_ns_d) + ...
    (clim_dis_gh_d/86400 .* clim_tp_gh_ns_d) + (clim_dis_gh2_d/86400 .* clim_tp_gh2_ns_d) + (clim_dis_d/86400 .* clim_tp_w_2nd_ns_d) +...
    (clim_dis_y_d/86400 .* clim_tp_y_ns_d) + (clim_dis_j_d/86400 .* clim_tp_j_2nd_ns_d) + (clim_dis_ye_d/86400 .* clim_tp_ye_ns_d) + ...
    (clim_dis_jw_d/86400 .* clim_tp_jw_ns_d);

total_dis_tp_3rd=(clim_dis_g_d/86400 .* clim_tp_g_ns_d) + (clim_dis_gc_d/86400 .*  clim_tp_gc_ns_d) + ...
    (clim_dis_gh_d/86400 .* clim_tp_gh_ns_d) + (clim_dis_gh2_d/86400 .* clim_tp_gh2_ns_d) + (clim_dis_d/86400 .* clim_tp_w_3rd_ns_d) +...
    (clim_dis_y_d/86400 .* clim_tp_y_ns_d) + (clim_dis_j_d/86400 .* clim_tp_j_3rd_ns_d) + (clim_dis_ye_d/86400 .* clim_tp_ye_ns_d) + ...
    (clim_dis_jw_d/86400 .* clim_tp_jw_ns_d);

figure;

grid on; xlim([1 365]); ylim([0 55]);

%sewer
figure;
plot(total_dis_tn'.*(0.3) ,'r'); hold on; plot(total_dis_tn_2nd.*(0.3) ,'g'); plot(total_dis_tn_3rd.*(0.3) ,'b'); 
grid on; xlim([1 365]); ylim([0 90]); ylabel('load NH4-N (TN*0.3) (g/s)'); xlabel('days'); 
text(10,50,num2str(mean(sj_tn_load_1.*(0.3)),'%0.2f'),'color','r');
text(10,45,num2str(mean(sj_tn_load_2.*(0.3)),'%0.2f'),'color','g');
text(10,40,num2str(mean(sj_tn_load_3.*(0.3)),'%0.2f'),'color','b');

figure;
plot(total_dis_tn'.*(0.6) ,'r'); hold on; plot(total_dis_tn_2nd.*(0.6) ,'g'); plot(total_dis_tn_3rd.*(0.6) ,'b'); 
grid on; xlim([1 365]); ylim([0 1000]); ylabel('load NO3-N (TN*0.6) (g/s)'); xlabel('days'); 
text(10,500,num2str(mean(sj_tn_load_1.*(0.6)),'%0.2f'),'color','r');
text(10,450,num2str(mean(sj_tn_load_2.*(0.6)),'%0.2f'),'color','g');
text(10,400,num2str(mean(sj_tn_load_3.*(0.6)),'%0.2f'),'color','b');

figure;
plot(total_dis_tn'.*(140/1000000) ,'r'); hold on; plot(total_dis_tn_2nd.*(140/1000000) ,'g'); plot(total_dis_tn_3rd.*(140/1000000) ,'b'); 
grid on; xlim([1 365]); ylim([0 90]); ylabel('load NH4-N (TN *140ppm) (g/s)'); xlabel('days'); 
text(10,50,num2str(mean(sj_tn_load_1.*(140/1000000)),'%0.2f'),'color','r');
text(10,45,num2str(mean(sj_tn_load_2.*(140/1000000)),'%0.2f'),'color','g');
text(10,40,num2str(mean(sj_tn_load_3.*(140/1000000)),'%0.2f'),'color','b');

figure;
plot(total_dis_tn'.*(136/1000000) ,'r'); hold on; plot(total_dis_tn_2nd.*(136/1000000) ,'g'); plot(total_dis_tn_3rd.*(136/1000000) ,'b'); 
grid on; xlim([1 365]); ylim([0 1000]); ylabel('load NO3-N (TN *136ppm) (g/s)'); xlabel('days'); 
text(10,500,num2str(mean(sj_tn_load_1.*(136/1000000)),'%0.2f'),'color','r');
text(10,450,num2str(mean(sj_tn_load_2.*(136/1000000)),'%0.2f'),'color','g');
text(10,400,num2str(mean(sj_tn_load_3.*(136/1000000)),'%0.2f'),'color','b');


%sewer vs sj       
figure; hold on;
plot(sj_tn_load_1,'r'); plot(sj_tn_load_2,'g'); 
plot(sj_tn_load_3,'b'); grid on;
plot(total_dis_tn,'r--'); hold on; plot(total_dis_tn_2nd,'g--'); plot(total_dis_tn_3rd,'b--'); 
text(10,1500,num2str(mean(sj_tn_load_1),'%0.2f'),'color','r');
text(10,1300,num2str(mean(sj_tn_load_2),'%0.2f'),'color','g');
text(10,1100,num2str(mean(sj_tn_load_3),'%0.2f'),'color','b');
text(50,1500,num2str(mean(total_dis_tn),'%0.2f'),'color','r');
text(50,1300,num2str(mean(total_dis_tn_2nd),'%0.2f'),'color','g');
text(50,1100,num2str(mean(total_dis_tn_3rd),'%0.2f'),'color','b');
ylabel('load TN (g/s)'); xlabel('days'); xlim([1 365])


figure; hold on;
plot(sj_tp_load_1,'r'); plot(sj_tp_load_2,'g'); 
plot(sj_tp_load_3,'b'); grid on;
plot(total_dis_tp,'r--'); hold on; plot(total_dis_tp_2nd,'g--'); plot(total_dis_tp_3rd,'b--'); 
text(10,35,num2str(mean(sj_tp_load_1),'%0.2f'),'color','r');
text(10,33,num2str(mean(sj_tp_load_2),'%0.2f'),'color','g');
text(10,31,num2str(mean(sj_tp_load_3),'%0.2f'),'color','b');
text(50,35,num2str(mean(total_dis_tp),'%0.2f'),'color','r');
text(50,33,num2str(mean(total_dis_tp_2nd),'%0.2f'),'color','g');
text(50,31,num2str(mean(total_dis_tp_3rd),'%0.2f'),'color','b');
ylabel('load TP (g/s)'); xlabel('days'); xlim([1 365])


