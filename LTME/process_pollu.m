clc;clear all;close all

raw=xlsread('���縸 �ڷ� �ۺ�_ed_ed.xlsx','������');
raw2=xlsread('���縸 �ڷ� �ۺ�_ed_ed.xlsx','������ �ϵ� ����, ����');
 % 1-year 3-month 4-discharge 5-BOD 6-COD 7-SS 8-TN 9-TP
for i=1:11
    data(i,:,:)=raw(1+(i-1)*108:108*i,:);  % year, data , 1~9 
end

for i=1:11
    data_sum(i,:,:)=raw2(1+(i-1)*12:12*i,:);
end

meandata=squeeze(nanmean(data));
meandatas=squeeze(nanmean(data_sum));

% 1.�����������������ó����
% 2.�����߾��ϼ�����ó����
% 3.�����ϼ�����ó����
% 4.�����ϼ�����ó����
% 5.����������������������ó����
% 6.�������̻�������������ó����
% 7.���������������������ó����
% 8.�����ϼ�����ó����
% 9.�����ϼ�����ó����

T=[0.5:1:365];
t=[-15:30:380];

for ii=1:9
    eval(['riv',num2str(ii),'_dis=meandata(1+(',num2str(ii),'-1)*12:12*',num2str(ii),',4);riv',num2str(ii),'_TN=meandata(1+(',num2str(ii),'-1)*12:12*',num2str(ii),',8);riv',num2str(ii),'_TP=meandata(1+(',num2str(ii),'-1)*12:12*',num2str(ii),',9);']);
    eval(['riv',num2str(ii),'_dis=[riv',num2str(ii),'_dis(12);riv',num2str(ii),'_dis;riv',num2str(ii),'_dis(1)];']);
    eval(['riv',num2str(ii),'_TN=[riv',num2str(ii),'_TN(12);riv',num2str(ii),'_TN;riv',num2str(ii),'_TN(1)];']);
    eval(['riv',num2str(ii),'_TP=[riv',num2str(ii),'_TP(12);riv',num2str(ii),'_TP;riv',num2str(ii),'_TP(1)];']);
    eval(['riv',num2str(ii),'_dis_f=interp1(t,','riv',num2str(ii),'_dis,T);']);
    eval(['riv',num2str(ii),'_TN_f=interp1(t,','riv',num2str(ii),'_TN,T);']);
    eval(['riv',num2str(ii),'_TP_f=interp1(t,','riv',num2str(ii),'_TP,T);'])
end

riv_s_dis=meandatas(1+(1-1)*12:12*1,3);riv_s_TN=meandatas(1+(1-1)*12:12*1,7);
riv_s_dis=[riv_s_dis(12);riv_s_dis;riv_s_dis(1)];
riv_s_TN=[riv_s_TN(12);riv_s_TN;riv_s_TN(1)];
riv_s_dis_f=interp1(t,riv_s_dis,T);riv_s_TN_f=interp1(t,riv_s_TN,T);

save gy_pollution.mat riv_s_dis_f riv_s_TN_f riv1_dis_f riv2_dis_f riv3_dis_f riv4_dis_f riv5_dis_f riv6_dis_f riv7_dis_f riv8_dis_f riv9_dis_f  ...
    riv1_TN_f riv2_TN_f riv3_TN_f riv4_TN_f riv5_TN_f riv6_TN_f riv7_TN_f riv8_TN_f riv9_TN_f riv1_TP_f riv2_TP_f riv3_TP_f riv4_TP_f riv5_TP_f riv6_TP_f riv7_TP_f riv8_TP_f riv9_TP_f