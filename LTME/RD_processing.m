clc; clear all; close all;

% mat ���Ϸ� ������� ���� �б�

st = 2000; en = 2015;
% in_path =[num2str(yid),'\hourly_mean\'];
% out_path=[num2str(yid),'\out\sst\']; % ����� ������
% filename = [in_path,mid_name,midd,'.mat'];
t_RD = [];

for i = st:en
    mid=[num2char(i,4)];
    temp = load(['RD',mid,'.mat']);    
    data = temp.temp2;
    t_RD = [t_RD data];
end

% save t_RD t_RD
%% ���ذ��� �˾ƺ���
mean(nanmean(t_RD))
mean(nanmedian(t_RD))
mean(mode(t_RD))

%%
% [a b]= size(t_RD);
% figure()
