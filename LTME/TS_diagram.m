% TS?diagram , DO?S and DO?T Diagram
% Input file: E03.dat
% m-file: TS_diagram1.m
% 
clc;clear all;close all
figure();      % T-S diagram �ۼ�
% input_file=input(' Input file  : ','s');   % �̹���� ����ϸ� ���� data�� �ѹ��� �ۼ�
% output_file=input(' Outout file : ','s');
% data_raw=load(input_file);

%---- 2015. ���� ����  -- ������ �Ǵ� �� 

% data_raw=load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_summer_ys\CTD\150827_neap\am_ebb\am_ebb_all.dat');    % �ϴ� �� �ڷḸ�� ���� �����ϰ� �Ѵ�.
% T=data_raw(:,4);
% S=data_raw(:,6);
data_raw=load('D:\mepl\data\01_Sumjin_yeosu_camp\2015_summer_ys\CTD\150827_neap\pm_flood\pm_flood_all.dat');
T=data_raw(:,4);
S=data_raw(:,6);
% T=[T;data_raw(:,4)];
% S=[S;data_raw(:,6)];
P=[-10 1000];
dens=0:1:30;            % �е� Grid����, ���� �� ����
xy=[0 35 0 35];           % ���� ���� & ���¹���
% tsdiagram(S,T,P,dens,'.b',xy);
tsdiagram(S,T,P,dens,'gx',xy,'markersize',20);     % tidiagram�� �ۼ��ϱ� ���� ���� (blue ?dot)
% print -dpsc TSdiagram1.ps     % TSdiagram1.ps�� ����

data_raw=load('D:\mepl\data\01_Sumjin_yeosu_camp\2015_summer_ys\CTD\150816_spring\am_flood\am_flood_all.dat');    % �ϴ� �� �ڷḸ�� ���� �����ϰ� �Ѵ�.
T=data_raw(:,4);
S=data_raw(:,6);
plot(S,T,'g.','markersize',20)

% data_raw=load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_summer_ys\CTD\150816_spring\pm_ebb\pm_ebb_all.dat');
% T=[T;data_raw(:,4)];
% S=[S;data_raw(:,6)];

%--- 2013. ����
% data_raw=load('D:\SilverStar\data\01_SumJin\2013_08_13_ys\High\high_all.dat');    
% T=data_raw(:,4);
% S=data_raw(:,6);
% plot(S,T,'g.')
% data_raw=load('D:\SilverStar\data\01_SumJin\2013_08_13_ys\Low\low_all.dat');
% T=[T;data_raw(:,4)];
% S=[S;data_raw(:,6)];


%--- 2012. ����
% data_raw=load('D:\SilverStar\data\01_SumJin\2012_08_07_ys\high_total_data.txt');    
% T=data_raw(:,3);
% S=data_raw(:,5);
% plot(S,T,'k.')
% data_raw=load('D:\SilverStar\data\01_SumJin\2012_08_07_ys\low_total_data.txt');
% T=[T;data_raw(:,3)];
% S=[S;data_raw(:,5)];




%--- 2015. �ܿ� ����
data_raw=load('D:\mepl\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150226_ys_neap tide\am_flood\am_flood_all.dat');    % �ϴ� �� �ڷḸ�� ���� �����ϰ� �Ѵ�.
T=data_raw(:,4);
S=data_raw(:,6);
plot(S,T,'bx','markersize',15)
% data_raw=load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150226_ys_neap tide\pm_ebb\pm_ebb_all.dat');
% T=[T;data_raw(:,4)];
% S=[S;data_raw(:,6)];



% fclose('all');
data_raw=load('D:\mepl\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150305 ys_spring tide\am_flood\am_flood_all.dat');    % �ϴ� �� �ڷḸ�� ���� �����ϰ� �Ѵ�.
T=data_raw(:,4);
S=data_raw(:,6);
plot(S,T,'b.')
% data_raw=load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150305 ys_spring tide\pm_ebb\pm_ebb_all.dat');
% T=[T;data_raw(:,4)];
% S=[S;data_raw(:,6)];


%% �� ���� ���� ��

%---- 2015. �� ����

% data_raw=load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_spring_ys\CTD\150513_CTD_neap\am_ebb\am_ebb_all.dat');    % �ϴ� �� �ڷḸ�� ���� �����ϰ� �Ѵ�.
% T=data_raw(:,4);
% S=data_raw(:,6);

data_raw=load('D:\meplr\data\01_Sumjin_yeosu_camp\2015_spring_ys\CTD\150513_CTD_neap\pm_flood\pm_flood_all.dat');
T=data_raw(:,4);
S=data_raw(:,6);
plot(S,T,'rx','markersize',15)
data_raw=load('D:\mepl\data\01_Sumjin_yeosu_camp\2015_spring_ys\CTD\150519_CTD_spring\am_flood\am_flood_all.dat');    % �ϴ� �� �ڷḸ�� ���� �����ϰ� �Ѵ�.
T=data_raw(:,4);
S=data_raw(:,6);
plot(S,T,'r.','markersize',20)
%{
%---- ������ �Ǵ� ���� �����ϰ� �־�� �� �κ� -------------------
P=[-10 1000];
dens=0:1:30;            % �е� Grid����, ���� �� ����
xy=[0 35 0 35];           % ���� ���� & ���¹���
tsdiagram(S,T,P,dens,'k+',xy);     % tidiagram�� �ۼ��ϱ� ���� ���� (blue ?dot)
print -dpsc TSdiagram1.ps     % TSdiagram1.ps�� ����
%---------------------------------------------------------------
% data_raw=load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_spring_ys\CTD\150519_CTD_spring\pm_ebb\pm_ebb_all.dat');
% T=[T;data_raw(:,4)];
% S=[S;data_raw(:,6)];
% plot(S,T,'k.')

%--- 2012�� ��-------------------------------------------------------------
data_raw=load('E:\11.���\������_1�ܰ�\mtkwak\������_gwangyang\[01] ProcessedData\sort_data\ts_diagram_all\20120620_high_ts_dia.dat');    % �ϴ� �� �ڷḸ�� ���� �����ϰ� �Ѵ�.
T=data_raw(:,2);
S=data_raw(:,3);
plot(S,T,'g*')

% data_raw=load('E:\11.���\������_1�ܰ�\mtkwak\������_gwangyang\[01] ProcessedData\sort_data\ts_diagram_all\20120620_low_ts_dia.dat');
% T=[T;data_raw(:,2)];
% S=[S;data_raw(:,3)];
% plot(S,T,'y.')

%--- 2013�� ��-------------------------------------------------------------
data_raw=load('E:\11.���\������_1�ܰ�\mtkwak\������_gwangyang\[01] ProcessedData\sort_data\ts_diagram_all\20130603_high_ts_dia.dat');    % �ϴ� �� �ڷḸ�� ���� �����ϰ� �Ѵ�.
T=data_raw(:,2);
S=data_raw(:,3);
plot(S,T,'r.')

% data_raw=load('E:\11.���\������_1�ܰ�\mtkwak\������_gwangyang\[01] ProcessedData\sort_data\ts_diagram_all\20130603_low_ts_dia.dat');
% T=[T;data_raw(:,2)];
% S=[S;data_raw(:,3)];
% plot(S,T,'b.')

%}
fclose('all');