% TS?diagram , DO?S and DO?T Diagram
% Input file: E03.dat
% m-file: TS_diagram1.m
% 
clc;clear all;close all
figure();      % T-S diagram 작성
% input_file=input(' Input file  : ','s');   % 이방법을 사용하면 여러 data를 한번에 작성
% output_file=input(' Outout file : ','s');
% data_raw=load(input_file);

%---- 2015. 여름 관측  -- 기준이 되는 것 

% data_raw=load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_summer_ys\CTD\150827_neap\am_ebb\am_ebb_all.dat');    % 일단 한 자료만을 직접 선택하게 한다.
% T=data_raw(:,4);
% S=data_raw(:,6);
data_raw=load('D:\mepl\data\01_Sumjin_yeosu_camp\2015_summer_ys\CTD\150827_neap\pm_flood\pm_flood_all.dat');
T=data_raw(:,4);
S=data_raw(:,6);
% T=[T;data_raw(:,4)];
% S=[S;data_raw(:,6)];
P=[-10 1000];
dens=0:1:30;            % 밀도 Grid설정, 범위 및 간격
xy=[0 35 0 35];           % 염분 범위 & 수온범위
% tsdiagram(S,T,P,dens,'.b',xy);
tsdiagram(S,T,P,dens,'gx',xy,'markersize',20);     % tidiagram을 작성하기 위한 내용 (blue ?dot)
% print -dpsc TSdiagram1.ps     % TSdiagram1.ps로 저장

data_raw=load('D:\mepl\data\01_Sumjin_yeosu_camp\2015_summer_ys\CTD\150816_spring\am_flood\am_flood_all.dat');    % 일단 한 자료만을 직접 선택하게 한다.
T=data_raw(:,4);
S=data_raw(:,6);
plot(S,T,'g.','markersize',20)

% data_raw=load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_summer_ys\CTD\150816_spring\pm_ebb\pm_ebb_all.dat');
% T=[T;data_raw(:,4)];
% S=[S;data_raw(:,6)];

%--- 2013. 여름
% data_raw=load('D:\SilverStar\data\01_SumJin\2013_08_13_ys\High\high_all.dat');    
% T=data_raw(:,4);
% S=data_raw(:,6);
% plot(S,T,'g.')
% data_raw=load('D:\SilverStar\data\01_SumJin\2013_08_13_ys\Low\low_all.dat');
% T=[T;data_raw(:,4)];
% S=[S;data_raw(:,6)];


%--- 2012. 여름
% data_raw=load('D:\SilverStar\data\01_SumJin\2012_08_07_ys\high_total_data.txt');    
% T=data_raw(:,3);
% S=data_raw(:,5);
% plot(S,T,'k.')
% data_raw=load('D:\SilverStar\data\01_SumJin\2012_08_07_ys\low_total_data.txt');
% T=[T;data_raw(:,3)];
% S=[S;data_raw(:,5)];




%--- 2015. 겨울 관측
data_raw=load('D:\mepl\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150226_ys_neap tide\am_flood\am_flood_all.dat');    % 일단 한 자료만을 직접 선택하게 한다.
T=data_raw(:,4);
S=data_raw(:,6);
plot(S,T,'bx','markersize',15)
% data_raw=load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150226_ys_neap tide\pm_ebb\pm_ebb_all.dat');
% T=[T;data_raw(:,4)];
% S=[S;data_raw(:,6)];



% fclose('all');
data_raw=load('D:\mepl\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150305 ys_spring tide\am_flood\am_flood_all.dat');    % 일단 한 자료만을 직접 선택하게 한다.
T=data_raw(:,4);
S=data_raw(:,6);
plot(S,T,'b.')
% data_raw=load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150305 ys_spring tide\pm_ebb\pm_ebb_all.dat');
% T=[T;data_raw(:,4)];
% S=[S;data_raw(:,6)];


%% 봄 관측 끼지 비교

%---- 2015. 봄 관측

% data_raw=load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_spring_ys\CTD\150513_CTD_neap\am_ebb\am_ebb_all.dat');    % 일단 한 자료만을 직접 선택하게 한다.
% T=data_raw(:,4);
% S=data_raw(:,6);

data_raw=load('D:\meplr\data\01_Sumjin_yeosu_camp\2015_spring_ys\CTD\150513_CTD_neap\pm_flood\pm_flood_all.dat');
T=data_raw(:,4);
S=data_raw(:,6);
plot(S,T,'rx','markersize',15)
data_raw=load('D:\mepl\data\01_Sumjin_yeosu_camp\2015_spring_ys\CTD\150519_CTD_spring\am_flood\am_flood_all.dat');    % 일단 한 자료만을 직접 선택하게 한다.
T=data_raw(:,4);
S=data_raw(:,6);
plot(S,T,'r.','markersize',20)
%{
%---- 기준이 되는 것이 포함하고 있어야 할 부분 -------------------
P=[-10 1000];
dens=0:1:30;            % 밀도 Grid설정, 범위 및 간격
xy=[0 35 0 35];           % 염분 범위 & 수온범위
tsdiagram(S,T,P,dens,'k+',xy);     % tidiagram을 작성하기 위한 내용 (blue ?dot)
print -dpsc TSdiagram1.ps     % TSdiagram1.ps로 저장
%---------------------------------------------------------------
% data_raw=load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_spring_ys\CTD\150519_CTD_spring\pm_ebb\pm_ebb_all.dat');
% T=[T;data_raw(:,4)];
% S=[S;data_raw(:,6)];
% plot(S,T,'k.')

%--- 2012년 봄-------------------------------------------------------------
data_raw=load('E:\11.사업\장기생태_1단계\mtkwak\장기생태_gwangyang\[01] ProcessedData\sort_data\ts_diagram_all\20120620_high_ts_dia.dat');    % 일단 한 자료만을 직접 선택하게 한다.
T=data_raw(:,2);
S=data_raw(:,3);
plot(S,T,'g*')

% data_raw=load('E:\11.사업\장기생태_1단계\mtkwak\장기생태_gwangyang\[01] ProcessedData\sort_data\ts_diagram_all\20120620_low_ts_dia.dat');
% T=[T;data_raw(:,2)];
% S=[S;data_raw(:,3)];
% plot(S,T,'y.')

%--- 2013년 봄-------------------------------------------------------------
data_raw=load('E:\11.사업\장기생태_1단계\mtkwak\장기생태_gwangyang\[01] ProcessedData\sort_data\ts_diagram_all\20130603_high_ts_dia.dat');    % 일단 한 자료만을 직접 선택하게 한다.
T=data_raw(:,2);
S=data_raw(:,3);
plot(S,T,'r.')

% data_raw=load('E:\11.사업\장기생태_1단계\mtkwak\장기생태_gwangyang\[01] ProcessedData\sort_data\ts_diagram_all\20130603_low_ts_dia.dat');
% T=[T;data_raw(:,2)];
% S=[S;data_raw(:,3)];
% plot(S,T,'b.')

%}
fclose('all');