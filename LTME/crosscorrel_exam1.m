% Cross correlation  Ȱ���Ͽ� ��õ ���������� Residual tide�� ������ ������
% �ڷ���� ������踦 �м�
% crosscorrel_exam1.m
% input: weather_tide.dat 
%
clear all;close all;clc;
[yy,mm,dd,hr,min,w_deg,w_speed,air_pressure,wu_comp,wv_comp,rd_tide]...,
    = textread('weather_tide.dat ','','headerlines',1);
n=length(yy); 
% ��ü�ڷḦ �����ؼ� ���ϱ�, Lags�� 'Series���̺��� �۾ƾ� �Ѵ�.(length - 1)
% Bounds=2�� 95 percent confidence interval

% ǳ�Ӱ� Residual tide
subplot(4,1,1)
[XCF,Lags,Bounds] = crosscorr(w_speed,rd_tide,n-1,2); 
crosscorr(w_speed,rd_tide,n-1,2)
title('Wind speed VS Residual tide');ylim([-1, 1]);xlim([-110, 110]);

% ǳ�� U���а� Residual tide
subplot(4,1,2)
[XCF,Lags,Bounds] = crosscorr(wu_comp,rd_tide,n-1,2); 
crosscorr(wu_comp,rd_tide,n-1,2)
title('Wind U comp VS Residual tide');ylim([-1, 1]);xlim([-110, 110]);

% ǳ�� V���а� Residual tide
subplot(4,1,3)
[XCF,Lags,Bounds] = crosscorr(wv_comp,rd_tide,n-1,2); 
crosscorr(wv_comp,rd_tide,n-1,2)
title('Wind V comp VS Residual tide');ylim([-1, 1]);xlim([-110, 110]);

% ��а� Residual tide 
subplot(4,1,4)
[XCF,Lags,Bounds] = crosscorr(air_pressure,rd_tide,n-1,2); 
crosscorr(air_pressure,rd_tide,n-1,2)
title('Air Pressure VS Residual tide');ylim([-1, 1]);xlim([-110, 110]);