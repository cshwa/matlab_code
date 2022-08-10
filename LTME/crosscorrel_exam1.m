% Cross correlation  활용하여 인천 조위관측소 Residual tide와 서수도 기상관측
% 자료와의 상관관계를 분석
% crosscorrel_exam1.m
% input: weather_tide.dat 
%
clear all;close all;clc;
[yy,mm,dd,hr,min,w_deg,w_speed,air_pressure,wu_comp,wv_comp,rd_tide]...,
    = textread('weather_tide.dat ','','headerlines',1);
n=length(yy); 
% 전체자료를 적용해서 구하기, Lags는 'Series길이보다 작아야 한다.(length - 1)
% Bounds=2는 95 percent confidence interval

% 풍속과 Residual tide
subplot(4,1,1)
[XCF,Lags,Bounds] = crosscorr(w_speed,rd_tide,n-1,2); 
crosscorr(w_speed,rd_tide,n-1,2)
title('Wind speed VS Residual tide');ylim([-1, 1]);xlim([-110, 110]);

% 풍속 U성분과 Residual tide
subplot(4,1,2)
[XCF,Lags,Bounds] = crosscorr(wu_comp,rd_tide,n-1,2); 
crosscorr(wu_comp,rd_tide,n-1,2)
title('Wind U comp VS Residual tide');ylim([-1, 1]);xlim([-110, 110]);

% 풍속 V성분과 Residual tide
subplot(4,1,3)
[XCF,Lags,Bounds] = crosscorr(wv_comp,rd_tide,n-1,2); 
crosscorr(wv_comp,rd_tide,n-1,2)
title('Wind V comp VS Residual tide');ylim([-1, 1]);xlim([-110, 110]);

% 기압과 Residual tide 
subplot(4,1,4)
[XCF,Lags,Bounds] = crosscorr(air_pressure,rd_tide,n-1,2); 
crosscorr(air_pressure,rd_tide,n-1,2)
title('Air Pressure VS Residual tide');ylim([-1, 1]);xlim([-110, 110]);