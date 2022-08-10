% m-file: fft_exam3.m
% 
clc;clear all;close all;
[YY,MM,DD,HR,TH]= textread('jeju0510.dat','');
xx=TH-mean(TH);
Fs = 1000;   
periodogram(xx,[],'onesided',1024,Fs)