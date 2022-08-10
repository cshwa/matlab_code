% 이 코드는 ARGOdata의 Pressure, Temp., Salinity data를 이용하여 밀도를 계산
% by peter
% M-file: argo_cal_density.m 
% seawater toolbox(ver2_0_1)필요
% 이 toolbox는 'ftp://ftp.marine.csiro.au/pub/morgan/seawater/' 다운로드 가능
% 
clc; clear all 
[header,data] = hdrload('E:\MATLAB7\work\ARGO_kma\2900300_096.txt'); 
P = data(:,1);       % Pressure 
T = data(:,2);       % Temperature 
S = data(:,3);       % Salinity 
Density = sw_dens(S,T,P);  % calculate the density 

% save as Prssure,Temperature,Salinity,Density 
Ndata = [P,T,S,Density]; 
fid = fopen('2900300_096_den.dat','w'); 
fprintf(fid,'%6.3f%8.3f%10.4f%11.4f\n',Ndata') 
fclose(fid)