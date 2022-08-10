
% check each regime's input 
close all;clear; clc;

%% songjung
v5=load('sumjin(songjung)_polynomial_climate_to2004_advanced(v4)_3sig.mat','yp_w_*'); % advanced regime fitting for 1st regime
cd D:\장기생태\Dynamic\06_river\환경과학원
default=load('sumjin(songjung)_polynomial_climate_to2004_3sig.mat','yp_w_*'); % ~2004 : 1st regime, 2004~2018 : 2nd regime

figure; hold on
plot(v5.yp_w_no3_04'*1000/14,'r');
plot(default.yp_w_no3_04'*1000/14,'b+');
alpha(0.3)

% plot(v5.yp_w_no3_af'*1000/14,'r');
plot(default.yp_w_no3_af'*1000/14,'b+');  %2nd reigme
alpha(0.3)



figure; hold on
plot(v5.yp_w_nh4_04'*1000/14,'r');
plot(default.yp_w_nh4_04'*1000/14,'b+');
alpha(0.3)

% plot(v5.yp_w_nh4_af'*1000/14,'r');
plot(default.yp_w_nh4_af'*1000/14,'b+');  %2nd reigme
alpha(0.3)


po4_case = load('D:\장기생태\Dynamic\06_river\환경과학원\sumjin(songjung)_polynomial_climate_to2011_3sig_po4.mat');
% ~2010 : 1st regime(yp_w_po4_04), 2011~2018 : 2nd(yp_w_po4_af)

figure; hold on
plot(po4_case.yp_w_po4_04'.*1000./30.973762,'r');
plot(po4_case.yp_w_po4_af'.*1000./30.973762,'b+');
alpha(0.3)

close all; clear; clc;
%% Namgang
v5=load('Namgang_polynomial_climate_to2004_advanced(v4)_3sig.mat','yp_w_*'); % advanced regime fitting for 1st regime
cd D:\장기생태\Dynamic\06_river\환경과학원
default=load('Namgang_polynomial_climate_to2004_3sig.mat','yp_w_*'); % ~2004 : 1st regime, 2004~2018 : 2nd regime

figure; hold on
plot(v5.yp_w_no3_04'*1000/14,'r');
plot(default.yp_w_no3_04'*1000/14,'b+');
alpha(0.3)

% plot(v5.yp_w_no3_af'*1000/14,'r');
plot(default.yp_w_no3_af'*1000/14,'b+');  %2nd reigme
alpha(0.3)



figure; hold on
plot(v5.yp_w_nh4_04'*1000/14,'r');
plot(default.yp_w_nh4_04'*1000/14,'b+');
alpha(0.3)

% plot(v5.yp_w_nh4_af'*1000/14,'r');
plot(default.yp_w_nh4_af'*1000/14,'b+');  %2nd reigme
alpha(0.3)


po4_case = load('D:\장기생태\Dynamic\06_river\환경과학원\Namgang_polynomial_climate_to2011_3sig_po4.mat');
% ~2010 : 1st regime(yp_w_po4_04), 2011~2018 : 2nd(yp_w_po4_af)

figure; hold on
plot(po4_case.yp_w_po4_04'.*1000./30.973762,'r');
plot(po4_case.yp_w_po4_af'.*1000./30.973762,'b+');
alpha(0.3)


