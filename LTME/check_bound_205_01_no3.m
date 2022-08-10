close all
clear; 
clc;

in_2003=load('check_boundary_no3_input_2003_v2.mat');
in_2001=load('check_boundary_no3_input_2001.mat');
figure; hold on;
plot(in_2001.y_s_no3_205,'b');
plot(in_2003.y_s_no3_205,'r');

figure; hold on;
plot(in_2001.y_b_no3_205,'b');
plot(in_2003.y_b_no3_205,'r');





close all
clear; 
clc;

in_2003=load('check_boundary_no3_input_2003_v2.mat');
in_2001=load('check_boundary_no3_input_2001_v2.mat');
figure; hold on;
plot(in_2001.y_s_no3_205,'b');
plot(in_2003.y_s_no3_205,'r');

figure; hold on;
plot(in_2001.y_b_no3_205,'b');
plot(in_2003.y_b_no3_205,'r');


close all
clear; 
clc;

in_2003=load('no3_koem_input_400_16.mat');
in_2001=load('po4_koem_input_400_16.mat');
figure; hold on;
plot(in_2001.yp_w_no3_2,'b');
plot(in_2003.yp_w_no3_2,'r');

figure; hold on;
plot(in_2001.yp_w_no3_b_2,'b');
plot(in_2003.yp_w_no3_b_2,'r');