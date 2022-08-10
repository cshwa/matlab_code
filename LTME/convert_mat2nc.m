close all; clear all; clc; 
load('2009_v10.mat');
create_roms_forcing_V(fname,inp_V_daily,time,yy,yy)

close all; clear all; clc; 
load('2009_u10.mat');
create_roms_forcing_U(fname,inp_U_daily,time,yy,yy)

close all; clear all; clc; 
load('2009_Tair.mat');
create_roms_forcing_T(fname,inp_T_daily,time,yy,yy)

close all; clear all; clc; 
load('2009_msl.mat');
create_roms_forcing_P(fname,inp_P_daily,time,yy,yy)

close all; clear all; clc; 
load('2009_dewt.mat');
create_roms_forcing_D(fname,inp_D_daily,time,yy,yy)

close all; clear all; clc; 
load('2009_swradd.mat');
create_roms_forcing_SWD(fname,inp_SW_daily,time,yy,yy)