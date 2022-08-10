%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make boundary file
% % %  1 micro gram to mili gram  = 0.001
%%% mg/L to mM/m^3;
% sumjin_do = yp_w_do'*0.7*44.661;
% sumjin_chl = yp_w_chl';
% sumjin_no3 = yp_w_no3'*1000/14;
% sumjin_nh4 = yp_w_nh4'*1000/14;
close all; clear all; clc;

for t_i = 1:3  % regime

clearvars -except t_i

%% KOEM
%  from polynomial_fitting_on_ecological_variable_v3_3sig(boundary)_chl_2019.m
% there are 2 regimes for nh4 (KOEM 7st.)
koem=load('koem_input_fix_to06to15_3sig(3regime).mat'); 
if t_i == 1
    y_s_nh4 = koem.mon_clim_sur_nh4_1_in ./ koem.N_MW; %%% ug/L to mM/m^3;
    y_b_nh4 = koem.mon_clim_bot_nh4_1_in ./ koem.N_MW; %%% ug/L to mM/m^3;
elseif t_i == 2
    y_s_nh4 = koem.mon_clim_sur_nh4_2_in ./ koem.N_MW; %%% ug/L to mM/m^3;
    y_b_nh4 = koem.mon_clim_bot_nh4_2_in ./ koem.N_MW; %%% ug/L to mM/m^3;
elseif t_i == 3
    y_s_nh4 = koem.mon_clim_sur_nh4_3_in ./ koem.N_MW; %%% ug/L to mM/m^3;
    y_b_nh4 = koem.mon_clim_bot_nh4_3_in ./ koem.N_MW; %%% ug/L to mM/m^3;
end
 
if t_i == 1
    y_s_chl = koem.mon_clim_sur_chl_1_in; %%% ug/L;
    y_b_chl = koem.mon_clim_bot_chl_1_in; %%% ug/L;
elseif t_i == 2
    y_s_chl = koem.mon_clim_sur_chl_2_in; %%% ug/L;
    y_b_chl = koem.mon_clim_bot_chl_2_in; %%% ug/L;
elseif t_i == 3
    y_s_chl = koem.mon_clim_sur_chl_3_in; %%% ug/L;
    y_b_chl = koem.mon_clim_bot_chl_3_in; %%% ug/L;
end

%% it's kodc
%% from plot_KODC_6_20501_polynomial_to06to15.m
%  there are 3 regimes for no3 (KODC)
%  1997~2006, 2007~2015, 2016~2019
kodc205=load('kodc_input_fix_to06to15_3sig(3regime)_20501.mat');

if t_i == 1
    y_s_no3_205 = kodc205.mon_clim_sur_no3_1_in./ koem.N_MW; %%% mg/L to mM/m^3;
    y_b_no3_205 = kodc205.mon_clim_bot_no3_1_in./ koem.N_MW; %%% mg/L to mM/m^3;
elseif t_i == 2
    y_s_no3_205 = kodc205.mon_clim_sur_no3_2_in./ koem.N_MW; %%% mg/L to mM/m^3;
    y_b_no3_205 = kodc205.mon_clim_bot_no3_2_in./ koem.N_MW; %%% mg/L to mM/m^3;
elseif t_i == 3
    y_s_no3_205 = kodc205.mon_clim_sur_no3_3_in./ koem.N_MW; %%% mg/L to mM/m^3;
    y_b_no3_205 = kodc205.mon_clim_bot_no3_3_in./ koem.N_MW; %%% mg/L to mM/m^3;
end

% DO is reconstructed
clearvars do_re_b
load do_b_kodc_input_205_01.mat % do_re_b (day)
if t_i == 1
    y_s_do_205 = kodc205.mon_clim_sur_do_1_in .*0.7.*44.661; %%% mg/L to mM/m^3;
    y_b_do_205 = kodc205.mon_clim_bot_do_1_in .*0.7.*44.661; %%% mg/L to mM/m^3;
elseif t_i == 2
    y_s_do_205 = kodc205.mon_clim_sur_do_2_in .*0.7.*44.661; %%% mg/L to mM/m^3;
    y_b_do_205 = kodc205.mon_clim_bot_do_2_in .*0.7.*44.661; %%% mg/L to mM/m^3;
elseif t_i == 3
    y_s_do_205 = kodc205.mon_clim_sur_do_3_in .*0.7.*44.661; %%% mg/L to mM/m^3;
    y_b_do_205 = kodc205.mon_clim_bot_do_3_in .*0.7.*44.661; %%% mg/L to mM/m^3;
end

% PO4
if t_i == 1
    y_s_po4_205 = kodc205.mon_clim_sur_po4_1_in; %%% mg/L to mM/m^3;
    y_b_po4_205 = kodc205.mon_clim_bot_po4_1_in; %%% mg/L to mM/m^3;
elseif t_i == 2
    y_s_po4_205 = kodc205.mon_clim_sur_po4_2_in; %%% mg/L to mM/m^3;
    y_b_po4_205 = kodc205.mon_clim_bot_po4_2_in; %%% mg/L to mM/m^3;
elseif t_i == 3
    y_s_po4_205 = kodc205.mon_clim_sur_po4_3_in; %%% mg/L to mM/m^3;
    y_b_po4_205 = kodc205.mon_clim_bot_po4_3_in; %%% mg/L to mM/m^3;
end

%%  plot_KODC_6_40016_polynomial_to06to15.m
%  there are 3 regimes for no3 (KODC)
%  1997~2006, 2007~2015, 2016~2019
kodc400=load('kodc_input_fix_to06to15_3sig(3regime)_40016.mat');

if t_i == 1
    y_s_no3_400 = kodc400.mon_clim_sur_no3_1_in./ koem.N_MW; %%% mg/L to mM/m^3;
    y_b_no3_400 = kodc400.mon_clim_bot_no3_1_in./ koem.N_MW; %%% mg/L to mM/m^3;
elseif t_i == 2
    y_s_no3_400 = kodc400.mon_clim_sur_no3_2_in./ koem.N_MW; %%% mg/L to mM/m^3;
    y_b_no3_400 = kodc400.mon_clim_bot_no3_2_in./ koem.N_MW; %%% mg/L to mM/m^3;
elseif t_i == 3
    y_s_no3_400 = kodc400.mon_clim_sur_no3_3_in./ koem.N_MW; %%% mg/L to mM/m^3;
    y_b_no3_400 = kodc400.mon_clim_bot_no3_3_in./ koem.N_MW; %%% mg/L to mM/m^3;
end

% DO
if t_i == 1
    y_s_do_400 = kodc400.mon_clim_sur_do_1_in .*0.7.*44.661; %%% mg/L to mM/m^3;
    y_b_do_400 = kodc400.mon_clim_bot_do_1_in .*0.7.*44.661; %%% mg/L to mM/m^3;
elseif t_i == 2
    y_s_do_400 = kodc400.mon_clim_sur_do_2_in .*0.7.*44.661; %%% mg/L to mM/m^3;
    y_b_do_400 = kodc400.mon_clim_bot_do_2_in .*0.7.*44.661; %%% mg/L to mM/m^3;
elseif t_i == 3
    y_s_do_400 = kodc400.mon_clim_sur_do_3_in .*0.7.*44.661; %%% mg/L to mM/m^3;
    y_b_do_400 = kodc400.mon_clim_bot_do_3_in .*0.7.*44.661; %%% mg/L to mM/m^3;
end

% PO4
if t_i == 1
    y_s_po4_400 = kodc400.mon_clim_sur_po4_1_in; %%% mg/L to mM/m^3;
    y_b_po4_400 = kodc400.mon_clim_bot_po4_1_in; %%% mg/L to mM/m^3;
elseif t_i == 2
    y_s_po4_400 = kodc400.mon_clim_sur_po4_2_in; %%% mg/L to mM/m^3;
    y_b_po4_400 = kodc400.mon_clim_bot_po4_2_in ; %%% mg/L to mM/m^3;
elseif t_i == 3
    y_s_po4_400 = kodc400.mon_clim_sur_po4_3_in; %%% mg/L to mM/m^3;
    y_b_po4_400 = kodc400.mon_clim_bot_po4_3_in; %%% mg/L to mM/m^3;
end

%%
y_s_no3_205_m=y_s_no3_205;
y_b_no3_205_m=y_b_no3_205;
y_s_no3_400_m=y_s_no3_400;
y_b_no3_400_m=y_b_no3_400;
y_s_po4_205_m=y_s_po4_205;
y_b_po4_205_m=y_b_po4_205;
y_s_po4_400_m=y_s_po4_400;
y_b_po4_400_m=y_b_po4_400;
do_re_b_205_m = y_b_do_205;
do_re_b_400_m = y_b_do_400;
do_re_s_205_m = y_s_do_205;
do_re_s_400_m = y_s_do_400;   


%% vertical interp.  %time s_rho lon_rho
obc_path ='D:\장기생태\Dynamic\07_boundary_ts\Gwangyang_bry_ykang\Gwangyang_add_uv_bry\';
% ncdump([obc_path, 'Gwangyang_new_Y1980.nc'])

lon_rho = ncread('grid_sumjin_v1970_fix_3m.nc','lon_rho');
lat_rho = ncread('grid_sumjin_v1970_fix_3m.nc','lat_rho');
h = ncread('grid_sumjin_v1970_fix_3m.nc','h');
lon_rho_1 = repmat(lon_rho(:,1),1,20); lon_rho_end = lon_rho(:,end);
lat_rho_1 = repmat(lat_rho(1,:),1,20); lat_rho_end = lon_rho(end,:);
z_1 = NaN(length(lon_rho(:,1)),20);
zz = [20:-1:1];
for i = 1:20
    z_1(:,i) = zz(i);
end

t_w=ncread([obc_path, 'Gwangyang_new_Y1980.nc'],'temp_west');
t_e=ncread([obc_path, 'Gwangyang_new_Y1980.nc'],'temp_east');
t_n=ncread([obc_path, 'Gwangyang_new_Y1980.nc'],'temp_north');
t_s=ncread([obc_path, 'Gwangyang_new_Y1980.nc'],'temp_south');

%% chl
chl_in = NaN(12,20); chl_in(1:12,20)=y_s_chl; % 2nd reigme
chl_in(1:12,1)=y_b_chl;

t = 1:20;
for i = 1:12
    clearvars temp_s
temp_s = chl_in(i,:);
temp_s(isnan(temp_s)) = interp1( t(~isnan(temp_s)), temp_s(~isnan(temp_s)), t(isnan(temp_s)) );  
chl_in(i,:) = temp_s;
end
chl_input_sn = ones(size(t_s,1),size(t_s,2),size(t_s,3));
chl_input_we = ones(size(t_w,1),size(t_w,2),size(t_w,3));

for i = 1:12
    for j = 1:20
        chl_input_sn(:,j,i) = chl_input_sn(:,j,i) .* chl_in(i,j);
        chl_input_we(:,j,i) = chl_input_we(:,j,i) .* chl_in(i,j);
    end
end

%% nh4
nh4_in = NaN(12,20); nh4_in(1:12,20)=y_s_nh4; nh4_in(1:12,1)=y_b_nh4;

for i = 1:12
    clearvars temp_s
temp_s = nh4_in(i,:);
temp_s(isnan(temp_s)) = interp1( t(~isnan(temp_s)), temp_s(~isnan(temp_s)), t(isnan(temp_s)) );  
nh4_in(i,:) = temp_s;
end
nh4_input_sn = ones(size(t_s,1),size(t_s,2),size(t_s,3));
nh4_input_we = ones(size(t_w,1),size(t_w,2),size(t_w,3));

for i = 1:12
    for j = 1:20
        nh4_input_sn(:,j,i) = nh4_input_sn(:,j,i) .* nh4_in(i,j);
        nh4_input_we(:,j,i) = nh4_input_we(:,j,i) .* nh4_in(i,j);
    end
end


% 20501 34.3717	127.8083
% 40016 34.5167	128.2033

%% do 
do_in_s = NaN(2,12); do_in_s(1,:) = do_re_s_205_m; do_in_s(2,:) = do_re_s_400_m;
do_in_b = NaN(2,12); do_in_b(1,:) = do_re_b_205_m; do_in_b(2,:) = do_re_b_400_m;
do_in = NaN(2,20,12); do_in(:,1,:) = do_in_b; do_in(:,20,:) = do_in_s;
x_in = NaN(2,20); x_in(1,:) = 127.8083; x_in(2,:) = 128.2033;
z_in = NaN(2,20); 
zz = [20:-1:1];
for i = 1:20
    z_in(:,i) = zz(i);
end

% vertical interp on KODC st.
for i = 1:2
        for k = 1:12
            clearvars temp_s 
            temp_s = squeeze(do_in(i,:,k));
            temp_z = zz;
                temp_s(isnan(temp_s)) = interp1( temp_z(~isnan(temp_s)), temp_s(~isnan(temp_s)), temp_z(isnan(temp_s))); 
            do_in(i,:,k) = temp_s;
        end
end


for i = 1:12
    clearvars tem tem_interp
    tem = squeeze(do_in(:,:,i));
    tem_interp = griddata(x_in,z_in,tem,lon_rho_1, z_1,'nearest');
    do_input(:,:,i) = tem_interp;
end


find(isnan(do_input)==1)

% do_input_sn = ones(size(t_s,1),size(t_s,2),size(t_s,3));
clearvars do_input_we
do_input_we = ones(size(t_w,1),size(t_w,2),size(t_w,3));

for i = 1:12
    for j = 1:20
        do_input_we(:,j,i) = squeeze(do_input_we(:,j,i)) .* squeeze(do_in(2,j,i));
    end
end

%% no3
no3_in_s = NaN(2,12); no3_in_s(1,:) = y_s_no3_205_m; no3_in_s(2,:) = y_s_no3_400_m;
no3_in_b = NaN(2,12); no3_in_b(1,:) = y_b_no3_205_m; no3_in_b(2,:) = y_b_no3_400_m;
no3_in = NaN(2,20,12); no3_in(:,1,:) = no3_in_b; no3_in(:,20,:) = no3_in_s;

% vertical interp on KODC st.
for i = 1:2
        for k = 1:12
            clearvars temp_s 
            temp_s = squeeze(no3_in(i,:,k));
            temp_z = zz;
            temp_s(isnan(temp_s)) = interp1( temp_z(~isnan(temp_s)), temp_s(~isnan(temp_s)), temp_z(isnan(temp_s))); 
            no3_in(i,:,k) = temp_s;
        end
end

for i = 1:12
    clearvars tem tem_interp
    tem = squeeze(no3_in(:,:,i));
    tem_interp = griddata(x_in,z_in,tem,lon_rho_1, z_1,'nearest');
    no3_input(:,:,i) = tem_interp;
end

find(isnan(no3_input)==1)

% no3_input_sn = ones(size(t_s,1),size(t_s,2),size(t_s,3));
clearvars no3_input_we
no3_input_we = ones(size(t_w,1),size(t_w,2),size(t_w,3));

for i = 1:12
    for j = 1:20
        no3_input_we(:,j,i) = squeeze(no3_input_we(:,j,i)) .* squeeze(no3_in(2,j,i));
    end
end

%% no3
no3_in_s = NaN(2,12); no3_in_s(1,:) = y_s_no3_205_m; no3_in_s(2,:) = y_s_no3_400_m;
no3_in_b = NaN(2,12); no3_in_b(1,:) = y_b_no3_205_m; no3_in_b(2,:) = y_b_no3_400_m;
no3_in = NaN(2,20,12); no3_in(:,1,:) = no3_in_b; no3_in(:,20,:) = no3_in_s;

% vertical interp on KODC st.
for i = 1:2
        for k = 1:12
            clearvars temp_s 
            temp_s = squeeze(no3_in(i,:,k));
            temp_z = zz;
            temp_s(isnan(temp_s)) = interp1( temp_z(~isnan(temp_s)), temp_s(~isnan(temp_s)), temp_z(isnan(temp_s))); 
            no3_in(i,:,k) = temp_s;
        end
end

for i = 1:12
    clearvars tem tem_interp
    tem = squeeze(no3_in(:,:,i));
    tem_interp = griddata(x_in,z_in,tem,lon_rho_1, z_1,'nearest');
    no3_input(:,:,i) = tem_interp;
end

find(isnan(no3_input)==1)

% no3_input_sn = ones(size(t_s,1),size(t_s,2),size(t_s,3));
clearvars no3_input_we
no3_input_we = ones(size(t_w,1),size(t_w,2),size(t_w,3));

for i = 1:12
    for j = 1:20
        no3_input_we(:,j,i) = squeeze(no3_input_we(:,j,i)) .* squeeze(no3_in(2,j,i));
    end
end

%% po4
po4_in_s = NaN(2,12); po4_in_s(1,:) = y_s_po4_205_m; po4_in_s(2,:) = y_s_po4_400_m;
po4_in_b = NaN(2,12); po4_in_b(1,:) = y_b_po4_205_m; po4_in_b(2,:) = y_b_po4_400_m;
po4_in = NaN(2,20,12); po4_in(:,1,:) = po4_in_b; po4_in(:,20,:) = po4_in_s;

% vertical interp on KODC st.
for i = 1:2
        for k = 1:12
            clearvars temp_s 
            temp_s = squeeze(po4_in(i,:,k));
            temp_z = zz;
            temp_s(isnan(temp_s)) = interp1( temp_z(~isnan(temp_s)), temp_s(~isnan(temp_s)), temp_z(isnan(temp_s))); 
            po4_in(i,:,k) = temp_s;
        end
end

for i = 1:12
    clearvars tem tem_interp
    tem = squeeze(po4_in(:,:,i));
    tem_interp = griddata(x_in,z_in,tem,lon_rho_1, z_1,'nearest');
    po4_input(:,:,i) = tem_interp;
end

find(isnan(po4_input)==1)

% po4_input_sn = ones(size(t_s,1),size(t_s,2),size(t_s,3));
clearvars po4_input_we
po4_input_we = ones(size(t_w,1),size(t_w,2),size(t_w,3));

for i = 1:12
    for j = 1:20
        po4_input_we(:,j,i) = squeeze(po4_input_we(:,j,i)) .* squeeze(po4_in(2,j,i));
    end
end

save(['boundary_bio_input_to06to15_',num2str(t_i),'regime.mat'],'*_input*');
save(['boundary_bio_input_for_ini_to06to15_',num2str(t_i),'regime.mat'],'*_in');
save(['check_boundary_no3_input_to06to15_',num2str(t_i),'regime.mat'],'y_s_no3_205','y_b_no3_205');
end
