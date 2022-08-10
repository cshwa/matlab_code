%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make boundary file
% % %  1 micro gram to mili gram  = 0.001
%%% mg/L to mM/m^3;
% sumjin_do = yp_w_do'*0.7*44.661;
% sumjin_chl = yp_w_chl';
% sumjin_no3 = yp_w_no3'*1000/14;
% sumjin_nh4 = yp_w_nh4'*1000/14;
close all; clear all; clc;
%%
% from polynomial_fitting_on_ecological_variable_v3_3sig(boundary)_chl.m
%  s: t_regime = length(1997:2005)*12 + 4;
%  no_interp-> s: 2006-05, b: no regime
load chl_koem_input.mat %  y_s_chl, y_b_chl (mon) 

%% 
% from polynomial_fitting_on_ecological_variable_v3_3sig(boundary).m
% there are 2 regimes for nh4 (KOEM 7st.)
% t_regime = length(1997:2008)*12 + 5;
% no_interp-> (s: 2010-02, b: 2010-02)
clearvars *_nh4
load nh4_koem_input.mat % y_s_nh4 & y_b_nh4 (mon) which means 1st regime
y_s_nh4 = y_s_nh4 .*0.001.*1000./14; %%% ug/L to mM/m^3;
y_b_nh4 = y_b_nh4 .*0.001.*1000./14; %%% ug/L to mM/m^3;
% extract negative value to 0
 minus_nh4_s=find(y_s_nh4 < 0); 
 minus_nh4_b=find(y_b_nh4 < 0); 

 y_s_nh4(minus_nh4_s) = 0.01;
 y_b_nh4(minus_nh4_b) = 0.01;
 
 figure; hold on;
 plot(y_s_nh4)
hold on; plot(y_b_nh4,'r')
 
%% it's kodc
% from plot_KODC_6_20501_polynomial.m
% DO is reconstructed
clearvars do_re_b
load do_b_kodc_input_205_01.mat % do_re_b (day)
do_re_b_205 = do_re_b .*0.7.*44.661; %%% mg/L to mM/m^3;

% from plot_KODC_6_40616_polynomial.m
clearvars do_re_b
load do_b_kodc_input_400_16.mat % do_re_b (day)
do_re_b_400 = do_re_b .*0.7.*44.661; %%% mg/L to mM/m^3;
%%
% from make_do_input.m
% DO is reconstructed
clearvars do_re_s
load do_s_kodc_input_205_01.mat % do_re_b (day)
do_re_s_205 = do_re_s .*0.7.*44.661; %%% mg/L to mM/m^3;

% from make_do_input.m
clearvars do_re_s
load do_s_kodc_input_400_16.mat % do_re_s (day)
do_re_s_400 = do_re_s .*0.7.*44.661; %%% mg/L to mM/m^3;


%% 
% 20501 : PO4 is no regime on surf and bot.
clearvars yp_w_po4_04 yp_w_po4_04_b
load po4_koem_input_204_01.mat % 'yp_w_po4_04','yp_w_po4_04_b'
y_s_po4_205 = yp_w_po4_04 .*0.001.*1000./30.973762; %%% ug/L to mM/m^3;
y_b_po4_205 = yp_w_po4_04_b .*0.001.*1000./30.973762; %%% ug/L to mM/m^3;

% 40016 : PO4 is no regime for surf but bot had it (shift 2011-03) 
% load KODC_data_monthly_v5_40016.mat
clearvars yp_w_po4_04* yp_w_po4_04_b*
load po4_koem_input_400_16.mat % yp_w_po4, yp_w_po4_b yp_w_po4_b_1 yp_w_po4_b_2
y_s_po4_400 = yp_w_po4 .*0.001.*1000./30.973762; %%% ug/L to mM/m^3;
y_b_po4_400= yp_w_po4_b_1 .*0.001.*1000./30.973762; %%% ug/L to mM/m^3;

%%
clearvars  y_s_no3* y_b_no3*
%  from plot_KODC_4_20501_polynomial_partial.m
%  there are 3 regimes for no3 (KODC)
%  % no3-kodc   s: 2000-10-01, 2010-04-01 no_interp->(2000-10-01, 2010-04-01)
%   | b: 2000-09-01, 2010-03-01  no_interp-> (2000-12, 2011-03)
load no3_koem_input_205_01_fix.mat % y_s_no3 &  y_b_no3 (366) 
% y_s_no3_205 = y_s_no3 .*0.001.*1000./14; %%% ug/L to mM/m^3;
% y_b_no3_205 = y_b_no3 .*0.001.*1000./14; %%% ug/L to mM/m^3;

y_s_no3_205 = yp_w_no3_1 .*0.001.*1000./14; %%% ug/L to mM/m^3;
y_b_no3_205 = yp_w_no3_b_1 .*0.001.*1000./14; %%% ug/L to mM/m^3;

%  from plot_KODC_4_40616_polynomial_partial.m
% s: 2002-06, 2010-04 no_interp-> same 
% | b: 2000-12, 2011-03 no_interp-> same 
clearvars  y_s_no3 y_b_no3 y_s_no3_2 y_b_no3_2 y_s_no3_3 y_b_no3_3
load no3_koem_input_400_16.mat % y_s_no3 &  y_b_no3 (366)
% y_s_no3_400 = y_s_no3 .*0.001.*1000./14; %%% ug/L to mM/m^3;
% y_b_no3_400 = y_b_no3 .*0.001.*1000./14; %%% ug/L to mM/m^3;
y_s_no3_400 = y_s_no3 .*0.001.*1000./14; %%% ug/L to mM/m^3;
y_b_no3_400 = yp_w_no3_b_1 .*0.001.*1000./14; %%% ug/L to mM/m^3;


% make 1980~present
k=0
for i = 1980:1980
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

clearvars em em_seanson k kk

kk=[1:12];
for i = 1:length(kk)
    if i == 1
        em_season{i} = 1:sum(eom_d(1:kk(i)));
    else
        em_season{i} = sum(eom_d(1:kk(i)-1))+1:sum(eom_d(1:kk(i)));
    end
end


for i = 1:12
        y_s_no3_205_m(i) = nanmean(squeeze(y_s_no3_205(1,em_season{i})));
        y_b_no3_205_m(i) = nanmean(squeeze(y_b_no3_205(1,em_season{i})));
        y_s_no3_400_m(i) = nanmean(squeeze(y_s_no3_400(1,em_season{i})));
        y_b_no3_400_m(i) = nanmean(squeeze(y_b_no3_400(1,em_season{i})));
        y_s_po4_205_m(i) = nanmean(squeeze(y_s_po4_205(1,em_season{i})));
        y_b_po4_205_m(i) = nanmean(squeeze(y_b_po4_205(1,em_season{i})));
        y_s_po4_400_m(i) = nanmean(squeeze(y_s_po4_400(1,em_season{i})));
        y_b_po4_400_m(i) = nanmean(squeeze(y_b_po4_400(1,em_season{i})));
end

clearvars eom_d
% make 1980~present
    k=0
    for i = 1980:2019
                k=k+1;
        for j = 1:12
            eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
        end
    end
    
for j = 1:length(1980:2019)
    for i = 1:12
        em(j,i) = sum(eom_d(j,1:i))
    if j== 1 & i == 12
        eom_dd(j,:) = em(j,:);
    end
    if j ~= 1 & i == 12
        clearvars pre_d
        pre_d = sum(sum(em(1:j-1,end)));
        eom_dd(j,:)= em(j,:) + pre_d; 
    end
    end
end

ed = sort(reshape(eom_dd,1,12*length(1980:2019)))      
for i = 1:length(ed)
    if i== 1
        ed_f(i) = 1;
    else
        ed_f(i)=ed(i-1)+1;
    end
end

for xx = 1:length(1980:2019)
    for i = 1:12 % xx = year, i = 1:12
            do_re_b_205_m(xx,i) = nanmean(squeeze(do_re_b_205(ed_f(1,(xx-1)*12+i):ed(1,(xx-1)*12+i))));
            do_re_b_400_m(xx,i) = nanmean(squeeze(do_re_b_400(ed_f(1,(xx-1)*12+i):ed(1,(xx-1)*12+i))));
            do_re_s_205_m(xx,i) = nanmean(squeeze(do_re_s_205(ed_f(1,(xx-1)*12+i):ed(1,(xx-1)*12+i))));
            do_re_s_400_m(xx,i) = nanmean(squeeze(do_re_s_400(ed_f(1,(xx-1)*12+i):ed(1,(xx-1)*12+i))));   
    end
end

%filling nan
do_re_s_205_m(1,1)=do_re_s_205_m(1,2);
do_re_b_205_m(1,1)=do_re_b_205_m(1,2);

do_re_s_400_m(1,1)=do_re_s_400_m(1,2);
do_re_b_400_m(1,1)=do_re_b_400_m(1,2);



%% vertical interp.  %time s_rho lon_rho
obc_path ='D:\장기생태\Dynamic\07_boundary_ts\Gwangyang_bry_ykang\Gwangyang_add_uv_bry\';
ncdump([obc_path, 'Gwangyang_new_Y1980.nc'])

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
chl_in = NaN(12,20); chl_in(1:12,20)=y_s_chl; chl_in(1:12,1)=y_b_chl;

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
do_in_s = NaN(2,length(1980:2019),12); do_in_s(1,:,:) = do_re_s_205_m; do_in_s(2,:,:) = do_re_s_400_m;
do_in_b = NaN(2,length(1980:2019),12); do_in_b(1,:,:) = do_re_b_205_m; do_in_b(2,:,:) = do_re_b_400_m;
do_in = NaN(2,20,length(1980:2019),12); do_in(:,1,:,:) = do_in_b; do_in(:,20,:,:) = do_in_s;
x_in = NaN(2,20); x_in(1,:) = 127.8083; x_in(2,:) = 128.2033;
z_in = NaN(2,20); 
zz = [20:-1:1];
for i = 1:20
    z_in(:,i) = zz(i);
end

% vertical interp on KODC st.
for i = 1:2
    for j = 1:length(1980:2019)
        for k = 1:12
            clearvars temp_s 
            temp_s = squeeze(do_in(i,:,j,k));
            temp_z = zz;
            if j < length(1980:2019)
                temp_s(isnan(temp_s)) = interp1( temp_z(~isnan(temp_s)), temp_s(~isnan(temp_s)), temp_z(isnan(temp_s))); 
            elseif j == length(1980:2019)
                temp_s;
            end
            do_in(i,:,j,k) = temp_s;
        end
    end
end

for j = 1:length(1980:2019)
for i = 1:12
    clearvars tem tem_interp
    tem = squeeze(do_in(:,:,j,i));
    tem_interp = griddata(x_in,z_in,tem,lon_rho_1, z_1,'nearest');
    do_input(:,:,j,i) = tem_interp;
end
end

find(isnan(do_input)==1)

% do_input_sn = ones(size(t_s,1),size(t_s,2),size(t_s,3));
clearvars do_input_we
do_input_we = ones(size(t_w,1),size(t_w,2),length(1980:2019),size(t_w,3));

for i = 1:12
    for k = 1:length(1980:2019)
    for j = 1:20
        do_input_we(:,j,k,i) = squeeze(do_input_we(:,j,k,i)) .* squeeze(do_in(2,j,k,i));
    end
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

save('boundary_bio_input_1994_v2.mat','*_input*');
save('boundary_bio_input_for_ini_1994_v2.mat','*_in');

save('check_boundary_no3_input_1994_v2.mat','y_s_no3_205','y_b_no3_205')


