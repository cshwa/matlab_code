close all; clc; clear;
load interp_dep_1970_fine.mat; h_1970=h(:,1:1017);
load interp_dep_present_fine2.mat; h_present=h*-1;
load interp_dep_present_fine.mat; h_silver=h(:,1:1017);
load coast_line_1980s_final.mat
load coast_line_present.mat
load coast_line_present_island.mat
load coast_line_1980s_island.mat


Xi=[127.5700:0.0005:128.1900]; % 0.0005 - 100m resolution
Yi=[34.5920:0.0005:35.1000];

[Y X]=meshgrid(Yi,Xi);

h_1970=h_1970+2.713; %recover depth to add MSL
diff=h_present-h_1970; 

% abs index
index_10=find(abs(diff) < 10);
index_1=find(abs(diff)<1);
index_5=find(abs(diff)<5);
index_3=find(abs(diff)<3);
index_2=find(abs(diff)<2);

% no abs index
no_abs_1 = find(diff < 1);
no_abs_0 = find(diff > 0);
no_abs_0_p = find(diff < 0);
no_abs_1_p = find(diff > 1);

diff_10=diff; diff_1=diff; diff_5=diff; diff_3=diff; diff_2=diff;
no_diff_1=diff; no_diff_0=diff; no_diff_0_p=diff;

diff_10(index_10) = NaN; diff_5(index_5) = NaN; diff_1(index_1) = NaN;
diff_3(index_3) = NaN; diff_2(index_2) = NaN; no_diff_1(no_abs_1) = NaN;
no_diff_0(no_abs_0) = NaN; no_diff_0_p(no_abs_0_p) = NaN;

% plot diff bigger than 1m (no_abs)
figure; pcolor(X,Y,no_diff_1);shading flat; h=colorbar;; hold on;
% for i=1:length(p_island_lon)
% patch(p_island_lon{i},p_island_lat{i},[192/255 192/255 192/255]);
% end
for i=1:length(island_lon)
patch(island_lon{i},island_lat{i},[222/255 184/255 135/255]);
end
ylabel(h,'diff. depth(m)');set(gca,'fontsize',17,'fontweight','bold'); 
hold off;

return

% picking index for replacing dep.
index_repl=find(X(no_abs_1_p) >= 127.7904); % 127.78645 is criteria of replacing dep.   

diff_r = diff; diff_r(no_abs_1_p(index_repl))=NaN; 

% no_abs_1_p(index_repl) would be index for replacing.
pcolor(X,Y,diff_r); shading flat; 

h_merge_1 = h_1970;
h_merge_1(no_abs_1_p(index_repl)) = h_present(no_abs_1_p(index_repl));

% replacing dep diff.
diff_merg=h_merge_1 - h_present;

% plot replacing dep diff.
figure; pcolor(X,Y,diff_merg);shading flat; h=colorbar;; hold on;
for i=1:length(p_island_lon)
patch(p_island_lon{i},p_island_lat{i},[192/255 192/255 192/255]);
end
for i=1:length(island_lon)
patch(island_lon{i},island_lat{i},[222/255 184/255 135/255]);
end
ylabel(h,'diff. depth(m)');set(gca,'fontsize',17,'fontweight','bold'); 
hold off;

% picking index for replacing dep. from present.(non_nan)
nan_idx_merg1 = find(isnan(h_merge_1)==1);

h_merge_1(nan_idx_merg1)=h_present(nan_idx_merg1);

%fill the sumjin dep. from silver's dep
nan_idx_merg2 = find(isnan(h_merge_1)==1);
h_merge_1(nan_idx_merg2)=h_silver(nan_idx_merg2);

% plot replacing dep.from present.(non_nan)
figure; pcolor(X,Y,h_merge_1);shading flat; h=colorbar;; hold on;
for i=1:length(p_island_lon)
patch(p_island_lon{i},p_island_lat{i},[192/255 192/255 192/255]);
end
for i=1:length(island_lon)
patch(island_lon{i},island_lat{i},[222/255 184/255 135/255]);
end
ylabel(h,'diff. depth(m)');set(gca,'fontsize',17,'fontweight','bold'); 
hold off;


save('merge_1970_depth_gjb.mat','h_merge_1','X','Y');


%%%%
% plot raw.
figure; pcolor(X,Y,h_present);shading flat; h=colorbar;; hold on;
for i=1:length(p_island_lon)
patch(p_island_lon{i},p_island_lat{i},[192/255 192/255 192/255]);
end
% for i=1:length(island_lon)
% patch(island_lon{i},island_lat{i},[222/255 184/255 135/255]);
% end
ylabel(h,'diff. depth(m)');set(gca,'fontsize',17,'fontweight','bold'); 
hold off;

% plot raw.
figure; pcolor(X,Y,h_1970);shading flat; h=colorbar;; hold on;
% for i=1:length(p_island_lon)
% patch(p_island_lon{i},p_island_lat{i},[192/255 192/255 192/255]);
% end
for i=1:length(island_lon)
patch(island_lon{i},island_lat{i},[222/255 184/255 135/255]);
end
ylabel(h,'diff. depth(m)');set(gca,'fontsize',17,'fontweight','bold'); 
hold off;

%%% check model's dep

lon_mo = ncread('grid_gy_v11_s.nc','lon_rho');
lat_mo = ncread('grid_gy_v11_s.nc','lat_rho');
dep_mo = ncread('grid_gy_v11_s.nc','h');

figure; pcolor(lon_mo,lat_mo,dep_mo);shading flat; h=colorbar;; hold on;
for i=1:length(p_island_lon)
patch(p_island_lon{i},p_island_lat{i},[192/255 192/255 192/255]);
end
for i=1:length(island_lon)
patch(island_lon{i},island_lat{i},[222/255 184/255 135/255]);
end
ylabel(h,'diff. depth(m)');set(gca,'fontsize',17,'fontweight','bold'); 
hold off;


min(min(lon_mo))
max(max(lon_mo))
min(min(lat_mo))
max(max(lat_mo))









