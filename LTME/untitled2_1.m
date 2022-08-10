close all; clc; clear;
load interp_dep_1970_fine.mat; h_1970=h(:,1:1017);
load interp_dep_present_fine2.mat; h_present=h*-1;
load coast_line_1980s_final.mat
load coast_line_present.mat
load coast_line_present_island.mat
load coast_line_1980s_island.mat

Xi=[127.5700:0.0005:128.1900]; % 0.0005 - 100m resolution
Yi=[34.5920:0.0005:35.1000];

[Y X]=meshgrid(Yi,Xi);

 diff=h_present-(h_1970+2.713);
% diff=(h_present+2.759)-(h_1970+2.713); %recover depth to add MSL

figure; pcolor(X,Y,diff);shading flat;hold on;h=colorbar;
% plot(lon,lat,'color','r','linew',2);
% plot(p_lon,p_lat,'color','k','linew',1.5);
for i=1:length(p_island_lon)
patch(p_island_lon{i},p_island_lat{i},[192/255 192/255 192/255]);
end
for i=1:length(island_lon)
patch(island_lon{i},island_lat{i},[222/255 184/255 135/255]);
end
ylabel(h,'diff. depth(m)');set(gca,'fontsize',17,'fontweight','bold'); 
hold off;
return

% abs index
index_10=find(abs(diff) < 10);
index_1=find(abs(diff)<1);
index_5=find(abs(diff)<5);
index_3=find(abs(diff)<3);
index_2=find(abs(diff)<2);

% no abs index
no_abs_2 = find(diff < 2);
no_abs_0 = find(diff > 0);
no_abs_0_p = find(diff < 0);

diff_10=diff; diff_1=diff; diff_5=diff; diff_3=diff; diff_2=diff;
no_diff_2=diff; no_diff_0=diff; no_diff_0_p=diff;

diff_10(index_10) = NaN; diff_5(index_5) = NaN; diff_1(index_1) = NaN;
diff_3(index_3) = NaN; diff_2(index_2) = NaN; no_diff_2(no_abs_2) = NaN;
no_diff_0(no_abs_0) = NaN; no_diff_0_p(no_abs_0_p) = NaN;

find(isnan(diff)==0)


% plot diff bigger than 10m
figure; pcolor(X,Y,diff_10);shading flat; h=colorbar;; hold on;
for i=1:length(p_island_lon)
patch(p_island_lon{i},p_island_lat{i},[192/255 192/255 192/255]);
end
for i=1:length(island_lon)
patch(island_lon{i},island_lat{i},[222/255 184/255 135/255]);
end
ylabel(h,'diff. depth(m)');set(gca,'fontsize',17,'fontweight','bold'); 
hold off;


% plot diff bigger than 5m
figure; pcolor(X,Y,diff_5);shading flat; h=colorbar;; hold on;
for i=1:length(p_island_lon)
patch(p_island_lon{i},p_island_lat{i},[192/255 192/255 192/255]);
end
for i=1:length(island_lon)
patch(island_lon{i},island_lat{i},[222/255 184/255 135/255]);
end
ylabel(h,'diff. depth(m)');set(gca,'fontsize',17,'fontweight','bold'); 
hold off;


% plot diff bigger than 3m
figure; pcolor(X,Y,diff_3);shading flat; h=colorbar;; hold on;
for i=1:length(p_island_lon)
patch(p_island_lon{i},p_island_lat{i},[192/255 192/255 192/255]);
end
for i=1:length(island_lon)
patch(island_lon{i},island_lat{i},[222/255 184/255 135/255]);
end
ylabel(h,'diff. depth(m)');set(gca,'fontsize',17,'fontweight','bold'); 
hold off;


% plot diff bigger than 2m
figure; pcolor(X,Y,diff_2);shading flat; h=colorbar;; hold on;
for i=1:length(p_island_lon)
patch(p_island_lon{i},p_island_lat{i},[192/255 192/255 192/255]);
end
for i=1:length(island_lon)
patch(island_lon{i},island_lat{i},[222/255 184/255 135/255]);
end
ylabel(h,'diff. depth(m)');set(gca,'fontsize',17,'fontweight','bold'); 
hold off;


% plot diff bigger than 1m
figure; pcolor(X,Y,diff_1);shading flat; h=colorbar;; hold on;
for i=1:length(p_island_lon)
patch(p_island_lon{i},p_island_lat{i},[192/255 192/255 192/255]);
end
for i=1:length(island_lon)
patch(island_lon{i},island_lat{i},[222/255 184/255 135/255]);
end
ylabel(h,'diff. depth(m)');set(gca,'fontsize',17,'fontweight','bold'); 
hold off;
return

%%%% no abs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot diff bigger than 2m (no_abs)
figure; pcolor(X,Y,no_diff_2);shading flat; h=colorbar;; hold on;
% for i=1:length(p_island_lon)
% patch(p_island_lon{i},p_island_lat{i},[192/255 192/255 192/255]);
% end
for i=1:length(island_lon)
patch(island_lon{i},island_lat{i},[222/255 184/255 135/255]);
end
ylabel(h,'diff. depth(m)');set(gca,'fontsize',17,'fontweight','bold'); 
hold off;

% plot diff smaller than 0m (no_abs)
figure; pcolor(X,Y,no_diff_0);shading flat; h=colorbar;; hold on;
% for i=1:length(p_island_lon)
% patch(p_island_lon{i},p_island_lat{i},[192/255 192/255 192/255]);
% end
for i=1:length(island_lon)
patch(island_lon{i},island_lat{i},[222/255 184/255 135/255]);
end
ylabel(h,'diff. depth(m)');set(gca,'fontsize',17,'fontweight','bold'); 
hold off;

% plot diff bigger than 0m (no_abs)
figure; pcolor(X,Y,no_diff_0_p);shading flat; h=colorbar;; hold on;
% for i=1:length(p_island_lon)
% patch(p_island_lon{i},p_island_lat{i},[192/255 192/255 192/255]);
% end
for i=1:length(island_lon)
patch(island_lon{i},island_lat{i},[222/255 184/255 135/255]);
end
ylabel(h,'diff. depth(m)');set(gca,'fontsize',17,'fontweight','bold'); 
hold off;



