close all; clear; clc;

x_ref = 115.00 : 2 :164.00;
y_ref = 15.00 : 2 : 52.00;

[y_ref_m x_ref_m]=meshgrid (y_ref, x_ref);

load NIFS_bio_trend_1976_2018.mat

%% numyear,numline,numsta,nummon,1:length(stddepth))
line=[107 106 105 104 103 102 209 208 ...
    207 206 400 205 204 203 317 316 315 314 313 ...
    312 311 310 309 308 307]; 


size(NIFS6.no3)

size(NIFS6.lon)
size(NIFS6.lat)


t_year= 1976:2018;

%  
% NIFS6.lat(NIFS6.lat==0)=NaN;
% NIFS6.lon(NIFS6.lon==0)=NaN;

x_nfis= squeeze(nanmean(nanmean(nanmean(NIFS6.lon,1),4),5));
y_nfis= squeeze(nanmean(nanmean(nanmean(NIFS6.lat,1),4),5));

% no_match = find(y_nfis ~= squeeze(NIFS6.lon(1,:,:,1,1)))
% first=squeeze(NIFS6.lon(1,:,:,1,1))
% x_nfis(no_match)


% nearest point
% k = 0
for i = 1:size(y_nfis,1)
    for j=1:size(y_nfis,2)
clearvars temp_2d near_temp
% k=k+1
if isnan(x_nfis(i,j)) == 0 & isnan(y_nfis(i,j)) == 0
    temp_2d = sqrt((x_ref_m - x_nfis(i,j)).^2 + (y_ref_m - y_nfis(i,j)).^2);
%     temp_2d = temp_2d .* (mask./mask);
    near_temp = find(nanmin(nanmin(temp_2d))==temp_2d);
    near_point(i,j) = near_temp(1);
%     dist_2d(:,:,k) = temp_2d;
elseif isnan(x_nfis(i,j)) == 1 | isnan(y_nfis(i,j))==1
    near_point(i,j) = NaN;
end
    end
end

return

for i = 1:size(y_nfis,1)
    for j=1:size(y_nfis,2)
    clearvars row col
    if isnan(near_point(i,j))==0
    [row,col]=ind2sub(size(y_nfis),near_point(i,j)); %% 1d to 2d index
    near_point_x(i,j) = row; 
    near_point_y(i,j) = col; 
    else isnan(near_point(i,j))==1
    near_point_x(i,j) = NaN; 
    near_point_y(i,j) = NaN; 
    end
    end
end


NIFS6_po4 = NaN(length(t_year), size(x_ref_m,1), size(x_ref_m,2), length(2:2:12), size(NIFS6.lon,5)); 
NIFS6_no3 = NaN(length(t_year), size(x_ref_m,1), size(x_ref_m,2), length(2:2:12), size(NIFS6.lon,5)); 
NIFS6_si = NaN(length(t_year), size(x_ref_m,1), size(x_ref_m,2), length(2:2:12), size(NIFS6.lon,5)); 
NIFS6_do = NaN(length(t_year), size(x_ref_m,1), size(x_ref_m,2), length(2:2:12), size(NIFS6.lon,5)); 
% 
% near_point_x=reshape(near_point_2d(:,1),size(x_nfis,1),size(x_nfis,2));
% near_point_y=reshape(near_point_2d(:,2),size(x_nfis,1),size(x_nfis,2));
% n=0;
% n = n+1;
for i = 1:length(t_year) %year
     for l = 1:length(2:2:12) %month
         for j = 1:size(NIFS6.lon,2) %numline
             for k = 1:size(NIFS6.lon,3) %numsta
                  for m = 1:size(NIFS6.lon,5) %depth  
                    
                if isnan(near_point_x(j,k)) == 0 & isnan(near_point_y(j,k)) == 0
%                     disp(l)
                    NIFS6_po4(i,near_point_x(j,k),near_point_y(j,k),l,m) = nanmean([squeeze(NIFS6_po4(i,near_point_x(j,k),near_point_y(j,k),l,m)); squeeze(NIFS6.po4(i,j,k,l,m));]);
                    NIFS6_no3(i,near_point_x(j,k),near_point_y(j,k),l,m) = nanmean([squeeze(NIFS6_no3(i,near_point_x(j,k),near_point_y(j,k),l,m)); squeeze(NIFS6.no3(i,j,k,l,m));]);
                    NIFS6_si(i,near_point_x(j,k),near_point_y(j,k),l,m) =nanmean([squeeze(NIFS6_si(i,near_point_x(j,k),near_point_y(j,k),l,m)); squeeze(NIFS6.si(i,j,k,l,m));]);
                    NIFS6_do(i,near_point_x(j,k),near_point_y(j,k),l,m) = nanmean([squeeze(NIFS6_do(i,near_point_x(j,k),near_point_y(j,k),l,m)); squeeze(NIFS6.do(i,j,k,l,m));]);
%                 elseif  isnan(near_point_2d(n,1)) == 1 & isnan(near_point_2d9s(n,2)) == 1
%                     NIFS6_po4(i,near_point_2d(n,1),near_point_2d(n,2),l,m) = NaN;
%                     NIFS6_no3(i,near_point_2d(n,1),near_point_2d(n,2),l,m) = NaN;
%                     NIFS6_si(i,near_point_2d(n,1),near_point_2d(n,2),l,m) = NaN;
%                     NIFS6_do(i,near_point_2d(n,1),near_point_2d(n,2),l,m) = NaN;
                end
                
                end
            end
        end
    end
end

return

for i = 1:length(t_year) %year
    for l = 1:length(2:2:12) %month
        if size(find(isnan(NIFS6_po4(i,:,:,l,1))==0)) > 0
        disp(['i =' num2str(i) ', l=' num2str(l)])
        end
    end
end

for i = 1:length(t_year) %year
    for l = 1:length(2:2:12) %month
        if length(find(isnan(NIFS6_no3(i,:,:,l,1))==0)) > 0
        disp(['i =' num2str(i) ', l=' num2str(l)])
        end
    end
end

for i = 1:length(t_year) %year
    for l = 1:length(2:2:12) %month
        if length(find(isnan(NIFS6.no3(i,:,:,l,1))==0)) > 0
        disp(['i =' num2str(i) ', l=' num2str(l)])
        end
    end
end

for i = 1:length(t_year) %year
    for l = 1:length(2:2:12) %month
        if size(find(isnan(NIFS6.temp(i,:,:,l,1))==0)) > 0
        disp(['i =' num2str(i) ', l=' num2str(l)])
        end
    end
end

find(isnan(NIFS6_no3)==0)

pcolor(x_ref_m, y_ref_m, squeeze(NIFS6_no3(19,:,:,3,1))); shading flat;

pcolor(x_ref_m, y_ref_m, squeeze(NIFS6_no3(14,:,:,1,1))); shading flat;


pcolor(x_ref_m, y_ref_m, squeeze(NIFS6_no3(21,:,:,1,1))); shading flat;

pcolor(x_ref_m, y_ref_m, squeeze(NIFS6_po4(1,:,:,3,1))); shading flat;

pcolor(x_ref_m, y_ref_m, squeeze(NIFS6_po4(1,:,:,1,1))); shading flat;
