%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc;

list = dir('*.nc');

for i = 1:length(list)
clearvars chl_a chla H3
    if i == 1
        lon = ncread(list(i).name,'lon');
        lat = ncread(list(i).name,'lat');

        [Y X]=meshgrid(lat,lon);

        %south only;
        % find(X >= 127.5700 & X <= 128.1900 & Y >= 34.5920 & Y <= 35.1000);
%         in_co=find(X < 127.5700 | X > 128.1900 | Y < 34.5920 | Y > 35.1000);
%         chl_a(in_co)=NaN;

        % yellow & south 
        % pcolor(X(7150:7500,1100:1500),Y(7150:7500,1100:1500),chl_a(7150:7500,1100:1500));
        % shading flat; colorbar;

        % new
        x_lo = 7383-1:7398;
        y_lo = 1321-3:1331;

        x1=X(x_lo,y_lo);
        y1=Y(x_lo,y_lo);

    end

chl_a = ncread(list(i).name,'chl_ocx');


%check boundary
% chla = chl_a;
% chla(in_co)=NaN;
% % old 
% x_lo=7383:7398;
% y_lo=1321:1331;


c1(:,:,i)=chl_a(x_lo,y_lo);

% pcolor(X(7383-1:7398+1,1321-1:1331+1),Y(7383-1:7398+1,1321-1:1331+1),chl_a(7383-1:7398+1,1321-1:1331+1));
% shading flat; colorbar;
% 
% pcolor(X(7383-1:7398,1321-1:1331+1),Y(7383-1:7398+1,1321-1:1331+1),chla(7383-1:7398+1,1321-1:1331+1));
% shading flat; colorbar;
% 
% figure;
% pcolor(X(x_lo,y_lo),Y(x_lo,y_lo),chla(x_lo,y_lo));
% shading flat; colorbar;

figure; hold on;
pcolor(x1,y1,squeeze(c1(:,:,i)));
shading flat; grid on; set(gca,'fontsize',13);
H3=colorbar;
set(get(H3,'title'),'string','Chl.a(mg/m^3)','fontsize',15,'fontweight','bold');
xlabel('lon','fontsize',13);
ylabel('lat','fontsize',13);
title([list(i).name],'fontsize',13)
set(gca,'fontweight','bold'); caxis([0 10]);
xlim([min(min(x1)) max(max(x1))]); ylim([min(min(y1)) max(max(y1))]);
hold off;


end
% figure;
% pcolor(chl_a); xlim([x_lo(1)-10 x_lo(end)+10]); ylim([y_lo(1)-10 y_lo(end)+10]); 
% pcolor(chl_a(7383:7398,1321:1331));