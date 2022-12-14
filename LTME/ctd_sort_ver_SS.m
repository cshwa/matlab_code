% matlab 으로 density 계산과 함께 sorting 된 파일 matlab으로 그림그려보기 위해 파일 처리
% CTD_processing_ver_SS.mat 파일 실행 후 다음 과정에 해당함

clc; clear all; close all;

data = load('pm_ebb_all.dat');

% 관측한 station 수
station = 30;
st = 01; en = 30; % flood
k = 1;
n = 2;
m = abs(min(data(:,3)));
ctd_en = 36;
m = ctd_en;
sal = zeros(m+2,station);
den = zeros(m+2,station);
temp = zeros(m+2,station);
for i = 1:length(data(:,3))-1
    
    if (data(i,8) == data(i+1,8))        
        temp(n,k) = data(i,4);
        sal(n,k) = data(i,6);
        den(n,k) = data(i,7); 
        n = n + 1;
    else 
        temp(n,k) = data(i,4);
        sal(n,k) = data(i,6);
        den(n,k) = data(i,7);         
        sal(1,k) = data(i,8);
        den(1,k) = data(i,8);
        bottom(k,2) = sal(n,k);
        bottom(k,3) = den(n,k);
        bottom(k,1) = temp(n,k);
        k = k + 1;
        n = 2;
    end 
    
end
temp(n,k) = data(i+1,4);
sal(n,k) = data(i+1,6);
den(n,k) = data(i+1,7); 

temp(1,k) = data(i+1,4);
sal(1,k) = data(i+1,8);
den(1,k) = data(i+1,8);        
        
temp(temp == 0) = NaN;
sal(sal == 0) = NaN;
den(den == 0) = NaN;

% ms(2,:) = nanmean(sal(2:end,:),1);
% md(2,:) = nanmean(den(2:end,:),1);

% 첫번째 행은 각 관측점 간 거리 
temp(1,1) = 0;
sal(1,1) = 0;
den(1,1) = 0;
% ms(1,:) = sal(1,:);
% md(1,:) = den(1,:);
%% surfer 그림용 bottom 값
 bottom(k,2) =data(i+1,6);
 bottom(k,3) = data(i+1,7);
 bottom(k,1) = data(i+1,4);

directory = 'D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_winter_ys';
filename = 'location.xlsx';
filename = fullfile(directory,filename);
loc =  xlsread(filename);
loc = loc(st:en,7:8);

bottom = [loc(:,2), loc(:,1), bottom];
save bottom.dat bottom -ascii

surface = [loc(:,2), loc(:,1), temp(2,:)',sal(2,:)',den(2,:)'];
save surface.dat surface -ascii
%%
save neap_flood_den.dat den -ascii
save neap_flood_sal.dat sal -ascii
save neap_flood_temp.dat temp -ascii
%%
ctd_en = 36;
figure()
subplot(311)
contourf(flipud(temp(2:ctd_en,:))); shading flat; 
colorbar; 
caxis([7 10]);
text(st+1,5,'Temperature','fontsize',12);
ylabel('Depth (m)','fontsize',12);
set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[0 35]);
set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
% grid on; 

subplot(312)
contourf(flipud(sal(2:ctd_en,:))); shading flat;
colorbar; colorbar; caxis([0 30]);
text(st+1,5,'Salinity','fontsize',12);
ylabel('Depth (m)','fontsize',12);
set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[0 35]);
set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);

subplot(313)
contourf(flipud(den(2:ctd_en,:))); shading flat;
colorbar; colorbar; caxis([1000 1027]);
text(st+1,5,'Density','fontsize',12);
ylabel('Depth (m)','fontsize',12);
xlabel('Station','fontsize',12);
set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[0 35]);
set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);

%% temp, sal, den depth mean with variation
ctd_en = 36;
figure()
subplot(3,1,1)
y =  nanmean(temp(2:ctd_en,:)); 
x = 1:station; 
e = nanstd(temp(2:ctd_en,:)); 
errorbar(x, y,e,'b.','markersize',20');
set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
set(gca,'ylim',[7 10]);
title('Depth mean temperature ','fontsize',14);
ylabel('℃','fontsize',14);

subplot(312)
y =  nanmean(sal(2:ctd_en,:)); 
x = 1:station; 
e = nanstd(sal(2:ctd_en,:)); 
errorbar(x, y,e,'b.','markersize',20');
set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
set(gca,'ylim',[-2 40]);
title('Depth mean salinity ','fontsize',14);
ylabel('‰','fontsize',14);


subplot(313)
y =  nanmean(den(2:ctd_en,:)); 
x = 1:station; 
e = nanstd(den(2:ctd_en,:)); 
errorbar(x, y,e,'b.','markersize',20');
set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
set(gca,'ylim',[990 1030]);
title('Depth mean density ','fontsize',14);
ylabel('\rho','fontsize',14);
%%
% 수온의 분산, 수직평균, contour
figure()
% 수심평균과 errorbar
subplot(5,1,1:2)
y =  nanmean(temp(2:end,:)); 
x = 1:station; 
e = nanstd(temp(2:end,:)); 
errorbar(x, y,e,'b.','markersize',20');
set(gca,'xtick',[1:3:25],'xticklabel',[6:3:35],'xlim',[0 25]);
set(gca,'ylim',[7 10]);
title('Temperature (Neap flood)','fontsize',14);
ylabel('Mean T','fontsize',14);

subplot(5,1,3:5)
contourf(flipud(temp(2:end,:))); shading flat;
% colorbar; 
caxis([5 10]);
text(2,5,'Temperature','fontsize',12);
ylabel('Depth (m)','fontsize',14);
xlabel('Number of station','fontsize',12);
set(gca,'xtick',[1:3:25],'xticklabel',[6:3:35],'xlim',[0 25]);
grid on;







