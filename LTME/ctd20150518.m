% castaway CTD casting st102-100, st1-11
% date: 2015 - 03 - 06
% 참고: 연속데이터 불러오는 형식
% for count=1:100
%     eval(['load t' num2str(count) '.dat']);
% end

clc; clear all; close all;
% 
% directory = 'D:\SilverStar\data\01_SumJin\2015_03_06\excel\processed\high';
% filename = 'h1.csv';
% filename = fullfile(directory,filename);
% data =  xlsread(filename);
% lat = data(3,2);
% lon = data(4,2);
% dep = data(23:end,2);
% temp = data(23:end,3);
% sal = data(23:end,6);
% den = data(23:end,8);


dl = 12; % 관측한 데이터 수
dd = 54; % depth 설정
% high 시기 data load
h_temp = 32768*ones(dd,dl);
h_sal = 32768*ones(dd,dl);
h_den = 32768*ones(dd,dl);
dep = ones(10,1);
for count=1:dl
    eval(['load D:\SilverStar\data\01_SumJin\2015_05_18\high/h' num2str(count) '.mat']);
    loc(count,1) = LatitudeStart;
    loc(count,2) = LongitudeStart;
    l = length(Depth);
    h_temp(1:l,count) = Temperature;
    h_sal(1:l, count) = Salinity;
    h_den(1:l, count) = Density;
    if length(Depth) >= length(dep)
        dep = Depth;
    end
end
h_temp(find(h_temp==32768)) = NaN;
h_sal(find(h_sal==32768)) = NaN;
h_den(find(h_den==32768)) = NaN;

spring_h_den = h_den;
spring_h_temp = h_temp;
spring_h_sal = h_sal;

%% figure - depth mean density
% dgd: density gradient depth
dgd = 54;

figure()
% high
[a b] = size(spring_h_den(1:dgd,:));
plot(1:b, nanmean(spring_h_den(1:dgd,:)),'r.','markersize',20);
hold on;
x = [1:b];
fd = nanmean(spring_h_den(1:dgd,:));
p = polyfit(x,fd,1);
y = p(1)*x+p(2);
plot(x,y,'r-');
axis([0 13 1005 1023]);
ylabel('Density','fontsize',14);
xlabel('Station','fontsize',14);
legend('Neap','High','Spring','High');

%% draw temp, sal, den during high/ebb/low tide
hfig = figure()
set(hfig, 'Position', [100 100 1400 600])

% high
subplot(3,4,1)
contourf(flipud(h_temp));
shading flat;
colorbar;
title('high','fontsize',14);
ylabel('depth (m)');
caxis([15 20]);
text(2,15,'temperature');
set(gca, 'ylim',[10 50],'ytick',[10:10:50],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);
subplot(3,4,5)
contourf(flipud(h_sal));
shading flat;
colorbar;caxis([10 35]);
text(2,15,'salinity');
set(gca, 'ylim',[10 50],'ytick',[10:10:50],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);
subplot(3,4,9)
contourf(flipud(h_den));
shading flat;
xlabel('station');
colorbar;caxis([1000 1030]);
text(2,15,'density');
set(gca, 'ylim',[10 50],'ytick',[10:10:50],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);

% ebb
subplot(3,4,2)
contourf(flipud(e_temp));
shading flat;
colorbar;
title('ebb','fontsize',14);
ylabel('depth (m)');
caxis([15 20]);
text(2,15,'temperature');
set(gca, 'ylim',[10 50],'ytick',[10:10:50],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);
subplot(3,4,6)
contourf(flipud(e_sal));
shading flat;
colorbar;caxis([10 35]);
text(2,15,'salinity');
set(gca, 'ylim',[10 50],'ytick',[10:10:50],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);
subplot(3,4,10)
contourf(flipud(e_den));
shading flat;
xlabel('station');
colorbar;caxis([1000 1030]);
text(2,15,'density');
set(gca, 'ylim',[10 50],'ytick',[10:10:50],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);

% low
subplot(3,4,3)
contourf(flipud(l_temp));
shading flat;
colorbar;
title('low','fontsize',14);
ylabel('depth (m)');
caxis([15 20]);
text(2,15,'temperature');
set(gca, 'ylim',[10 50],'ytick',[10:10:50],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);
subplot(3,4,7)
contourf(flipud(l_sal));
shading flat;
colorbar;caxis([10 35]);
text(2,15,'salinity');
set(gca, 'ylim',[10 50],'ytick',[10:10:50],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);
subplot(3,4,11)
contourf(flipud(l_den));
shading flat;
xlabel('station');
colorbar;caxis([1000 1030]);
text(2,15,'density');
set(gca, 'ylim',[10 50],'ytick',[10:10:50],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);

% flood
subplot(3,4,4)
contourf(flipud(f_temp));
shading flat;
colorbar;
title('Flood','fontsize',14);
ylabel('Depth (m)');
caxis([15 20]);
text(2,15,'temperature');
set(gca, 'ylim',[10 50],'ytick',[10:10:50],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);
subplot(3,4,8)
contourf(flipud(f_sal));
shading flat;
colorbar;caxis([10 35]);
text(2,15,'salinity');
set(gca, 'ylim',[10 50],'ytick',[10:10:50],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);
subplot(3,4,12)
contourf(flipud(f_den));
shading flat;
xlabel('station');
colorbar;caxis([1000 1030]);
text(2,15,'density');
set(gca, 'ylim',[10 50],'ytick',[10:10:50],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);


%% stratification parameter
h = (max(h_sal)-min(h_sal))./nanmean(h_sal);
e = (max(e_sal)-min(e_sal))./nanmean(e_sal);
l = (max(l_sal)-min(l_sal))./nanmean(l_sal);
f = (max(f_sal)-min(f_sal))./nanmean(f_sal);


hfig = figure();
set(hfig, 'Position', [100 100 600 300])

plot(h,'k.','markersize',15);
hold on;
plot(e,'b.','markersize',15);
plot(l,'c.','markersize',15);
plot(f,'r.','markersize',15);
text(0.1,0.15,'0.15');
text(0.1,0.32,'0.32');
ylabel('\deltaS/<S>','fontsize',14);
xlabel('seaward                      station                     landward','fontsize',14);
legend('high','ebb','low','flood');
plot(1:0.01:15,0.32,'k');
plot(1:0.01:15,0.15,'k');

sp(4,:) = h;
sp(1,:) = e;
sp(2,:) = l;
sp(3,:) = f;
save('sp20150306.dat','sp','-ascii');

%% horizontal density gradient
diff_den(1) = nanmean(h_den(:,1)) - nanmean(h_den(:,end));
diff_den(2) = nanmean(e_den(:,1)) - nanmean(e_den(:,end));
diff_den(3) = nanmean(l_den(:,1)) - nanmean(l_den(:,end));
diff_den(3) = nanmean(f_den(:,1)) - nanmean(f_den(:,end));

distance = 11.90*1000;
NdiffDen = [7.96404415584425,7.81162261904751,5.86997285714301];
NdiffDen = NdiffDen/distance;
SdiffDen = diff_den/distance;

figure()
plot(NdiffDen,'b.','markersize',20);
hold on;
plot(SdiffDen,'r.','markersize',20);
legend('Neap','Spring');
ylabel('\delta\rho/\deltax','fontsize',14);
set(gca,'xtick',[1:1:3],'xlim',[0 4],'xticklabel',{'high','ebb','low'},'fontsize',14);

% save density gradient for comparing to other obseved data
hdg(4,:) = nanmean(h_den);
hdg(1,:) = nanmean(e_den);
hdg(2,:) = nanmean(l_den);
save('hdg20150306.dat','hdg','-ascii');

distance = [0 1.261076124 2.40147083 3.535927832 4.609061223 ...
     5.578879746 6.554433813 7.63283546 8.821608365 9.959095513 ...
     10.86165037 11.68774175 12.64185767 13.50797955 ];
save('distance20150306.dat','distance','-ascii');




