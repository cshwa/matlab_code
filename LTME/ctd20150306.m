% castaway CTD casting st102-100, st1-11
% date: 2015 - 03 - 06
% 참고: 연속데이터 불러오는 형식
% for count=1:100
%     eval(['load t' num2str(count) '.dat']);
% end

% clc; clear all; close all;
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


dl = 14; % 관측한 데이터 수
% high 시기 data load
h_temp = 32768*ones(50,dl);
h_sal = 32768*ones(50,dl);
h_den = 32768*ones(50,dl);
dep = ones(10,1);
for count=1:dl
    eval(['load D:\SilverStar\data\01_SumJin\2015_03_06\mat\processed\high/h' num2str(count) '.mat']);
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

s_h_temp = h_temp;
s_h_sal = h_sal;
s_h_den = h_den;

% ebb 시기 data load
e_temp = 32768*ones(50,dl);
e_sal = 32768*ones(50,dl);
e_den = 32768*ones(50,dl);
dep = ones(10,1);

for count=1:dl
    eval(['load D:\SilverStar\data\01_SumJin\2015_03_06\mat\processed\ebb/e' num2str(count) '.mat']);
    loc(count,1) = LatitudeStart;
    loc(count,2) = LongitudeStart;
    l = length(Depth);
    e_temp(1:l,count) = Temperature;
    e_sal(1:l, count) = Salinity;
    e_den(1:l, count) = Density;
    if length(Depth) >= length(dep)
        dep = Depth;
    end
end
e_temp(find(e_temp==32768)) = NaN;
e_sal(find(e_sal==32768)) = NaN;
e_den(find(e_den==32768)) = NaN;

s_e_temp = e_temp;
s_e_sal = e_sal;
s_e_den = e_den;

% low 시기 data load
l_temp = 32768*ones(50,dl);
l_sal = 32768*ones(50,dl);
l_den = 32768*ones(50,dl);
dep = ones(10,1);

for count=1:dl
    eval(['load D:\SilverStar\data\01_SumJin\2015_03_06\mat\processed\low/l' num2str(count) '.mat']);
    loc(count,1) = LatitudeStart;
    loc(count,2) = LongitudeStart;
    l = length(Depth);
    l_temp(1:l,count) = Temperature;
    l_sal(1:l, count) = Salinity;
    l_den(1:l, count) = Density;
    if length(Depth) >= length(dep)
        dep = Depth;
    end
end
l_temp(find(l_temp==32768)) = NaN;
l_sal(find(l_sal==32768)) = NaN;
l_den(find(l_den==32768)) = NaN;

s_l_temp = l_temp;
s_l_sal = l_sal;
s_l_den = l_den;
%% figure - depth mean density
figure()

subplot(121)
[a b] = size(h_den(:,4:end));
plot(1:b, nanmean(h_den(:,4:end)),'r.','markersize',20);
hold on;
x = [1:b];
fd = nanmean(h_den(:,4:end));
p = polyfit(x,fd,1);
y = p(1)*x+p(2);
plot(x,y,'r-');
[a b] = size(l_den(:,4:end));
plot(1:b, nanmean(l_den(:,4:end)),'b.','markersize',20);
% linear poly fit
x = [1:b];
fd = nanmean(h_den(:,4:end));
p = polyfit(x,fd,1);
y = p(1)*x+p(2);
plot(x,y,'r-');
x = [1:b];
ed = nanmean(l_den(:,4:end));
p = polyfit(x,ed,1);
y = p(1)*x+p(2);
plot(x,y,'b-');

legend('high','low');

subplot(122)
[a b] = size(e_den(:,4:end));
plot(1:b, nanmean(e_den(:,4:end)),'r.','markersize',20);
hold on;
x = [1:b];
hd = nanmean(e_den(:,4:end));
p = polyfit(x,hd,1);
y = p(1)*x+p(2);
plot(x,y,'r-');
[a b] = size(l_den(:,4:end));
% plot(1:b, nanmean(l_den(:,4:end)),'b.','markersize',20);
% x = [1:b];
% ld = nanmean(l_den(:,4:end));
% p = polyfit(x,ld,1);
% y = p(1)*x+p(2);
% plot(x,y,'b-');
% legend('high','low');
% plot(6.5,1010:0.1:1024,'k');
% ylim([1010 1024]);
xlim([0 13]);



%% draw temp, sal, den during high/ebb/low tide
hfig = figure()
set(hfig, 'Position', [100 100 1400 600])

% high
subplot(331)
contourf(fliplr(flipud(h_temp)));
shading flat;
% colorbar;
title('High','fontsize',14);
ylabel('Depth (m)','fontsize',14);
caxis([5 10]);
text(2,15,'Temperature','fontsize',14);
set(gca, 'ylim',[10 50],'ytick',[10:10:50],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);

subplot(334)
contourf(fliplr(flipud(h_sal)));
shading flat;
% colorbar;
caxis([10 35]);
text(2,15,'Salinity','fontsize',14);
set(gca, 'ylim',[10 50],'ytick',[10:10:50],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);

subplot(337)
contourf(fliplr(flipud(h_den)));
shading flat;
xlabel('Station','fontsize',14);
% colorbar;
caxis([1000 1030]);
text(2,15,'Density','fontsize',14);
set(gca, 'ylim',[10 50],'ytick',[10:10:50],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);

% ebb
subplot(332)
contourf(fliplr(flipud(e_temp)));
shading flat;
% colorbar;
title('Ebb','fontsize',14);
% ylabel('depth (m)');
caxis([5 10]);
text(2,15,'Temperature','fontsize',14);
set(gca, 'ylim',[10 50],'ytick',[10:10:50],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);

subplot(335)
contourf(fliplr(flipud(e_sal)));
shading flat;
% colorbar;
caxis([10 35]);
text(2,15,'Salinity','fontsize',14);
set(gca, 'ylim',[10 50],'ytick',[10:10:50],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);

subplot(338)
contourf(fliplr(flipud(e_den)));
shading flat;
xlabel('Station','fontsize',14);
% colorbar;
caxis([1000 1030]);
text(2,15,'Density','fontsize',14);
set(gca, 'ylim',[10 50],'ytick',[10:10:50],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);


% low
subplot(333)
contourf(fliplr(flipud(l_temp)));
shading flat;
colorbar;
title('low','fontsize',14);
% ylabel('Depth (m)','fontsize',14));
caxis([5 10]);
text(2,15,'Temperature','fontsize',14);
set(gca, 'ylim',[10 50],'ytick',[10:10:50],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);
subplot(336)
contourf(fliplr(flipud(l_sal)));
shading flat;
colorbar;caxis([10 35]);
text(2,15,'Salinity','fontsize',14);
set(gca, 'ylim',[10 50],'ytick',[10:10:50],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);
subplot(339)
contourf(flipud(fliplr(l_den)));
shading flat;
xlabel('Station','fontsize',14);
colorbar;caxis([1000 1030]);
text(2,15,'Density','fontsize',14);
set(gca, 'ylim',[10 50],'ytick',[10:10:50],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);

%% stratification parameter
h = (max(h_sal)-min(h_sal))./nanmean(h_sal);
e = (max(e_sal)-min(e_sal))./nanmean(e_sal);
l = (max(l_sal)-min(l_sal))./nanmean(l_sal);
% f = (max(f_sal)-min(f_sal))./nanmean(f_sal);


hfig = figure();
set(hfig, 'Position', [100 100 600 300])

plot(h,'k.','markersize',15);
hold on;
plot(e,'b.','markersize',15);
plot(l,'c.','markersize',15);
% plot(f,'r.');
text(0.1,0.15,'0.15');
text(0.1,0.32,'0.32');
ylabel('\deltaS/<S>','fontsize',14);
xlabel('seaward                      station                     landward','fontsize',14);
legend('high','ebb','low');
plot(1:0.01:15,0.32,'k');
plot(1:0.01:15,0.15,'k');

sp(4,:) = h;
sp(1,:) = e;
sp(2,:) = l;
save('sp20150306.dat','sp','-ascii');

%% horizontal density gradient
diff_den(1) = nanmean(h_den(:,1)) - nanmean(h_den(:,end));
diff_den(2) = nanmean(e_den(:,1)) - nanmean(e_den(:,end));
diff_den(3) = nanmean(l_den(:,1)) - nanmean(l_den(:,end));

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




