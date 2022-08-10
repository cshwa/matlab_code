% castaway CTD casting st102-100, st1-11
% date: 2015 - 02 - 26~27
% st.8번은 관측이 안되어 idronaut 관측 데이터 참고
% 참고: 연속데이터 불러오는 형식
% for count=1:100
%     eval(['load t' num2str(count) '.dat']);
% end

clc; clear all; close all;


dl = 14; % 관측한 데이터 수
% high 시기 data load
h_temp = 32768*ones(50,dl);
h_sal = 32768*ones(50,dl);
h_den = 32768*ones(50,dl);
dep = ones(10,1);
for count=1:dl
    eval(['load D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\mat\processed\high/h' num2str(count) '.mat']);
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

% ebb 시기 data load
e_temp = 32768*ones(50,dl);
e_sal = 32768*ones(50,dl);
e_den = 32768*ones(50,dl);
dep = ones(10,1);

for count=1:dl
    eval(['load D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\mat\processed\ebb/e' num2str(count) '.mat']);
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


% low 시기 data load
l_temp = 32768*ones(50,dl);
l_sal = 32768*ones(50,dl);
l_den = 32768*ones(50,dl);
dep = ones(10,1);

for count=1:dl
    eval(['load D:\SilverStar\data\01_SumJin\2015_02_26\CTD_castaway\mat\processed\low/l' num2str(count) '.mat']);
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




hfig = figure()
set(hfig, 'Position', [100 100 1400 600])

% high
subplot(331)
contourf(flipud(h_temp));
shading flat;
colorbar;
title('high','fontsize',14);
ylabel('depth (m)');
caxis([5 10]);
text(2,15,'temperature');
set(gca, 'ylim',[10 50],'ytick',[10:10:50],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);
subplot(334)
contourf(flipud(h_sal));
shading flat;
colorbar;caxis([10 35]);
text(2,15,'salinity');
set(gca, 'ylim',[10 50],'ytick',[10:10:50],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);
subplot(337)
contourf(flipud(h_den));
shading flat;
xlabel('station');
colorbar;caxis([1000 1030]);
text(2,15,'density');
set(gca, 'ylim',[10 50],'ytick',[10:10:50],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);

% ebb
subplot(332)
contourf(flipud(e_temp));
shading flat;
colorbar;
title('ebb','fontsize',14);
ylabel('depth (m)');
caxis([5 10]);
text(2,15,'temperature');
set(gca, 'ylim',[10 50],'ytick',[10:10:50],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);
subplot(335)
contourf(flipud(e_sal));
shading flat;
colorbar;caxis([10 35]);
text(2,15,'salinity');
set(gca, 'ylim',[10 50],'ytick',[10:10:50],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);
subplot(338)
contourf(flipud(e_den));
shading flat;
xlabel('station');
colorbar;caxis([1000 1030]);
text(2,15,'density');
set(gca, 'ylim',[10 50],'ytick',[10:10:50],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);


% low
subplot(333)
contourf(flipud(l_temp));
shading flat;
colorbar;
title('low','fontsize',14);
ylabel('depth (m)');
caxis([5 10]);
text(2,15,'temperature');
set(gca, 'ylim',[10 50],'ytick',[10:10:50],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);
subplot(336)
contourf(flipud(l_sal));
shading flat;
colorbar;caxis([10 35]);
text(2,15,'salinity');
set(gca, 'ylim',[10 50],'ytick',[10:10:50],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);
subplot(339)
contourf(flipud(l_den));
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

figure()
plot(h,'m.','markersize',15);
hold on;
plot(e,'b.','markersize',15);
plot(l,'g.','markersize',15);
legend('high','ebb','low');
plot(1:0.01:15,0.32,'k');
plot(1:0.01:15,0.15,'k');

save('StrPar20150226.dat','h','e','l','-ascii');


% save density gradient for comparing to other obseved data
hdg(1,:) = nanmean(h_den);
hdg(2,:) = nanmean(e_den);
hdg(3,:) = nanmean(l_den);
save('hdg20150226.dat','hdg','-ascii');

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
hdg(1,:) = nanmean(h_den);
hdg(2,:) = nanmean(e_den);
hdg(3,:) = nanmean(l_den);
save('hdg20150226.dat','hdg','-ascii');
