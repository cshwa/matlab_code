% high 시기 low 시기 ctd -  idronaut 
% clc; clear all; close all;

% load CTD high data
directory = 'D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\1high';
filename = 'h_den.xlsx';
filename = fullfile(directory,filename);
h_den =  xlsread(filename);

filename = 'h_sal.xlsx';
filename = fullfile(directory,filename);
h_sal =  xlsread(filename);

filename = 'h_temp.xlsx';
filename = fullfile(directory,filename);
h_temp =  xlsread(filename);

% load CTD ebb data
directory = 'D:\SilverStar\data\01_SumJin\2015_02_26\ctd\ebb';
filename = 'i_den.xlsx';
filename = fullfile(directory,filename);
e_den =  xlsread(filename);

filename = 'i_sal.xlsx';
filename = fullfile(directory,filename);
e_sal =  xlsread(filename);

filename = 'i_temp.xlsx';
filename = fullfile(directory,filename);
e_temp =  xlsread(filename);

% load CTD low data
directory = 'D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\3low';
filename = 'l_den.xlsx';
filename = fullfile(directory,filename);
l_den =  xlsread(filename);

filename = 'l_sal.xlsx';
filename = fullfile(directory,filename);
l_sal =  xlsread(filename);

filename = 'l_temp.xlsx';
filename = fullfile(directory,filename);
l_temp =  xlsread(filename);

% load CTD flood data
directory = 'D:\SilverStar\data\01_SumJin\2015_02_26\ctd\flood';
filename = 'i_den.xlsx';
filename = fullfile(directory,filename);
f_den =  xlsread(filename);

filename = 'i_sal.xlsx';
filename = fullfile(directory,filename);
f_sal =  xlsread(filename);

filename = 'i_temp.xlsx';
filename = fullfile(directory,filename);
f_temp =  xlsread(filename);

%%
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
text(2,3,'Temperature','fontsize',14);
set(gca, 'ylim',[0 12],'ytick',[0:3:12],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);

subplot(334)
contourf(fliplr(flipud(h_sal)));
shading flat;
% colorbar;
caxis([10 35]);
text(2,3,'Salinity','fontsize',14);
set(gca, 'ylim',[0 12],'ytick',[0:3:12],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);

subplot(337)
contourf(fliplr(flipud(h_den)));
shading flat;
xlabel('Station','fontsize',14);
% colorbar;
caxis([1000 1030]);
text(2,3,'Density','fontsize',14);
set(gca, 'ylim',[0 12],'ytick',[0:3:12],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);

% ebb
subplot(3,3,2)
contourf(fliplr(flipud(e_temp)));
shading flat;
title('Ebb','fontsize',14);
% colorbar;
caxis([5 10]);
text(2,3,'Temperature','fontsize',14);
set(gca, 'ylim',[0 12],'ytick',[0:3:12],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);

subplot(3,3,5)
contourf(fliplr(flipud(e_sal)));
shading flat;
% colorbar;
caxis([10 35]);
text(2,3,'Salinity','fontsize',14);
set(gca, 'ylim',[0 12],'ytick',[0:3:12],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);

subplot(3,3,8)
contourf(fliplr(flipud(e_den)));
shading flat;
% colorbar;
% colorbar;
caxis([1000 1030]);
text(2,3,'Density','fontsize',14);
xlabel('Station','fontsize',14);
set(gca, 'ylim',[0 12],'ytick',[0:3:12],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);

% low
subplot(333)
contourf(fliplr(flipud(l_temp)));
shading flat;
title('Low','fontsize',14);
colorbar;
caxis([5 10]);
text(2,3,'Temperature','fontsize',14);
set(gca, 'ylim',[0 12],'ytick',[0:3:12],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);

subplot(336)
contourf(fliplr(flipud(l_sal)));
shading flat;
colorbar;
caxis([10 35]);
text(2,3,'Salinity','fontsize',14);
set(gca, 'ylim',[0 12],'ytick',[0:3:12],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);

subplot(339)
contourf(fliplr(flipud(l_den)));
shading flat;
colorbar;
caxis([1000 1030]);
text(2,3,'Density','fontsize',14);
xlabel('Station','fontsize',14);
set(gca, 'ylim',[0 12],'ytick',[0:3:12],'yticklabel',[12:-3:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);
%%

hfig = figure();
set(hfig, 'Position', [100 100 600 300])
h = (max(h_sal)-min(h_sal))./nanmean(h_sal);
e = (max(e_sal)-min(e_sal))./nanmean(e_sal);
l = (max(l_sal)-min(l_sal))./nanmean(l_sal);
f = (max(f_sal)-min(f_sal))./nanmean(f_sal);

plot(h,'m.','markersize',15);
hold on;
% plot(e,'b.','markersize',15);
plot(l,'g.','markersize',15);
% plot(f,'r.','markersize',15);
% legend('high','ebb','low','flood');
plot(1:0.01:15,0.32,'k');
plot(1:0.01:15,0.15,'k');
text(0.1,0.15,'0.15');
text(0.1,0.32,'0.32');
ylabel('\deltaS/<S>','fontsize',14);
xlabel('seaward                      station                     landward','fontsize',14);

diff_den(1) = nanmean(h_den(:,1)) - nanmean(h_den(:,end));
diff_den(2) = nanmean(e_den(:,1)) - nanmean(e_den(:,end));
diff_den(3) = nanmean(l_den(:,1)) - nanmean(l_den(:,end));
diff_den(4) = nanmean(f_den(:,1)) - nanmean(f_den(:,end));
figure()
plot(diff_den,'k.','markersize',20);
