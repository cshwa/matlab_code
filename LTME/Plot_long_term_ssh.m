close all; clear; clc;
[raw txt]=xlsread('Yeosu_UST_TG.xlsx','sheet1','');
[raw1 txt1]=xlsread('Yeosu_UST_TG.xlsx','sheet2','');

SSH=raw(:,6);
SSH_yr=raw1(:,6);



plot(SSH)
hold on; 
 %%% plot trend %%
         x=(1:length(SSH));
         idx = find(isnan(SSH)==1)
         p1 = polyfit(x', SSH,1);
         y1 = polyval(p1,x');
plot(1:length(SSH),y1,'r','linew',2);

        title('TG - detrend');
        set(gca,'linew',2);
        xlabel('Time (date)');
        ylabel('elevation (m)');
        set(gca,'xtick',1:240:length(time));
        set(gca,'xticklabel',tt(1:240:length(time),:));
%         set(gca,'XTickLabelRotation',45);
        set(gca,'fontsize',20);
        
        
        
       figure; plot(SSH_yr); hold on;
 %%% plot trend %%
         x=(1:length(SSH_yr));
         idx = find(isnan(SSH_yr)==1)
         p1 = polyfit(x', SSH_yr,1);
         y1 = polyval(p1,x');
plot(1:length(SSH_yr),y1,'r','linew',2);

        title('TG - detrend');
        set(gca,'linew',2);
        xlabel('Time (date)');
        ylabel('elevation (m)');
        set(gca,'xtick',1:240:length(time));
        set(gca,'xticklabel',tt(1:240:length(time),:));
%         set(gca,'XTickLabelRotation',45);
        set(gca,'fontsize',20);
