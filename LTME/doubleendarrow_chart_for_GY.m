close all; clear; clc;

addpath(genpath('/data1/cshwa/make_obc_ROMS_NWP/matlab/romsplot'));

% figure
% plot(1:10)
% x = [0.3 0.5];
% y = [0.6 0.5];
% annotation('textarrow',x,y,'String','y = x ')

figPos = [0 0 5 4];
fontSizeTick = 12;
fontSizeLab  = 12;
fontSizeCLab = 12;

figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos)
xlim([1980 2020])
ylim([0 6]); yticks([1,3,5]); yticklabels({'transport','nutrients','Chl.a'}); grid on; hold on;
plot(1998:2020,[repmat([4.5],length(1998:2006),1); repmat([5.5],length(2007:2015),1); repmat([4.5],length(2016:2020),1)],'g','linew',2);
plot(1989:2020,[repmat([3.5],length(1989:2004),1); repmat([2.5],length(2005:2020),1); ],'r','linew',2);
plot(1980:2020,[repmat([1.5],length(1980:1992),1); repmat([1.0],length(1993:2015),1); repmat([0.5],length(2016:2020),1); ],'b','linew',2);


% myaxes = axes();
sizeline = 6;
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos)
xlim([1982 2020])
ylim([0 6]); yticks([1,3,5]); grid on; hold on;
xticks(1982:5:2020); 
x = [1982:1986];
y = [repmat([5],length(x),1)];
plot(1980, 0 ,'.')
[normx, normy] = coord2norm(gca, [x(1) x(end)], [y(1) y(end)]);
annotation('doublearrow',normx,normy,'Color','r','LineWidth',sizeline,'Head1style','cback1','Head2style','cback1')
x2 = [1987:1992];
y2 = [repmat([4.5],length(x2),1)];
[normx2, normy2] = coord2norm(gca, [x2(1) x2(end)], [y2(1) y2(end)]);
annotation('doublearrow',normx2,normy2,'Color','r','LineWidth',sizeline,'Head1style','cback1','Head2style','cback1')
x3 = [1993:2004];
y3 = [repmat([3],length(x3),1)];
[normx3, normy3] = coord2norm(gca, [x3(1) x3(end)], [y3(1) y3(end)]);
annotation('doublearrow',normx3,normy3,'Color','b','LineWidth',sizeline,'Head1style','cback1','Head2style','cback1')
x4 = [2007:2015];
y4 = [repmat([1.5],length(x4),1)];
[normx4, normy4] = coord2norm(gca, [x4(1) x4(end)], [y4(1) y4(end)]);
annotation('doublearrow',normx4,normy4,'Color','k','LineWidth',sizeline,'Head1style','cback1','Head2style','cback1')
x5 = [2016:2020];
y5 = [repmat([1.0],length(x5),1)];
[normx5, normy5] = coord2norm(gca, [x5(1) x5(end)], [y5(1) y5(end)]);
annotation('doublearrow',normx5,normy5,'Color','k','LineWidth',sizeline,'Head1style','cback1','Head2style','cback1')

text(x(1),y(1)+0.4,'82~86','color','r')
text(x2(1),y2(1)+0.4,'87~92','color','r')
text(x3(1),y3(1)+0.4,'93~04','color','b')
text(x4(1),y4(1)+0.4,'07~15','color','k')
text(x5(1),y5(1)+0.4,'16~20','color','k')

text(x3(1),y(1),'Gwangyang-POSCO','color','r')
text(x4(2),y3(1)+0.5,'Juam-Dam','color','b')
text(x4(2),y3(1),'TN restriction','color','b')
text(x3(5),y5(1),'Chl.a (KOEM)','color','k')
set(gca,'YTick',[])
set(gca,'fontweight','bold');
box on;
print(gcf,['double_arrow_set_case_on_GY','.png'],'-dpng','-r200');


%% subplot 4,1

% case
sizeline = 6;
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos); hold on;
subplot(4,1,1);
xlim([1982 2020])
ylim([0 2]); yticks([1,3,5]); grid on; hold on;
xticks(1982:5:2020); 
x = [1982:0.1:1986.9];
y = [repmat([1],length(x),1)];
plot(1980, 0 ,'.')
[normx, normy] = coord2norm(gca, [x(1) x(end)], [y(1) y(end)]);
annotation('doublearrow',normx,normy,'Color','r','LineWidth',sizeline,'Head1style','cback1','Head2style','cback1')
x2 = [1987:0.1:1992.9];
y2 = [repmat([1],length(x2),1)];
[normx2, normy2] = coord2norm(gca, [x2(1) x2(end)], [y2(1) y2(end)]);
annotation('doublearrow',normx2,normy2,'Color','r','LineWidth',sizeline,'Head1style','cback1','Head2style','cback1')
x3 = [1993:0.1:2004.9];
y3 = [repmat([1],length(x3),1)];
[normx3, normy3] = coord2norm(gca, [x3(1) x3(end)], [y3(1) y3(end)]);
annotation('doublearrow',normx3,normy3,'Color','b','LineWidth',sizeline,'Head1style','cback1','Head2style','cback1')
x4 = [2007:0.1:2015.9];
y4 = [repmat([1],length(x4),1)];
[normx4, normy4] = coord2norm(gca, [x4(1) x4(end)], [y4(1) y4(end)]);
annotation('doublearrow',normx4,normy4,'Color','k','LineWidth',sizeline,'Head1style','cback1','Head2style','cback1')
x5 = [2016:2020];
y5 = [repmat([1],length(x5),1)];
[normx5, normy5] = coord2norm(gca, [x5(1) x5(end)], [y5(1) y5(end)]);
annotation('doublearrow',normx5,normy5,'Color','k','LineWidth',sizeline,'Head1style','cback1','Head2style','cback1')

text(x(1),y(1)+0.5,'82~86','color','k')
text(x2(10),y2(1)+0.5,'87~92','color','k')
text(x3(40),y3(1)+0.5,'93~04','color','k')
text(x4(30),y4(1)+0.5,'07~15','color','k')
text(x5(1),y5(1)+0.5,'16~20','color','k')
set(gca,'YTick',[])
set(gca,'fontweight','bold'); box on;
xline(1987,'r--','linew',2); % posco construction
xline(1993,'b--','linew',2); % posco construction
xline(2005,'g--','linew',2); % posco construction
xline(2016,'k--','linew',2); % posco construction

text(x2(1)+0.4,y2(1)-0.6,'POSCO','color','r')
text(x3(1)+0.4,y3(1)-0.6,'Juam-Dam','color','b')
text(x3(end)+0.5,y4(1)-0.6,'TN','color','g')
text(x4(end)+0.5,y5(1)-0.6,'Chl','color','k')

% transport
subplot(4,1,2); 
load('songjung_discharge_1980to2020.mat','dis_pre_total');
t_year = 1982:2020
clearvars tempo_trans sj_trans_out_365
order_i = 0
for i = t_year(1):t_year(end)
    order_i = order_i+1
    clearvars tempo_temp tempo_trans
% tempo_temp=sumjin_re_w_c{t_year(order_i)-1989};
tempo_trans = dis_pre_total{t_year(order_i)-1979}; %dis_pre_total
if order_i == 1
%     temperature=tempo_temp;
    sj_trans_out_yr=mean(tempo_trans);
else
    sj_trans_out_yr=cat(1,sj_trans_out_yr,mean(tempo_trans));
end
end
plot(1982:2020, sj_trans_out_yr,'b','linew',2);
xlim([1982 2020])
ylim([0 200]); yticks([0:50:200]); grid on; hold on;
xticks(1982:5:2020); 
ylabel('Discharge(CMS)')
set(gca,'fontweight','bold');

plot(1982:0.1:1992.9,repmat(mean(sj_trans_out_yr(1:11)),1,110),'r:','linew',3)
plot(1993:0.1:2003.9,repmat(mean(sj_trans_out_yr(12:22)),1,110),'g:','linew',3)
plot(2004:0.1:2014.9,repmat(mean(sj_trans_out_yr(23:33)),1,110),'k:','linew',3)
plot(2015:0.1:2020.9,repmat(mean(sj_trans_out_yr(34:end)),1,60),':','color',[.5 .5 .5],'linew',3)
% text
text(1987,150,num2str(round(mean(sj_trans_out_yr(1:11)),2),'%0.2f'),'color','r','fontweight','bold')
text(1996,150,num2str(round(mean(sj_trans_out_yr(12:22)),2),'%0.2f'),'color','g','fontweight','bold')
text(2008,150,num2str(round(mean(sj_trans_out_yr(23:33)),2),'%0.2f'),'color','k','fontweight','bold')
text(2015.2,150,num2str(round(mean(sj_trans_out_yr(34:end)),2),'%0.2f'),'color',[.5 .5 .5],'fontweight','bold')
xline(1993,'b--','linew',2); % Juam-dam construction



% nutrients
subplot(4,1,3); 
load('songjung_yymm_monthly_data_89to19_3sig.mat','monthly_*');
s20=load('songjung_yymm_monthly_data_to06to15_3sig_2020.mat','monthly_*');
monthly_no3(end+1:end+12)= s20.monthly_no3(end-11:end); %2020
monthly_nh4(end+1:end+12)= s20.monthly_nh4(end-11:end); %2020


for i = 1:length(monthly_no3)/12
yr_no3(i) = mean(monthly_no3((i-1)*12+1:i*12),'omitnan');
yr_nh4(i) = mean(monthly_nh4((i-1)*12+1:i*12),'omitnan');
end

clearvars plt_nh4 plt_no3 plt_po4
plt_nh4=yr_nh4 ./1000 .*14.006720.* sj_trans_out_yr(8:end)' .* 31536000 ./ 10^6 
plt_no3=yr_no3 ./1000 .*14.006720.* sj_trans_out_yr(8:end)' .* 31536000 ./ 10^6 
plt_DIN = plt_nh4 + plt_no3;

% plot(plt_no3,'r','linew',2);
% plot(plt_nh4,'k','linew',2);
plot(1989:2020,plt_DIN,'b','linew',2); ylabel('DIN(ton/yr)','color','k');
xlim([1982 2020])
ylim([0 7050]); yticks([0:2000:7050]); grid on; hold on;
xticks(1982:5:2020); set(gca,'fontweight','bold');

plot(1995:0.1:2004.9,repmat(mean(plt_DIN(7:16)),1,100),'g:','linew',3)
plot(2005:0.1:2014.9,repmat(mean(plt_DIN(17:26)),1,100),'k:','linew',3)
plot(2015:0.1:2020.9,repmat(mean(plt_DIN(27:end)),1,60),':','color',[.5 .5 .5],'linew',3)
% text
text(1997,1800,num2str(round(mean(plt_DIN(7:16)),2),'%0.2f'),'color','g','fontweight','bold')
text(2008,6000,num2str(round(mean(plt_DIN(17:26)),2),'%0.2f'),'color','k','fontweight','bold')
text(2015,5000,num2str(round(mean(plt_DIN(27:end)),2),'%0.2f'),'color','[.5 .5 .5]','fontweight','bold')
xline(2005,'g--','linew',2); % Juam-dam construction


% chlorophyll
subplot(4,1,4); 
yoonja=load('yoonjakangs_koem_data_monthly_v2_16points_2020.mat');
chl_sur_clim = yoonja.chl_sur;

chl_sur_clim = mean(yoonja.chl_sur,[1],'omitnan');
% [res I]=sort([1,5,4,3,2,6,8,9,7]); 
for i = 1:length(chl_sur_clim)/12
chl_sur_clim_yr(i) = mean(chl_sur_clim((i-1)*12+1:i*12),'omitnan');
end

        
sp_gy = size(chl_sur_clim_yr,1);
 for j=1:sp_gy %% num of spatial point
     for sig=1:6 %% sigma
                clearvars data data_ref
                data_ref = chl_sur_clim(j,:);
                data = chl_sur_clim(j,:);
                data(data_ref > nanmean(data_ref) + sig*nanstd(data_ref)) =NaN;
                data(data_ref < nanmean(data_ref) - sig*nanstd(data_ref)) =NaN;
                mtx_chl_s(j,:,sig) = data;  % extracted data OUT [sp_gy, T, sig]          
    end
 end

 chl_sur_clim_sp = squeeze(mean(mtx_chl_s(:,:,3),[1],'omitnan'));
%  chl_sur_clim_sp = squeeze(mean(mtx_chl_s(:,:,2),[1],'omitnan'));
% [res I]=sort([1,5,4,3,2,6,8,9,7]); 
for i = 1:length(chl_sur_clim_sp)/12
chl_sur_clim_yr(i) = mean(chl_sur_clim_sp((i-1)*12+1:i*12),'omitnan');
end


% plot(1997:2020,chl_sur_clim_yr,'b','linew',2);
plot(2007:2020,chl_sur_clim_yr(11:end),'b','linew',2);
xlim([1982 2020])
ylim([0 7]); yticks([0:2:7]); grid on; hold on;
xticks(1982:5:2020); 
ylabel('Chl.a(ug/L)')
set(gca,'fontweight','bold');

% plot(1999:0.1:2006.9,repmat(mean(chl_sur_clim_yr(1:10),'omitnan'),1,80),'g--','linew',3)
plot(2007:0.1:2015.9,repmat(mean(chl_sur_clim_yr(11:19),'omitnan'),1,90),'k:','linew',3)
plot(2016:0.1:2020.9,repmat(mean(chl_sur_clim_yr(20:24),'omitnan'),1,50),':','color',[.5 .5 .5],'linew',3)
% text
% text(2001,5,num2str(round(mean(chl_sur_clim_yr(1:10),'omitnan'),2),'%0.2f'),'color','g','fontweight','bold')
text(2009,2,num2str(round(mean(chl_sur_clim_yr(11:19),'omitnan'),2),'%0.2f'),'color','k','fontweight','bold')
text(2016.8,5,num2str(round(mean(chl_sur_clim_yr(20:24),'omitnan'),2),'%0.2f'),'color',[.5 .5 .5],'fontweight','bold')
xline(2016,'k--','linew',2); % posco construction
print(gcf,['double_arrow_set_case_on_GY_and_trans_nut_chl','.png'],'-dpng','-r200');






%% v2 
sizeline = 6;
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos)
xlim([1982 2020])
ylim([0 6]); yticks([0.5,1,2,2.5,3.5,4,5,5.5]); grid on; hold on;
xticks(1982:5:2020); 
plot(1980, 0 ,'.')
yticklabels({'low','high','after','before','after','before','after','before'});

x = [1982:1986];
y = [repmat([5.5],length(x),1)];
[normx, normy] = coord2norm(gca, [x(1) x(end)], [y(1) y(end)]);
annotation('doublearrow',normx,normy,'Color','r','LineWidth',sizeline,'Head1style','cback1','Head2style','cback1');
x1 = [1987:2020];
y1 = [repmat([5],length(x),1)];
[normx1, normy1] = coord2norm(gca, [x1(1) x1(end)], [y1(1) y1(end)]);
annotation('doublearrow',normx1,normy1,'Color','r','LineWidth',sizeline,'Head1style','cback1','Head2style','cback1');
x2 = [1982:1992];
y2 = [repmat([4],length(x2),1)];
[normx2, normy2] = coord2norm(gca, [x2(1) x2(end)], [y2(1) y2(end)]);
annotation('doublearrow',normx2,normy2,'Color','b','LineWidth',sizeline,'Head1style','cback1','Head2style','cback1');
x3 = [1993:2020];
y3 = [repmat([3.5],length(x3),1)];
[normx3, normy3] = coord2norm(gca, [x3(1) x3(end)], [y3(1) y3(end)]);
annotation('doublearrow',normx3,normy3,'Color','b','LineWidth',sizeline,'Head1style','cback1','Head2style','cback1');
x6 = [1982:2004];
y6 = [repmat([2.5],length(x6),1)];
[normx6, normy6] = coord2norm(gca, [x6(1) x6(end)], [y6(1) y6(end)]);
annotation('doublearrow',normx6,normy6,'Color','m','LineWidth',sizeline,'Head1style','cback1','Head2style','cback1');
x7 = [2005:2020];
y7 = [repmat([2],length(x7),1)];
[normx7, normy7] = coord2norm(gca, [x7(1) x7(end)], [y7(1) y7(end)]);
annotation('doublearrow',normx7,normy7,'Color','m','LineWidth',sizeline,'Head1style','cback1','Head2style','cback1');
x4 = [2007:2015];
y4 = [repmat([1],length(x4),1)];
[normx4, normy4] = coord2norm(gca, [x4(1) x4(end)], [y4(1) y4(end)]);
annotation('doublearrow',normx4,normy4,'Color','k','LineWidth',sizeline,'Head1style','cback1','Head2style','cback1');
x5 = [2016:2020];
y5 = [repmat([0.5],length(x5),1)];
[normx5, normy5] = coord2norm(gca, [x5(1) x5(end)], [y5(1) y5(end)]);
annotation('doublearrow',normx5,normy5,'Color','k','LineWidth',sizeline,'Head1style','cback1','Head2style','cback1');

% text(x(1),y(1)+0.4,'82~86','color','r')
% text(x2(1),y2(1)+0.4,'82~92','color','r')
% text(x3(1),y3(1)+0.4,'93~04','color','b')
% text(x4(1),y4(1)+0.4,'07~15','color','k')
% text(x5(1),y5(1)+0.4,'16~20','color','k')

text(1998,y(1)-0.1,'Gwangyang-POSCO','color','r')
text(1998,y2(1)-0.1,'Juam-Dam','color','b')
text(1987,y7(1)+0.1,'TN restriction','color','m')
text(1995,y4(1)-0.1,'Chl.a (KOEM)','color','k')
set(gca,'fontweight','bold');
box on;
print(gcf,['double_arrow_longterm_change_on_GY','.png'],'-dpng','-r200');


