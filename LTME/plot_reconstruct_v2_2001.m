close all; clear; clc;
load('TG_2021_G8_nondet.mat'); %pout_8(end-23:end)=[]; 
TG_tidestruc_8=tidestruc_8; TG_pout=pout_8;
load('model_2021_G8_nondet_t2.mat'); mod_tidestruc_8=tidestruc_8; mod_pout=pout_8;

figure()
plot(1:length(TG_pout), TG_pout./3.281,'linew',2); hold on;
% gap = 67+2189+12-4068+180-37;
gap = 0;
x = 1:length(mod_pout);
plot(x+gap, mod_pout./3.281,'color','r');
% line(time-datenum(2018,1,0),elevation1-pout,'linewi',2,'color','r');
xlabel('Time(Month)','fontsize',16,'fontweight','bold');
ylabel('Elevation (m)','fontsize',16,'fontweight','bold');
set(gca,'xtick',1:744:length(TG_pout));
set(gca,'xticklabel',[1:12]); grid('on');
legend({'TG-re','model-re'});
title('TG & Model reconstruction','fontsize',20,'fontweight','bold');  xlim([1 8761]); % ylim([-0.8 0.8]);
set(gca,'fontsize',16,'fontweight','bold');
xlim([1 length(TG_pout)])
%%%

figure;
plot(1:length(mod_pout), (TG_pout(1:length(mod_pout))-mod_pout)./3.281,'linew',2); hold on;
xlabel('Time(Month)','fontsize',16,'fontweight','bold');
ylabel('Elevation (m)','fontsize',16,'fontweight','bold');
% set(gca,'ytick',-0.8:0.2:0.8);
% set(gca,'yticklabel',-0.8:0.2:0.8);
% set(gca,'xtick',1:600:length(mod_t));
% set(gca,'xticklabel',tt_f(1:600:length(time_f),:));
% text(190,5.5,'Original Time series','color','b');
% text(190,4.75,'Tidal prediction from Analysis','color',[0 .5 0]);
% text(190,4.0,'Original time series minus Prediction','color','r');
% legend({'Original data','Tidal prediction', 'diff'})
xlabel('Time(Month)','fontsize',16,'fontweight','bold');
ylabel('Elevation (m)','fontsize',16,'fontweight','bold');
set(gca,'xtick',1:744:length(TG_pout));
set(gca,'xticklabel',[1:12]); grid('on');
legend({'diff','model-re'});
title('TG & Model reconstruction','fontsize',20,'fontweight','bold');  xlim([1 8761]); % ylim([-0.8 0.8]);
set(gca,'fontsize',16,'fontweight','bold');
xlim([1 length(TG_pout)])


%% 최대 창조 / 낙조
TG_part=TG_pout(1104:1272);

figure()
plot(1:length(TG_part), TG_part./3.281,'linew',2); hold on;
% gap = 67+2189+12-4068+180-37;
gap = 0;
x = 1:length(TG_part);
% plot(x+gap, mod_pout./3.281,'color','r');
% line(time-datenum(2018,1,0),elevation1-pout,'linewi',2,'color','r');
xlabel('2월 (일)','fontsize',16,'fontweight','bold');
ylabel('조위 (m)','fontsize',16,'fontweight','bold');
grid('on');
% legend({'TG-re','model-re'});
title('조위관측소 재구성 조위','fontsize',20,'fontweight','bold');  xlim([1 8761]); % ylim([-0.8 0.8]);
set(gca,'fontsize',16,'fontweight','bold');
xlim([1 length(TG_pout)])

% [31+15: 31+22].*24
xlim([1 length(TG_part)])
set(gca,'xtick',[1:24:length(TG_part).*24]);
set(gca,'xticklabel',[15:22]); 

zero_ind=find(TG_part == 0);


[pospks, poslocs] = findpeaks(TG_part);                                      % Positive Peaks & Locations (local maxima)
[negpks, neglocs] = findpeaks(-TG_part);                                     % Negative Peaks & Locations (local minima)
pkvct = [pospks negpks];                                                % Combine In One Vector FOr Sorting
[locs,idx] = sort([poslocs neglocs]);                                   % Sort
pks = pkvct(idx);
% locs_1d=reshape(locs,107*2,1);

figure; hold on;
plot(TG_part,'linew',2);
% plot(locs_1d,pout(locs_1d),'or');
plot(locs(:,1),pout(locs(:,1)),'or','linew',2);
plot(locs(:,2),pout(locs(:,2)),'oc','linew',2);
xlabel('Days in 2018','fontweight','bold');
ylabel('Elevation (m)','fontweight','bold');
% set(gca,'ytick',-0.8:0.2:0.8);
% set(gca,'yticklabel',-0.8:0.2:0.8);
set(gca,'xtick',1:60:length(tt));
set(gca,'xticklabel',tt(1:60:end,:));
legend({'Tidal reconstruction','Local maxima','Local minima'});grid('on');
ylim([-0.1 0.1]);
title('TD-pohang elevation');
set(gca,'fontsize',20,'fontweight','bold');
set(gca,'XTickLabelRotation',45);

plot(129*24,TG_pout(129*24)./3.281,'ro','linew',2)
plot(1123,TG_pout(1123)./3.281,'ro','linew',2)

all_indx=1:length(TG_part);

flood=[1:4, 11:17, 24:29, 35:41, 48:54, 60:66, 73:79, 84:90, 98:104, 109:115, 123:130, 134:140, 148:156, 161:166];

k=0;
for i = 1:length(all_indx)
if flood ~= all_indx(i)
    k=k+1;
    ebb(k) = i;
end
end





