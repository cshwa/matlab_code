close all; clear; clc;
load('TG_1980_G8_det_fix.mat'); %pout_8(end-23:end)=[]; 
TG_tidestruc_8=tidestruc_8; TG_pout=pout_8;
load('MODEL_1980_G8_det_fix.mat'); mod_tidestruc_8=tidestruc_8; mod_pout=pout_8;

figure()
plot(1:length(TG_pout), TG_pout./3.281,'linew',2); hold on;
% gap = 67+2189+12-4068+180-37;
gap = 0;
x = 1:length(mod_pout);
plot(x+gap, mod_pout./3.281,'color','r');
% line(time-datenum(2018,1,0),elevation1-pout,'linewi',2,'color','r');
xlabel('Time(Month)','fontsize',20,'fontweight','bold');
ylabel('Elevation (m)','fontsize',20,'fontweight','bold');
set(gca,'xtick',1:744:length(TG_pout));
set(gca,'xticklabel',[1:12]); grid('on');
legend({'TG-re','model-re'});
title('TG & Model reconstruction','fontsize',20,'fontweight','bold');  xlim([1 8761]); % ylim([-0.8 0.8]);
set(gca,'fontsize',20,'fontweight','bold');


    consti_name={'M2','S2','N2','K2','K1','O1','P1','Q1'};
    mod_em = mod_tidestruc_8.tidecon(:,1).*(30.48); mod_pha = mod_tidestruc_8.tidecon(:,3);
    TG_em = TG_tidestruc_8.tidecon(:,1).*(30.48); TG_pha = TG_tidestruc_8.tidecon(:,3);
    marker_spec = {'p','o','<','h','^','s','v','>'};
    
    
figure();hold on;
for i = 1:8
    plot(TG_em(i), mod_em(i),marker_spec{i},'linew',2,'MarkerEdgeColor',...
        [51/255 204/255 204/255],'MarkerFaceColor',[51/255 204/255 204/255],'MarkerSize',12); 
end
le=legend(consti_name); set(le,'fontsize',15,'location','southeast');
plot(0:110, 0:110,'linew',1,'color','b');
title('1980 Amplitude','fontsize',20,'fontweight','bold');
xlabel('Observation (cm)','fontsize',20,'fontweight','bold');
ylabel('Model (cm)','fontsize',20,'fontweight','bold');
set(gca,'fontsize',20,'fontweight','bold');
grid on; xlim([0 110]); ylim([0 110]);
print('-dpng',['1980em.png']); close;


figure();hold on;
for i = 1:8
    plot(TG_pha(i), mod_pha(i),marker_spec{i},'linew',2,'MarkerEdgeColor',...
        [51/255 204/255 204/255],'MarkerFaceColor',[51/255 204/255 204/255],'MarkerSize',12); 
end
le=legend(consti_name); set(le,'fontsize',15,'location','southeast');
plot(0:350, 0:350,'linew',1,'color','b');
title('1980 Phase','fontsize',20,'fontweight','bold');
xlabel('Observation (phase)','fontsize',20,'fontweight','bold');
ylabel('Model (phase)','fontsize',20,'fontweight','bold');
set(gca,'fontsize',20,'fontweight','bold');
grid on; xlim([0 350]); ylim([0 350]);
print('-dpng',['1980_pha.png']); close;


%% 1981
close all; clear; clc;
load('TG_1981_G8_det_fix.mat'); %pout_8(end-23:end)=[]; 
TG_tidestruc_8=tidestruc_8; TG_pout=pout_8;
load('MODEL_1981_G8_nondet_fix.mat'); mod_tidestruc_8=tidestruc_8; mod_pout=pout_8;

figure()
plot(1:length(TG_pout), TG_pout./3.281,'linew',2); hold on;
% gap = 67+2189+12-4068+180-37;
gap = 0;
x = 1:length(mod_pout);
plot(x+gap, mod_pout./3.281,'color','r');
% line(time-datenum(2018,1,0),elevation1-pout,'linewi',2,'color','r');
xlabel('Time(Month)','fontsize',20,'fontweight','bold');
ylabel('Elevation (m)','fontsize',20,'fontweight','bold');
set(gca,'xtick',1:744:length(TG_pout));
set(gca,'xticklabel',[1:12]); grid('on');
legend({'TG-re','model-re'});
title('TG & Model reconstruction','fontsize',20,'fontweight','bold');  xlim([1 8761]); % ylim([-0.8 0.8]);
set(gca,'fontsize',20,'fontweight','bold');


    consti_name={'M2','S2','N2','K2','K1','O1','P1','Q1'};
    mod_em = mod_tidestruc_8.tidecon(:,1).*(30.48); mod_pha = mod_tidestruc_8.tidecon(:,3);
    TG_em = TG_tidestruc_8.tidecon(:,1).*(30.48); TG_pha = TG_tidestruc_8.tidecon(:,3);
    marker_spec = {'p','o','<','h','^','s','v','>'};
    marker_col = [255/255 204/255 102/255];
    
figure();hold on;
for i = 1:8
    plot(TG_em(i), mod_em(i),marker_spec{i},'linew',2,'MarkerEdgeColor',...
        marker_col,'MarkerFaceColor',marker_col,'MarkerSize',12); 
end
le=legend(consti_name); set(le,'fontsize',15,'location','southeast');
plot(0:110, 0:110,'linew',1,'color','b');
title('1981 Amplitude','fontsize',20,'fontweight','bold');
xlabel('Observation (cm)','fontsize',20,'fontweight','bold');
ylabel('Model (cm)','fontsize',20,'fontweight','bold');
set(gca,'fontsize',20,'fontweight','bold');
grid on; xlim([0 110]); ylim([0 110]);
print('-dpng',['1981em.png']); close;

figure();hold on;
for i = 1:8
    plot(TG_pha(i), mod_pha(i),marker_spec{i},'linew',2,'MarkerEdgeColor',...
        marker_col,'MarkerFaceColor',marker_col,'MarkerSize',12); 
end
le=legend(consti_name); set(le,'fontsize',15,'location','southeast');
plot(0:350, 0:350,'linew',1,'color','b');
title('1981 Phase','fontsize',20,'fontweight','bold');
xlabel('Observation (phase)','fontsize',20,'fontweight','bold');
ylabel('Model (phase)','fontsize',20,'fontweight','bold');
set(gca,'fontsize',20,'fontweight','bold');
grid on; xlim([0 350]); ylim([0 350]);
print('-dpng',['1981_pha.png']);
close;

%% 1982
close all; clear; clc;

st_name = {'Gunsan', 'Heuksando', 'Yeosu', 'Tongyeong', 'Busan', ...
    'Jeju', 'Geomundo', 'Ulsan', 'Pohang', 'Mukho', 'Sokcho', 'Ulleungdo' };
for j = 1:length(st_name) 
clearvars -except st_name j

load(['TG_1982_G10_det_fix_',st_name{j},'.mat']); %pout_8(end-23:end)=[]; 
TG_tidestruc_10=tidestruc_10; TG_pout=pout_10;
load(['D:\장기생태\Dynamic\tide_plot\model_1982_6mins_G10_nondet_',st_name{j},'.mat']); mod_tidestruc_10=tidestruc_10; mod_pout=pout_10;

% figure()
% plot(1:length(TG_pout), TG_pout./3.281,'linew',2); hold on;
% % gap = 67+2189+12-4068+180-37;
% gap = 0;
% x = 1:0.1:length(TG_pout)+1;
% plot(x+gap, mod_pout./3.281,'color','r');
% % line(time-datenum(2018,1,0),elevation1-pout,'linewi',2,'color','r');
% xlabel('Time(Month)','fontsize',20,'fontweight','bold');
% ylabel('Elevation (m)','fontsize',20,'fontweight','bold');
% set(gca,'xtick',1:744:length(TG_pout));
% set(gca,'xticklabel',[1:12]); grid('on');
% legend({'TG-re','model-re'});
% title('TG & Model reconstruction','fontsize',20,'fontweight','bold');  xlim([1 8761]); % ylim([-0.8 0.8]);
% set(gca,'fontsize',20,'fontweight','bold');
% print('-dpng',['1982em_NWP_G10.png']); close;

    consti_name={'M2','S2','N2','K2','K1','O1','P1','Q1','MF','MM'};
    mod_em = mod_tidestruc_10.tidecon(:,1).*(30.48); mod_pha = mod_tidestruc_10.tidecon(:,3);
    TG_em = TG_tidestruc_10.tidecon(:,1).*(30.48); TG_pha = TG_tidestruc_10.tidecon(:,3);
    marker_spec = {'p','o','<','X','h','^','s','v','>','+'};
    marker_col =[255/255 153/255 204/255];
    
figure();hold on;
for i = 1:length(consti_name)
    plot(TG_em(i), mod_em(i),marker_spec{i},'linew',2,'MarkerEdgeColor',...
        marker_col,'MarkerFaceColor',marker_col,'MarkerSize',12); 
end
plot(0:110, 0:110,'linew',1,'color','b');
le=legend(consti_name); set(le,'fontsize',13,'location','northwest','NumColumns',3);
title(['1982 진폭 ',st_name{j}],'fontsize',20,'fontweight','bold');
xlabel('Observation (cm)','fontsize',20,'fontweight','bold');
ylabel('Model (cm)','fontsize',20,'fontweight','bold');
set(gca,'fontsize',20,'fontweight','bold');
grid on; xlim([0 110]); ylim([0 110]);
print('-dpng',['1982em_NWP_G10_',st_name{j},'.png']); close;

figure();hold on;
for i = 1:length(consti_name)
    plot(TG_pha(i), mod_pha(i),marker_spec{i},'linew',2,'MarkerEdgeColor',...
        marker_col,'MarkerFaceColor',marker_col,'MarkerSize',12); 
end
plot(0:350, 0:350,'linew',1,'color','b');
le=legend(consti_name); set(le,'fontsize',13,'location','northwest','NumColumns',3);
title(['1982 Phase ',st_name{j}],'fontsize',20,'fontweight','bold');
xlabel('Observation (phase)','fontsize',20,'fontweight','bold');
ylabel('Model (phase)','fontsize',20,'fontweight','bold');
set(gca,'fontsize',20,'fontweight','bold');
grid on; xlim([0 350]); ylim([0 350]);
print('-dpng',['1982_pha_NWP_G10_',st_name{j},'.png']); close;
end

close all; clear; clc;
load('TG_1982_G10_det_fix.mat'); %pout_8(end-23:end)=[]; 
TG_tidestruc_10=tidestruc_10; TG_pout=pout_10;
load('model_1982_6mins_G10_nondet_regrid_v3.mat'); mod_tidestruc_10=tidestruc_10; mod_pout=pout_10;

figure()
plot(1:length(TG_pout), TG_pout./3.281,'linew',2); hold on;
% gap = 67+2189+12-4068+180-37;
gap = 0;
x = 1:0.1:length(TG_pout)+1;
plot(x+gap, mod_pout./3.281,'color','r');
% line(time-datenum(2018,1,0),elevation1-pout,'linewi',2,'color','r');
xlabel('Time(Month)','fontsize',20,'fontweight','bold');
ylabel('Elevation (m)','fontsize',20,'fontweight','bold');
set(gca,'xtick',1:744:length(TG_pout));
set(gca,'xticklabel',[1:12]); grid('on');
legend({'TG-re','model-re'});
title('TG & Model reconstruction','fontsize',20,'fontweight','bold');  xlim([1 8761]); % ylim([-0.8 0.8]);
set(gca,'fontsize',20,'fontweight','bold');


    consti_name={'M2','S2','N2','K2','K1','O1','P1','Q1','MF','MM'};
    mod_em = mod_tidestruc_10.tidecon(:,1).*(30.48); mod_pha = mod_tidestruc_10.tidecon(:,3);
    TG_em = TG_tidestruc_10.tidecon(:,1).*(30.48); TG_pha = TG_tidestruc_10.tidecon(:,3);
    marker_spec = {'p','o','<','X','h','^','s','v','>','+'};
    marker_col =[255/255 153/255 204/255];
    
figure();hold on;
for i = 1:length(consti_name)
    plot(TG_em(i), mod_em(i),marker_spec{i},'linew',2,'MarkerEdgeColor',...
        marker_col,'MarkerFaceColor',marker_col,'MarkerSize',12); 
end
plot(0:110, 0:110,'linew',1,'color','b');
le=legend(consti_name); set(le,'fontsize',13,'location','northwest','NumColumns',3);
title('1982 Amplitude','fontsize',20,'fontweight','bold');
xlabel('Observation (cm)','fontsize',20,'fontweight','bold');
ylabel('Model (cm)','fontsize',20,'fontweight','bold');
set(gca,'fontsize',20,'fontweight','bold');
grid on; xlim([0 110]); ylim([0 110]);
print('-dpng',['1982em_NWP_G10.png']); close;

figure();hold on;
for i = 1:length(consti_name)
    plot(TG_pha(i), mod_pha(i),marker_spec{i},'linew',2,'MarkerEdgeColor',...
        marker_col,'MarkerFaceColor',marker_col,'MarkerSize',12); 
end
plot(0:350, 0:350,'linew',1,'color','b');
le=legend(consti_name); set(le,'fontsize',13,'location','northwest','NumColumns',3);
title('1982 Phase','fontsize',20,'fontweight','bold');
xlabel('Observation (phase)','fontsize',20,'fontweight','bold');
ylabel('Model (phase)','fontsize',20,'fontweight','bold');
set(gca,'fontsize',20,'fontweight','bold');
grid on; xlim([0 350]); ylim([0 350]);
print('-dpng',['1982_pha_NWP_G10.png']); close;

%% 1983
close all; clear; clc;
load('TG_1983_G8_det_fix.mat'); %pout_8(end-23:end)=[]; 
TG_tidestruc_8=tidestruc_8; TG_pout=pout_8;
load('MODEL_1983_G8_nondet_fix.mat'); mod_tidestruc_8=tidestruc_8; mod_pout=pout_8;

figure()
plot(1:length(TG_pout), TG_pout./3.281,'linew',2); hold on;
% gap = 67+2189+12-4068+180-37;
gap = 0;
x = 1:length(mod_pout);
plot(x+gap, mod_pout./3.281,'color','r');
% line(time-datenum(2018,1,0),elevation1-pout,'linewi',2,'color','r');
xlabel('Time(Month)','fontsize',20,'fontweight','bold');
ylabel('Elevation (m)','fontsize',20,'fontweight','bold');
set(gca,'xtick',1:744:length(TG_pout));
set(gca,'xticklabel',[1:12]); grid('on');
legend({'TG-re','model-re'});
title('TG & Model reconstruction','fontsize',20,'fontweight','bold');  xlim([1 8761]); % ylim([-0.8 0.8]);
set(gca,'fontsize',20,'fontweight','bold');


    consti_name={'M2','S2','N2','K2','K1','O1','P1','Q1'};
    mod_em = mod_tidestruc_8.tidecon(:,1).*(30.48); mod_pha = mod_tidestruc_8.tidecon(:,3);
    TG_em = TG_tidestruc_8.tidecon(:,1).*(30.48); TG_pha = TG_tidestruc_8.tidecon(:,3);
    marker_spec = {'p','o','<','h','^','s','v','>'};
    marker_col =[102/255 204/255 153/255];
    
figure();hold on;
for i = 1:8
    plot(TG_em(i), mod_em(i),marker_spec{i},'linew',2,'MarkerEdgeColor',...
        marker_col,'MarkerFaceColor',marker_col,'MarkerSize',12); 
end
le=legend(consti_name); set(le,'fontsize',15,'location','southeast');
plot(0:110, 0:110,'linew',1,'color','b');
title('1983 Amplitude','fontsize',20,'fontweight','bold');
xlabel('Observation (cm)','fontsize',20,'fontweight','bold');
ylabel('Model (cm)','fontsize',20,'fontweight','bold');
set(gca,'fontsize',20,'fontweight','bold');
grid on; xlim([0 110]); ylim([0 110]);
print('-dpng',['1983em.png']); close;

figure();hold on;
for i = 1:8
    plot(TG_pha(i), mod_pha(i),marker_spec{i},'linew',2,'MarkerEdgeColor',...
        marker_col,'MarkerFaceColor',marker_col,'MarkerSize',12); 
end
le=legend(consti_name); set(le,'fontsize',15,'location','southeast');
plot(0:350, 0:350,'linew',1,'color','b');
title('1983 Phase','fontsize',20,'fontweight','bold');
xlabel('Observation (phase)','fontsize',20,'fontweight','bold');
ylabel('Model (phase)','fontsize',20,'fontweight','bold');
set(gca,'fontsize',20,'fontweight','bold');
grid on; xlim([0 350]); ylim([0 350]);
print('-dpng',['1983_pha.png']); close;

%% 1984
close all; clear; clc;
load('TG_1984_G8_det_fix.mat'); %pout_8(end-23:end)=[]; 
TG_tidestruc_8=tidestruc_8; TG_pout=pout_8;
load('MODEL_1984_G8_nondet_fix.mat'); mod_tidestruc_8=tidestruc_8; mod_pout=pout_8;

figure()
plot(1:length(TG_pout), TG_pout./3.281,'linew',2); hold on;
% gap = 67+2189+12-4068+180-37;
gap = 0;
x = 1:length(mod_pout);
plot(x+gap, mod_pout./3.281,'color','r');
% line(time-datenum(2018,1,0),elevation1-pout,'linewi',2,'color','r');
xlabel('Time(Month)','fontsize',20,'fontweight','bold');
ylabel('Elevation (m)','fontsize',20,'fontweight','bold');
set(gca,'xtick',1:744:length(TG_pout));
set(gca,'xticklabel',[1:12]); grid('on');
legend({'TG-re','model-re'});
title('TG & Model reconstruction','fontsize',20,'fontweight','bold');  xlim([1 8761]); % ylim([-0.8 0.8]);
set(gca,'fontsize',20,'fontweight','bold');


    consti_name={'M2','S2','N2','K2','K1','O1','P1','Q1'};
    mod_em = mod_tidestruc_8.tidecon(:,1).*(30.48); mod_pha = mod_tidestruc_8.tidecon(:,3);
    TG_em = TG_tidestruc_8.tidecon(:,1).*(30.48); TG_pha = TG_tidestruc_8.tidecon(:,3);
    marker_spec = {'p','o','<','h','^','s','v','>'};
    marker_col =[51/255 153/255 51/255];
    
figure();hold on;
for i = 1:8
    plot(TG_em(i), mod_em(i),marker_spec{i},'linew',2,'MarkerEdgeColor',...
        marker_col,'MarkerFaceColor',marker_col,'MarkerSize',12); 
end
le=legend(consti_name); set(le,'fontsize',15,'location','southeast');
plot(0:110, 0:110,'linew',1,'color','b');
title('1984 Amplitude','fontsize',20,'fontweight','bold');
xlabel('Observation (cm)','fontsize',20,'fontweight','bold');
ylabel('Model (cm)','fontsize',20,'fontweight','bold');
set(gca,'fontsize',20,'fontweight','bold');
grid on; xlim([0 110]); ylim([0 110]);
print('-dpng',['1984em.png']); close;

figure();hold on;
for i = 1:8
    plot(TG_pha(i), mod_pha(i),marker_spec{i},'linew',2,'MarkerEdgeColor',...
        marker_col,'MarkerFaceColor',marker_col,'MarkerSize',12); 
end
le=legend(consti_name); set(le,'fontsize',15,'location','southeast');
plot(0:350, 0:350,'linew',1,'color','b');
title('1984 Phase','fontsize',20,'fontweight','bold');
xlabel('Observation (phase)','fontsize',20,'fontweight','bold');
ylabel('Model (phase)','fontsize',20,'fontweight','bold');
set(gca,'fontsize',20,'fontweight','bold');
grid on; xlim([0 350]); ylim([0 350]);
print('-dpng',['1984_pha.png']); close;



%% 1985
close all; clear; clc;
load('TG_1985_G8_det_fix.mat'); %pout_8(end-23:end)=[]; 
TG_tidestruc_8=tidestruc_8; TG_pout=pout_8;
load('MODEL_1985_G8_nondet_fix.mat'); mod_tidestruc_8=tidestruc_8; mod_pout=pout_8;

figure()
plot(1:length(TG_pout), TG_pout./3.281,'linew',2); hold on;
% gap = 67+2189+12-4068+180-37;
gap = 0;
x = 1:length(mod_pout);
plot(x+gap, mod_pout./3.281,'color','r');
% line(time-datenum(2018,1,0),elevation1-pout,'linewi',2,'color','r');
xlabel('Time(Month)','fontsize',20,'fontweight','bold');
ylabel('Elevation (m)','fontsize',20,'fontweight','bold');
set(gca,'xtick',1:744:length(TG_pout));
set(gca,'xticklabel',[1:12]); grid('on');
legend({'TG-re','model-re'});
title('TG & Model reconstruction','fontsize',20,'fontweight','bold');  xlim([1 8761]); % ylim([-0.8 0.8]);
set(gca,'fontsize',20,'fontweight','bold');


    consti_name={'M2','S2','N2','K2','K1','O1','P1','Q1'};
    mod_em = mod_tidestruc_8.tidecon(:,1).*(30.48); mod_pha = mod_tidestruc_8.tidecon(:,3);
    TG_em = TG_tidestruc_8.tidecon(:,1).*(30.48); TG_pha = TG_tidestruc_8.tidecon(:,3);
    marker_spec = {'p','o','<','h','^','s','v','>'};
    marker_col =[153/255 102/255 204/255];
    
figure();hold on;
for i = 1:8
    plot(TG_em(i), mod_em(i),marker_spec{i},'linew',2,'MarkerEdgeColor',...
        marker_col,'MarkerFaceColor',marker_col,'MarkerSize',12); 
end
le=legend(consti_name); set(le,'fontsize',15,'location','southeast');
plot(0:110, 0:110,'linew',1,'color','b');
title('1985 Amplitude','fontsize',20,'fontweight','bold');
xlabel('Observation (cm)','fontsize',20,'fontweight','bold');
ylabel('Model (cm)','fontsize',20,'fontweight','bold');
set(gca,'fontsize',20,'fontweight','bold');
grid on; xlim([0 110]); ylim([0 110]);
print('-dpng',['1985em.png']); close;

figure();hold on;
for i = 1:8
    plot(TG_pha(i), mod_pha(i),marker_spec{i},'linew',2,'MarkerEdgeColor',...
        marker_col,'MarkerFaceColor',marker_col,'MarkerSize',12); 
end
le=legend(consti_name); set(le,'fontsize',15,'location','southeast');
plot(0:350, 0:350,'linew',1,'color','b');
title('1985 Phase','fontsize',20,'fontweight','bold');
xlabel('Observation (phase)','fontsize',20,'fontweight','bold');
ylabel('Model (phase)','fontsize',20,'fontweight','bold');
set(gca,'fontsize',20,'fontweight','bold');
grid on; xlim([0 350]); ylim([0 350]);
print('-dpng',['1985_pha.png']); close;

%%
close all; clear; clc;
load('TG_1981_G8_det_fix.mat'); %pout_8(end-23:end)=[]; 
TG_tidestruc_8=tidestruc_8; TG_pout=pout_8;
load('MODEL_1981_G8_det_fix.mat'); mod_tidestruc_8=tidestruc_8; mod_pout=pout_8;

figure()
plot(1:length(TG_pout), TG_pout./3.281,'linew',2); hold on;
gap = 67+2189+12-4068+180-37;
% gap = 0;
x = 1:length(mod_pout);
plot(x+gap, mod_pout./3.281,'color','r');
% line(time-datenum(2018,1,0),elevation1-pout,'linewi',2,'color','r');
xlabel('Time(Month)','fontsize',20,'fontweight','bold');
ylabel('Elevation (m)','fontsize',20,'fontweight','bold');
set(gca,'xtick',1:744:length(TG_pout));
set(gca,'xticklabel',[1:12]); grid('on');
legend({'TG-re','model-re'});
title('TG & Model reconstruction','fontsize',20,'fontweight','bold');  xlim([1 8761]); % ylim([-0.8 0.8]);
set(gca,'fontsize',20,'fontweight','bold');

%%%

figure;
plot(1:length(mod_pout), (TG_pout(1:length(mod_pout))-mod_pout)./3.281,'linew',2); hold on;
xlabel('Time(Month)','fontsize',20,'fontweight','bold');
ylabel('Elevation (m)','fontsize',20,'fontweight','bold');
% set(gca,'ytick',-0.8:0.2:0.8);
% set(gca,'yticklabel',-0.8:0.2:0.8);
% set(gca,'xtick',1:600:length(mod_t));
% set(gca,'xticklabel',tt_f(1:600:length(time_f),:));
% text(190,5.5,'Original Time series','color','b');
% text(190,4.75,'Tidal prediction from Analysis','color',[0 .5 0]);
% text(190,4.0,'Original time series minus Prediction','color','r');
% legend({'Original data','Tidal prediction', 'diff'})
xlabel('Time(Month)','fontsize',20,'fontweight','bold');
ylabel('Elevation (m)','fontsize',20,'fontweight','bold');
set(gca,'xtick',1:744:length(TG_pout));
set(gca,'xticklabel',[1:12]); grid('on');
legend({'diff','model-re'});
title('TG & Model reconstruction','fontsize',20,'fontweight','bold');  xlim([1 8761]); % ylim([-0.8 0.8]);
set(gca,'fontsize',20,'fontweight','bold');
return
%% 6th order detrend
opol = 6;
[p,s,mu] = polyfit(TG_t',TG_pout,opol);
f_y = polyval(p,TG_t',[],mu);
plot(1:length(TG_t), (TG_pout- f_y)./3.281,'r'); hold on;

%% 1993
close all; clear; clc;
load('TG_1993_G8_nondet.mat'); %pout_8(end-23:end)=[]; 
TG_tidestruc_8=tidestruc_8; TG_pout=pout_8;
load('MODEL_1993_G8_nondet.mat'); mod_tidestruc_8=tidestruc_8; mod_pout=pout_8;

figure()
plot(1:length(TG_pout), TG_pout./3.281,'linew',2); hold on;
% gap = 67+2189+12-4068+180-37;
gap = 0;
x = 1:length(mod_pout);
plot(x+gap, mod_pout./3.281,'color','r');
% line(time-datenum(2018,1,0),elevation1-pout,'linewi',2,'color','r');
xlabel('Time(Month)','fontsize',20,'fontweight','bold');
ylabel('Elevation (m)','fontsize',20,'fontweight','bold');
set(gca,'xtick',1:744:length(TG_pout));
set(gca,'xticklabel',[1:12]); grid('on');
legend({'TG-re','model-re'});
title('TG & Model reconstruction','fontsize',20,'fontweight','bold');  xlim([1 8761]); % ylim([-0.8 0.8]);
set(gca,'fontsize',20,'fontweight','bold');


    consti_name={'M2','S2','N2','K2','K1','O1','P1','Q1'};
    mod_em = mod_tidestruc_8.tidecon(:,1).*(30.48); mod_pha = mod_tidestruc_8.tidecon(:,3);
    TG_em = TG_tidestruc_8.tidecon(:,1).*(30.48); TG_pha = TG_tidestruc_8.tidecon(:,3);
    marker_spec = {'p','o','<','h','^','s','v','>'};
    marker_col = [255/255 204/255 102/255];
    
figure();hold on;
for i = 1:8
    plot(TG_em(i), mod_em(i),marker_spec{i},'linew',2,'MarkerEdgeColor',...
        marker_col,'MarkerFaceColor',marker_col,'MarkerSize',12); 
end
le=legend(consti_name); set(le,'fontsize',15,'location','southeast');
plot(0:110, 0:110,'linew',1,'color','b');
title('1993 Amplitude','fontsize',20,'fontweight','bold');
xlabel('Observation (cm)','fontsize',20,'fontweight','bold');
ylabel('Model (cm)','fontsize',20,'fontweight','bold');
set(gca,'fontsize',20,'fontweight','bold');
grid on; xlim([0 110]); ylim([0 110]);
print('-dpng',['1993em.png']); close;

figure();hold on;
for i = 1:8
    plot(TG_pha(i), mod_pha(i),marker_spec{i},'linew',2,'MarkerEdgeColor',...
        marker_col,'MarkerFaceColor',marker_col,'MarkerSize',12); 
end
le=legend(consti_name); set(le,'fontsize',15,'location','southeast');
plot(0:350, 0:350,'linew',1,'color','b');
title('1993 Phase','fontsize',20,'fontweight','bold');
xlabel('Observation (phase)','fontsize',20,'fontweight','bold');
ylabel('Model (phase)','fontsize',20,'fontweight','bold');
set(gca,'fontsize',20,'fontweight','bold');
grid on; xlim([0 350]); ylim([0 350]);
print('-dpng',['1993_pha.png']);
close;

%% 1994
close all; clear; clc;
load('TG_1994_G8_nondet.mat'); %pout_8(end-23:end)=[]; 
TG_tidestruc_8=tidestruc_8; TG_pout=pout_8;
load('MODEL_1994_G8_nondet.mat'); mod_tidestruc_8=tidestruc_8; mod_pout=pout_8;

figure()
plot(1:length(TG_pout), TG_pout./3.281,'linew',2); hold on;
% gap = 67+2189+12-4068+180-37;
gap = 0;
x = 1:length(mod_pout);
plot(x+gap, mod_pout./3.281,'color','r');
% line(time-datenum(2018,1,0),elevation1-pout,'linewi',2,'color','r');
xlabel('Time(Month)','fontsize',20,'fontweight','bold');
ylabel('Elevation (m)','fontsize',20,'fontweight','bold');
set(gca,'xtick',1:744:length(TG_pout));
set(gca,'xticklabel',[1:12]); grid('on');
legend({'TG-re','model-re'});
title('TG & Model reconstruction','fontsize',20,'fontweight','bold');  xlim([1 8761]); % ylim([-0.8 0.8]);
set(gca,'fontsize',20,'fontweight','bold');


    consti_name={'M2','S2','N2','K2','K1','O1','P1','Q1'};
    mod_em = mod_tidestruc_8.tidecon(:,1).*(30.48); mod_pha = mod_tidestruc_8.tidecon(:,3);
    TG_em = TG_tidestruc_8.tidecon(:,1).*(30.48); TG_pha = TG_tidestruc_8.tidecon(:,3);
    marker_spec = {'p','o','<','h','^','s','v','>'};
    marker_col = [255/255 204/255 102/255];
    
figure();hold on;
for i = 1:8
    plot(TG_em(i), mod_em(i),marker_spec{i},'linew',2,'MarkerEdgeColor',...
        marker_col,'MarkerFaceColor',marker_col,'MarkerSize',12); 
end
le=legend(consti_name); set(le,'fontsize',15,'location','southeast');
plot(0:110, 0:110,'linew',1,'color','b');
title('1994 Amplitude','fontsize',20,'fontweight','bold');
xlabel('Observation (cm)','fontsize',20,'fontweight','bold');
ylabel('Model (cm)','fontsize',20,'fontweight','bold');
set(gca,'fontsize',20,'fontweight','bold');
grid on; xlim([0 110]); ylim([0 110]);
print('-dpng',['1994em.png']); close;

figure();hold on;
for i = 1:8
    plot(TG_pha(i), mod_pha(i),marker_spec{i},'linew',2,'MarkerEdgeColor',...
        marker_col,'MarkerFaceColor',marker_col,'MarkerSize',12); 
end
le=legend(consti_name); set(le,'fontsize',15,'location','southeast');
plot(0:350, 0:350,'linew',1,'color','b');
title('1994 Phase','fontsize',20,'fontweight','bold');
xlabel('Observation (phase)','fontsize',20,'fontweight','bold');
ylabel('Model (phase)','fontsize',20,'fontweight','bold');
set(gca,'fontsize',20,'fontweight','bold');
grid on; xlim([0 350]); ylim([0 350]);
print('-dpng',['1994_pha.png']);
close;

%% 1995
close all; clear; clc;
load('TG_1995_G8_nondet.mat'); %pout_8(end-23:end)=[]; 
TG_tidestruc_8=tidestruc_8; TG_pout=pout_8;
load('MODEL_1995_G8_nondet.mat'); mod_tidestruc_8=tidestruc_8; mod_pout=pout_8;

figure()
plot(1:length(TG_pout), TG_pout./3.281,'linew',2); hold on;
% gap = 67+2189+12-4068+180-37;
gap = 0;
x = 1:length(mod_pout);
plot(x+gap, mod_pout./3.281,'color','r');
% line(time-datenum(2018,1,0),elevation1-pout,'linewi',2,'color','r');
xlabel('Time(Month)','fontsize',20,'fontweight','bold');
ylabel('Elevation (m)','fontsize',20,'fontweight','bold');
set(gca,'xtick',1:744:length(TG_pout));
set(gca,'xticklabel',[1:12]); grid('on');
legend({'TG-re','model-re'});
title('TG & Model reconstruction','fontsize',20,'fontweight','bold');  xlim([1 8761]); % ylim([-0.8 0.8]);
set(gca,'fontsize',20,'fontweight','bold');


    consti_name={'M2','S2','N2','K2','K1','O1','P1','Q1'};
    mod_em = mod_tidestruc_8.tidecon(:,1).*(30.48); mod_pha = mod_tidestruc_8.tidecon(:,3);
    TG_em = TG_tidestruc_8.tidecon(:,1).*(30.48); TG_pha = TG_tidestruc_8.tidecon(:,3);
    marker_spec = {'p','o','<','h','^','s','v','>'};
    marker_col = [255/255 204/255 102/255];
    
figure();hold on;
for i = 1:8
    plot(TG_em(i), mod_em(i),marker_spec{i},'linew',2,'MarkerEdgeColor',...
        marker_col,'MarkerFaceColor',marker_col,'MarkerSize',12); 
end
le=legend(consti_name); set(le,'fontsize',15,'location','southeast');
plot(0:110, 0:110,'linew',1,'color','b');
title('1995 Amplitude','fontsize',20,'fontweight','bold');
xlabel('Observation (cm)','fontsize',20,'fontweight','bold');
ylabel('Model (cm)','fontsize',20,'fontweight','bold');
set(gca,'fontsize',20,'fontweight','bold');
grid on; xlim([0 110]); ylim([0 110]);
print('-dpng',['1995em.png']); close;

figure();hold on;
for i = 1:8
    plot(TG_pha(i), mod_pha(i),marker_spec{i},'linew',2,'MarkerEdgeColor',...
        marker_col,'MarkerFaceColor',marker_col,'MarkerSize',12); 
end
le=legend(consti_name); set(le,'fontsize',15,'location','southeast');
plot(0:350, 0:350,'linew',1,'color','b');
title('1995 Phase','fontsize',20,'fontweight','bold');
xlabel('Observation (phase)','fontsize',20,'fontweight','bold');
ylabel('Model (phase)','fontsize',20,'fontweight','bold');
set(gca,'fontsize',20,'fontweight','bold');
grid on; xlim([0 350]); ylim([0 350]);
print('-dpng',['1995_pha.png']);
close;

%% 1996
close all; clear; clc;
load('TG_1996_G8_nondet.mat'); %pout_8(end-23:end)=[]; 
TG_tidestruc_8=tidestruc_8; TG_pout=pout_8;
load('MODEL_1996_G8_nondet.mat'); mod_tidestruc_8=tidestruc_8; mod_pout=pout_8;

figure()
plot(1:length(TG_pout), TG_pout./3.281,'linew',2); hold on;
% gap = 67+2189+12-4068+180-37;
gap = 0;
x = 1:length(mod_pout);
plot(x+gap, mod_pout./3.281,'color','r');
% line(time-datenum(2018,1,0),elevation1-pout,'linewi',2,'color','r');
xlabel('Time(Month)','fontsize',20,'fontweight','bold');
ylabel('Elevation (m)','fontsize',20,'fontweight','bold');
set(gca,'xtick',1:744:length(TG_pout));
set(gca,'xticklabel',[1:12]); grid('on');
legend({'TG-re','model-re'});
title('TG & Model reconstruction','fontsize',20,'fontweight','bold');  xlim([1 8761]); % ylim([-0.8 0.8]);
set(gca,'fontsize',20,'fontweight','bold');


    consti_name={'M2','S2','N2','K2','K1','O1','P1','Q1'};
    mod_em = mod_tidestruc_8.tidecon(:,1).*(30.48); mod_pha = mod_tidestruc_8.tidecon(:,3);
    TG_em = TG_tidestruc_8.tidecon(:,1).*(30.48); TG_pha = TG_tidestruc_8.tidecon(:,3);
    marker_spec = {'p','o','<','h','^','s','v','>'};
    marker_col = [255/255 204/255 102/255];
    
figure();hold on;
for i = 1:8
    plot(TG_em(i), mod_em(i),marker_spec{i},'linew',2,'MarkerEdgeColor',...
        marker_col,'MarkerFaceColor',marker_col,'MarkerSize',12); 
end
le=legend(consti_name); set(le,'fontsize',15,'location','southeast');
plot(0:110, 0:110,'linew',1,'color','b');
title('1996 Amplitude','fontsize',20,'fontweight','bold');
xlabel('Observation (cm)','fontsize',20,'fontweight','bold');
ylabel('Model (cm)','fontsize',20,'fontweight','bold');
set(gca,'fontsize',20,'fontweight','bold');
grid on; xlim([0 110]); ylim([0 110]);
print('-dpng',['1996em.png']); close;

figure();hold on;
for i = 1:8
    plot(TG_pha(i), mod_pha(i),marker_spec{i},'linew',2,'MarkerEdgeColor',...
        marker_col,'MarkerFaceColor',marker_col,'MarkerSize',12); 
end
le=legend(consti_name); set(le,'fontsize',15,'location','southeast');
plot(0:350, 0:350,'linew',1,'color','b');
title('1996 Phase','fontsize',20,'fontweight','bold');
xlabel('Observation (phase)','fontsize',20,'fontweight','bold');
ylabel('Model (phase)','fontsize',20,'fontweight','bold');
set(gca,'fontsize',20,'fontweight','bold');
grid on; xlim([0 350]); ylim([0 350]);
print('-dpng',['1996_pha.png']);
close;

%% 1997
close all; clear; clc;
load('TG_1997_G8_nondet.mat'); %pout_8(end-23:end)=[]; 
TG_tidestruc_8=tidestruc_8; TG_pout=pout_8;
load('MODEL_1997_G8_nondet.mat'); mod_tidestruc_8=tidestruc_8; mod_pout=pout_8;

figure()
plot(1:length(TG_pout), TG_pout./3.281,'linew',2); hold on;
% gap = 67+2189+12-4068+180-37;
gap = 0;
x = 1:length(mod_pout);
plot(x+gap, mod_pout./3.281,'color','r');
% line(time-datenum(2018,1,0),elevation1-pout,'linewi',2,'color','r');
xlabel('Time(Month)','fontsize',20,'fontweight','bold');
ylabel('Elevation (m)','fontsize',20,'fontweight','bold');
set(gca,'xtick',1:744:length(TG_pout));
set(gca,'xticklabel',[1:12]); grid('on');
legend({'TG-re','model-re'});
title('TG & Model reconstruction','fontsize',20,'fontweight','bold');  xlim([1 8761]); % ylim([-0.8 0.8]);
set(gca,'fontsize',20,'fontweight','bold');


    consti_name={'M2','S2','N2','K2','K1','O1','P1','Q1'};
    mod_em = mod_tidestruc_8.tidecon(:,1).*(30.48); mod_pha = mod_tidestruc_8.tidecon(:,3);
    TG_em = TG_tidestruc_8.tidecon(:,1).*(30.48); TG_pha = TG_tidestruc_8.tidecon(:,3);
    marker_spec = {'p','o','<','h','^','s','v','>'};
    marker_col = [255/255 204/255 102/255];
    
figure();hold on;
for i = 1:8
    plot(TG_em(i), mod_em(i),marker_spec{i},'linew',2,'MarkerEdgeColor',...
        marker_col,'MarkerFaceColor',marker_col,'MarkerSize',12); 
end
le=legend(consti_name); set(le,'fontsize',15,'location','southeast');
plot(0:110, 0:110,'linew',1,'color','b');
title('1997 Amplitude','fontsize',20,'fontweight','bold');
xlabel('Observation (cm)','fontsize',20,'fontweight','bold');
ylabel('Model (cm)','fontsize',20,'fontweight','bold');
set(gca,'fontsize',20,'fontweight','bold');
grid on; xlim([0 110]); ylim([0 110]);
print('-dpng',['1997em.png']); close;

figure();hold on;
for i = 1:8
    plot(TG_pha(i), mod_pha(i),marker_spec{i},'linew',2,'MarkerEdgeColor',...
        marker_col,'MarkerFaceColor',marker_col,'MarkerSize',12); 
end
le=legend(consti_name); set(le,'fontsize',15,'location','southeast');
plot(0:350, 0:350,'linew',1,'color','b');
title('1997 Phase','fontsize',20,'fontweight','bold');
xlabel('Observation (phase)','fontsize',20,'fontweight','bold');
ylabel('Model (phase)','fontsize',20,'fontweight','bold');
set(gca,'fontsize',20,'fontweight','bold');
grid on; xlim([0 350]); ylim([0 350]);
print('-dpng',['1997_pha.png']);
close;
