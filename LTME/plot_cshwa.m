clear all; clc; close all;

% model = load('transp_07_5&1.mat');
model1 = load('transp_08_8&3_fix.mat');
model2 = load('transp_08_8&3_dig.mat');
% model3 = load('transp_08_8&3_dig_4yr.mat');
model4 = load('transp_08_8&3_dig_4yr_1.mat');
model5 = load('test9_1yr.mat');
model6 = load('test9_2yr.mat');
obs = load('transport_Nishimura.dat');

% model_kor = model.transp(1:12,1);  
model1_kor = model1.transp(1:12,1);
model2_kor = model2.transp(1:12,1);
% model3_kor = model3.transp(1:12,1);
model4_kor = model4.transp(1:12,1);
model5_kor = model5.transp(1:12,1);
model6_kor = model6.transp(1:12,1);
obs_07 = obs(373:384,2);

% model_tsugaru = model.transp(1:12,2); 
% model_soya = model.transp(1:12,3);

model1_tsugaru = model1.transp(1:12,2);
model1_soya = model1.transp(1:12,3);

model2_tsugaru = model2.transp(1:12,2); 
model2_soya = model2.transp(1:12,3);

% model3_tsugaru = model3.transp(1:12,2); 
% model3_soya = model3.transp(1:12,3);

model4_tsugaru = model4.transp(1:12,2); 
model4_soya = model4.transp(1:12,3);

model5_tsugaru = model5.transp(1:12,2); 
model5_soya = model5.transp(1:12,3);

model6_tsugaru = model6.transp(1:12,2); 
model6_soya = model6.transp(1:12,3);

model5_combine = model5_tsugaru + model5_soya;

figure; hold on;
% plot(1:12, model_kor,'linewidth',4,'color','r');
plot(1:12, model1_kor,'linewidth',4,'color','g');
plot(1:12, obs_07,'linewidth',4,'color','b');
plot(1:12, model2_kor,'linewidth',4,'color','k');
% plot(1:12, model5_combine,'linewidth',4,'color','y');
plot(1:12, model4_kor,'linewidth',4,'color','c');
plot(1:12, model5_kor,'linewidth',4,'color','r');
plot(1:12, model6_kor,'linewidth',4,'color','m')
% plot(1:12, model5_combine,'linewidth',4,'color','m');
set(gca,'fontsize',17);
% legend('test7(5&2)','test8(8&3)','Nishimura','test8-dig','test8-dig-4yr','test8-dig-4yr_1');
legend('test8(8&3)','Nishimura','test8-dig','test8-dig-4yr_1','test9-1yr','test9-2yr');
% set(gcf,'Position',[200 100 1000 400])
set(gca, 'XTick',(1:1:12))
ylabel('Transport(Sv)','fontsize', 17);xlabel('Month','fontsize',17);
title('Korea strait transport','fontsize',17);
set(gca,'fontsize',17);

figure; hold on;
% plot(1:12, model_tsugaru,'linewidth',4,'color','r');
plot(1:12, model1_tsugaru,'linewidth',4,'color','g');
plot(1:12, model2_tsugaru,'linewidth',4,'color','k');
% plot(1:12, model3_tsugaru,'linewidth',4,'color','b');
plot(1:12, model4_tsugaru,'linewidth',4,'color','c');
plot(1:12, model5_tsugaru,'linewidth',4,'color','r');
plot(1:12, model6_tsugaru,'linewidth',4,'color','b')
set(gca,'fontsize',17);
% legend('test7(5&2)','test8(8&3)','test8-dig','test8-dig-4yr','test8-dig-4yr_1','test8_combine');
legend('test8(8&3)','test8-dig','test8-dig-4yr','test8-dig-4yr_1', 'test9-1yr', 'test9-2yr');
% set(gcf,'Position',[200 100 1000 400])
set(gca, 'XTick',(1:1:12))
ylabel('Transport(Sv)','fontsize', 17);xlabel('Month','fontsize',17);
title('Tsugaru strait transport','fontsize',17);
set(gca,'fontsize',17);

figure; hold on;
% plot(1:12, model_soya,'linewidth',4,'color','r');
plot(1:12, model1_soya,'linewidth',4,'color','g');
plot(1:12, model2_soya,'linewidth',4,'color','k');
% plot(1:12, model3_soya,'linewidth',4,'color','b');
plot(1:12, model4_soya,'linewidth',4,'color','c');
plot(1:12, model5_soya,'linewidth',4,'color','r');
 plot(1:12, model6_soya,'linewidth',4,'color','b');
set(gca,'fontsize',17);
% legend('test7(5&2)','test8(8&3)','test8-dig','test8-dig-4yr','test8-dig-4yr_1');
legend('test8(8&3)','test8-dig','test8-dig-4yr','test8-dig-4yr_1','test9-1yr','test9-2yr');
% set(gcf,'Position',[200 100 1000 400])
set(gca, 'XTick',(1:1:12))
ylabel('Transport(Sv)','fontsize', 17);xlabel('Month','fontsize',17);
title('Soya strait transport','fontsize',17);
set(gca,'fontsize',17);

mean(model4_kor)
mean(model4_tsugaru)
mean(model4_soya)

%% ´ë¸¸ÇØÇù
clear all; clc; close all;

model = load('transp_08_8&3_dig_4yr_taibei.mat');

model_taibei = model.transp(1:12,1);  

figure; hold on;
plot(1:12, model_taibei,'linewidth',4,'color','r');
set(gca,'fontsize',17);
legend('Taiwon');
% set(gcf,'Position',[200 100 1000 400])
set(gca, 'XTick',(1:1:12))
ylabel('Transport(Sv)','fontsize', 17);xlabel('Month','fontsize',17);
title('Taiwon strait transport','fontsize',17);
set(gca,'fontsize',17);

