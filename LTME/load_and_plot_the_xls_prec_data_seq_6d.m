%% load xls data to mat   
    close all; clear all; clc; 
f_path = 'D:\장기생태\Dynamic\06_river\data\유역평균강수량\';
for i = 1980:2018
[raw txt]=xlsread([f_path,'구례군(송정리)_수위표_유역평균강수량.xlsx'],num2str(i),'');

    mth1=raw(1:31, 2);  mth1(isnan(mth1)) = [];           %row=days column=month
    mth2=raw(1:31, 3);  mth2(isnan(mth2)) = [];
    mth3=raw(1:31, 4);  mth3(isnan(mth3)) = [];
    mth4=raw(1:31, 5);  mth4(isnan(mth4)) = [];
    mth5=raw(1:31, 6);  mth5(isnan(mth5)) = [];
    mth6=raw(1:31, 7);  mth6(isnan(mth6)) = [];
    mth7=raw(1:31, 8);  mth7(isnan(mth7)) = [];
    mth8=raw(1:31, 9);  mth8(isnan(mth8)) = [];
    mth9=raw(1:31, 10); mth9(isnan(mth9)) = [];
    mth10=raw(1:31, 11);mth10(isnan(mth10)) = [];
    mth11=raw(1:31, 12);mth11(isnan(mth11)) = [];
    mth12=raw(1:31, 13);mth12(isnan(mth12)) = [];
    
    mth1_sum=raw(32, 2);  mth1(isnan(mth1_sum)) = [];           %row=days column=month
    mth2_sum=raw(32, 3);  mth2(isnan(mth2_sum)) = [];
    mth3_sum=raw(32, 4);  mth3(isnan(mth3_sum)) = [];
    mth4_sum=raw(32, 5);  mth4(isnan(mth4_sum)) = [];
    mth5_sum=raw(32, 6);  mth5(isnan(mth5_sum)) = [];
    mth6_sum=raw(32, 7);  mth6(isnan(mth6_sum)) = [];
    mth7_sum=raw(32, 8);  mth7(isnan(mth7_sum)) = [];
    mth8_sum=raw(32, 9);  mth8(isnan(mth8_sum)) = [];
    mth9_sum=raw(32, 10); mth9(isnan(mth9_sum)) = [];
    mth10_sum=raw(32, 11);mth10(isnan(mth10_sum)) = [];
    mth11_sum=raw(32, 12);mth11(isnan(mth11_sum)) = [];
    mth12_sum=raw(32, 13);mth12(isnan(mth12_sum)) = [];


tt=i-1979;
merg_prec{tt} = [mth1; mth2; mth3; mth4; mth5; mth6; mth7; mth8; mth9; ...
    mth10; mth11; mth12;];
merg_prec_sum{tt} = [mth1_sum; mth2_sum; mth3_sum; mth4_sum; mth5_sum; mth6_sum; mth7_sum; mth8_sum; mth9_sum; ...
    mth10_sum; mth11_sum; mth12_sum;];
end
pre_merg_prec=merg_prec;
pre_merg_prec_sum=merg_prec_sum;

% save('songjung_prec_data_1980-2018.mat','merg_prec','merg_prec_sum');

% %% plot x,y scatter
% t_year = 1980:2018
% clearvars sj_prec sj_trans_out
% order_i = 0
% for i = t_year(1):t_year(end)
%     order_i = order_i+1
%     clearvars tempo_temp tempo_trans
% tempo_prec=merg_prec{t_year(order_i)-1979};
% tempo_trans = sj_dis.dis_pre_total{t_year(order_i)-1979}; %dis_pre_total
% if order_i == 1
%     sj_prec=tempo_prec;
%     sj_trans_out=tempo_trans;
% else
%     sj_prec=cat(1,sj_prec,tempo_prec);
%     sj_trans_out=cat(1,sj_trans_out,tempo_trans);
% end
% end
% 
% plot(sj_prec,sj_trans_out,'.')
% xlabel('rain (mm/m^2)');
% ylabel('discharge (m^3/s)'); grid on;
% 
% mu = mean(sj_prec);
% sig = std(sj_prec,0,1);
% test = (sj_prec - mu) ./ sig;
% 
% mu = mean(sj_trans_out);
% sig = std(sj_trans_out,0,1);
% test2 = (sj_trans_out - mu) ./ sig;
% 
% plot(test,test2,'.')
% xlabel('rain (mm/m^2)');
% ylabel('discharge (m^3/s)'); grid on;


%river discharge 1980~1989, calculated from water lv. & discharge relation
%from ykang
sj_dis=load('D:\장기생태\Dynamic\06_river\data\sj_1980to1996\songjung_discharge_1980to2018.mat');  % pre_merg_dis is discharge


t_year = 1980:2010
clearvars sj_prec sj_trans_out
order_i = 0
for i = t_year(1):t_year(end)
    order_i = order_i+1
    clearvars tempo_temp tempo_trans
tempo_prec=merg_prec{t_year(order_i)-1979};
tempo_trans = sj_dis.dis_pre_total{t_year(order_i)-1979}; %dis_pre_total
tempo_trans = movmean(tempo_trans,5);
if order_i == 1
    sj_prec=tempo_prec;
    sj_trans_out=tempo_trans;
else
    sj_prec=cat(1,sj_prec,tempo_prec);
    sj_trans_out=cat(1,sj_trans_out,tempo_trans);
end
end

clearvars sj_prec_c sj_trans_out_c
%% seq to 6 days
for i = 1:round(length(sj_prec)/6)
    if i == round(length(sj_prec)/6)
        % take all remains
        sj_prec_c{i,1} = sj_prec((i-1)*6+1:end)';
        sj_trans_out_c{i,1} = sj_trans_out((i-1)*6+1:end)';
    else
        sj_prec_c{i,1} = sj_prec((i-1)*6+1:i*6)';
        sj_trans_out_c{i,1} = sj_trans_out((i-1)*6+1:i*6)';
    end
end

t_year = 2011:2018
clearvars sj_prec_2 sj_trans_out_2
order_i = 0
for i = t_year(1):t_year(end)
    order_i = order_i+1
    clearvars tempo_temp tempo_trans
tempo_prec=merg_prec{t_year(order_i)-1979};
tempo_trans = sj_dis.dis_pre_total{t_year(order_i)-1979}; %dis_pre_total
tempo_trans = movmean(tempo_trans,5);
if order_i == 1
    sj_prec_2=tempo_prec;
    sj_trans_out_2=tempo_trans;
else
    sj_prec_2=cat(1,sj_prec_2,tempo_prec);
    sj_trans_out_2=cat(1,sj_trans_out_2,tempo_trans);
end
end

clearvars sj_prec_2_c sj_trans_out_2_c
%% seq to 6 days
for i = 1:round(length(sj_prec_2)/6)
    if i == round(length(sj_prec)/6)
        % take all remains
        sj_prec_2_c{i,1} = sj_prec_2((i-1)*6+1:end)';
        sj_trans_out_2_c{i,1} = sj_trans_out_2((i-1)*6+1:end)';
    else
        sj_prec_2_c{i,1} = sj_prec_2((i-1)*6+1:i*6)';
        sj_trans_out_2_c{i,1} = sj_trans_out_2((i-1)*6+1:i*6)';
    end
end


%% LSTM-sequence-to-sequence exam (from https://kr.mathworks.com/help/deeplearning/ug/sequence-to-sequence-regression-using-deep-learning.html )

% clearvars train_prec train_trans
% order_i = 0
% t_year = 1980:2010
% % t_year = 1999:2018
% for i =t_year(1):t_year(end)
%      order_i = order_i+1
% train_prec{order_i,1} =merg_prec{t_year(order_i)-1979}';
% train_trans{order_i,1} = movmean(sj_dis.dis_pre_total{t_year(order_i)-1979}',5); %dis_pre_total
% end

% clearvars test_prec test_trans
% t_year = 2011:2018
% % t_year = 1997:1998
% order_i = 0
% for i =t_year(1):t_year(end)
%      order_i = order_i+1
% test_prec{order_i,1} =merg_prec{t_year(order_i)-1979}';
% test_trans{order_i,1} = movmean(sj_dis.dis_pre_total{t_year(order_i)-1979}',5); %dis_pre_total
% end

% for i = 1:numel(train_prec)
%     if length(train_prec{i}) == 366
%         train_prec{i}(60) = [];
%         train_trans{i}(60) = [];
%     end
% end
% 
% for i = 1:numel(test_prec)
%     if length(test_prec{i}) == 366
%         test_prec{i}(60) = [];
%         test_trans{i}(60) = [];
%     end
% end

clearvars XTrain YTrain XTest YTest

% XTrain = train_prec;
% YTrain = train_trans;
% 
% XTest = test_prec;
% YTest = test_trans;

XTrain = sj_prec_c;
YTrain = sj_trans_out_c;

XTest = sj_prec_2_c;
YTest = sj_trans_out_2_c;

% numFeatures = size(XTrain{1},1)
clearvars net YPred
numFeatures = 1

%% 훈련 예측 변수 정규화하기
mu = mean([XTrain{:}],2);
sig = std([XTrain{:}],0,2);

% mu = mean(sj_prec);
% sig = std(sj_prec,0,1);

for i = 1:numel(XTrain)
    XTrain{i} = (XTrain{i} - mu) ./ sig;
end

for i=1:numel(XTrain)
    sequence = XTrain{i};
    sequenceLengths(i) = size(sequence,2);
%     sequenceLengths(i) = 2;
end

miniBatchSize = 120;
numfully_connected_layers = 365;
numResponses = size(YTrain{1},1);
% numResponses = 366;
numHiddenUnits = 128;
% numHiddenUnits = 128; %% best
% numHiddenUnits = 365;

layers = [ ...
    sequenceInputLayer(numFeatures)
    lstmLayer(numHiddenUnits,'OutputMode','sequence','StateActivationFunction','tanh','GateActivationFunction','sigmoid')
    lstmLayer(numHiddenUnits,'OutputMode','sequence','StateActivationFunction','tanh','GateActivationFunction','sigmoid')
    fullyConnectedLayer(numfully_connected_layers)
    dropoutLayer(0.2)
    fullyConnectedLayer(numResponses)
    reluLayer
    regressionLayer];

%%% OLD, best
% layers = [ ...
%     sequenceInputLayer(numFeatures)
%     lstmLayer(numHiddenUnits,'OutputMode','sequence','StateActivationFunction','tanh','GateActivationFunction','sigmoid')
%     lstmLayer(numHiddenUnits,'OutputMode','sequence','StateActivationFunction','tanh','GateActivationFunction','sigmoid')
%     fullyConnectedLayer(numfully_connected_layers)
%     dropoutLayer(0.5)
%     fullyConnectedLayer(numResponses)
%     reluLayer
%     regressionLayer];

% options = trainingOptions('adam', ...
%     'MaxEpochs',maxEpochs, ...
%     'MiniBatchSize',miniBatchSize, ...
%     'InitialLearnRate',0.01, ...
%     'GradientThreshold',1, ...
%     'Shuffle','never', ...
%     'Plots','training-progress',...
%     'Verbose',0);

maxEpochs = 200;

options = trainingOptions('adam', ...
    'MaxEpochs',maxEpochs, ...
    'MiniBatchSize',miniBatchSize, ...
    'InitialLearnRate',0.01, ...
    'GradientThreshold',1, ...
    'Shuffle','never', ...
    'Plots','training-progress',...
    'Verbose',0);
%% load previous trains
% clearvars net net1 net2
% net1=load('machine_learned_200.mat','net');
% net2=load('machine_learned_1200.mat','net');

% net = trainNetwork(XTrain,YTrain,net1.net.Layers,options); % keep going

net = trainNetwork(XTrain,YTrain,layers,options); % origin

% clearvars XTest YTest
% XTest = test_prec;
% YTest = test_trans;

%% Test the Network
for i = 1:numel(XTest)
%     XTest{i}(idxConstant,:) = [];
    XTest{i} = (XTest{i} - mu) ./ sig;
%     YTest{i}(YTest{i} > thr) = thr;
end

clearvars YPred
YPred = predict(net,XTest,'MiniBatchSize',1);
% YPred = predict(net1.net,XTest,'MiniBatchSize',1);
% YPred = predict(net2.net,XTest,'MiniBatchSize',1);

clearvars Y_plot
for i = 1:61
    if i ==1
        Y_plot=YPred{i};
    else;
        Y_plot=cat(2,Y_plot,YPred{i});
    end
end

clearvars Y_test_plot
for i = 1:61
    if i ==1
        Y_test_plot=YTest{i};
    else;
        Y_test_plot=cat(2,Y_test_plot,YTest{i});
    end
end

plot(Y_plot,'b'); hold on; plot(Y_test_plot,'r');

idx = randperm(numel(YPred),8);
% idx = randperm(numel(YPred),2);
figure
for i = 1:numel(idx)
    subplot(2,4,i)
% subplot(2,1,i)
    
    plot(YTest{idx(i)},'-','linew',2)
    hold on
    plot(YPred{idx(i)},'--','linew',2)
    grid on;
    hold off  
    ylim([0 4500])
    xlim([1 365])
    title("Test Observation " + (idx(i)+2010))
    xlabel("days")
    ylabel("m^3/s")
end
legend(["Test Data" "Predicted"],'Location','southeast')

for i = 1:numel(YTest)
    YTestLast(i) = YTest{i}(end);
    YPredLast(i) = YPred{i}(end);
end
figure
rmse = sqrt(mean((YPredLast - YTestLast).^2))
histogram(YPredLast - YTestLast)
title("RMSE = " + rmse)
ylabel("Frequency")
xlabel("Error")

%% nonlinear regression

t_year = 1980:2018
clearvars sj_prec_tot sj_trans_out_tot
order_i = 0
for i = t_year(1):t_year(end)
    order_i = order_i+1
    clearvars tempo_temp tempo_trans
tempo_prec=merg_prec{t_year(order_i)-1979};
tempo_trans = sj_dis.dis_pre_total{t_year(order_i)-1979}; %dis_pre_total
if order_i == 1
    sj_prec_tot=tempo_prec;
    sj_trans_out_tot=tempo_trans;
else
    sj_prec_tot=cat(1,sj_prec_tot,tempo_prec);
    sj_trans_out_tot=cat(1,sj_trans_out_tot,tempo_trans);
end
end

t_year = 1980:2010
clearvars sj_prec sj_trans_out
order_i = 0
for i = t_year(1):t_year(end)
    order_i = order_i+1
    clearvars tempo_temp tempo_trans
tempo_prec=merg_prec{t_year(order_i)-1979};
tempo_trans = sj_dis.dis_pre_total{t_year(order_i)-1979}; %dis_pre_total
if order_i == 1
    sj_prec=tempo_prec;
    sj_trans_out=tempo_trans;
else
    sj_prec=cat(1,sj_prec,tempo_prec);
    sj_trans_out=cat(1,sj_trans_out,tempo_trans);
end
end


t_year = 2011:2018
clearvars sj_prec_2 sj_trans_out_2 tempo_prec tempo_trans
order_i = 0
for i = t_year(1):t_year(end)
    order_i = order_i+1
    clearvars tempo_temp tempo_trans
tempo_prec=merg_prec{t_year(order_i)-1979};
tempo_trans = sj_dis.dis_pre_total{t_year(order_i)-1979}; %dis_pre_total
if order_i == 1
    sj_prec_2=tempo_prec;
    sj_trans_out_2=tempo_trans;
else
    sj_prec_2=cat(1,sj_prec_2,tempo_prec);
    sj_trans_out_2=cat(1,sj_trans_out_2,tempo_trans);
end
end

mu = mean(sj_prec);
sig = std(sj_prec,0,1);
x = (sj_prec - mu) ./ sig;

% x=sj_prec;
y=sj_trans_out;
z = [ones(size(x)) x x.^2]; %[z]행렬을 생성함
a = polyfit(x,y,4); % polyfit 함수 % 사용
% a = z\y;

f = polyval(a,((sj_prec_2- mu) ./ sig));

clearvars Y_AI
for i = 1:length(YPred)
    if i== 1
       Y_AI =  YPred{i};
    else
        Y_AI =cat(2,Y_AI,YPred{i});
    end
end

sj_prec_tot_cms = (sj_prec_tot .* 1000) .* (61.89 * 10^6) ./ (3600*24); % km^2	* 10^6 = m^2, mm *1000 = m 

figure; 
plot(sj_prec_tot_cms,sj_trans_out_tot,'b.');
xlabel('CMS from rain ((A * rain)/86400)'); ylabel('OBS CMS'); grid on;
set(gca,'fontsize',13);

%% moving average input test
mv_trans_tot = movmean(sj_trans_out_tot,5);
figure; hold on; plot(sj_trans_out_tot,'r'); plot(mv_trans_tot,'b');

%%%%%%

figure; 
plot(sj_trans_out_2,f,'k.'); hold on;  
plot(sj_trans_out_2,Y_AI,'r.'); grid on; plot(sj_trans_out_2,sj_trans_out_2,'--','color',[0.5 0.5 0.5]);
xlim([0 500]); ylim([0 500])

%% RMSE compare
sqrt(mean((sj_trans_out_2-f).^2))
sqrt(mean((sj_trans_out_2-Y_AI').^2))

figure; plot((sj_trans_out_2-f),'k'); hold on; plot((sj_trans_out_2-Y_AI'),'r');


for i=1:numel(YTest)
    sequence = YTest{i};
    poly_sequenceLengths(i) = size(sequence,2);
end

for i = 1:numel(YTest)
yr_last(i)  =  sum(poly_sequenceLengths(1:i));
rmse(i) = sqrt(mean((YPred{i} - YTest{i}).^2));
yearly_YTest(i) = mean( YTest{i} );
end

clearvars poly_fy
for i = 1:numel(yr_last)
if i==1
  poly_fy{i,1} = f(1:yr_last(i))';
else
  poly_fy{i,1} = f(yr_last(i-1)+1:yr_last(i))';
end
end

for i = 1:numel(poly_fy)
rmse_poly(i) = sqrt(mean((poly_fy{i} - YTest{i}).^2));
end

figure; plot(rmse,'r'); hold on;  plot(rmse_poly,'k'); plot(yearly_YTest,'b');
xlabel('year'); ylabel('RMSE')
 grid on; set(gca,'fontsize',13);
% ylim([80 100]); yticks(80:2:100);
xticks(1:8); 
xticklabels(2011:2018)

clearvars Y_AI
for i = 1:length(YPred)
    if i== 1
       Y_AI =  YPred{i};
    else
        Y_AI =cat(2,Y_AI,YPred{i});
    end
end

figure; 
plot(f,'k--'); hold on;  plot(sj_trans_out_2,'b');
plot(Y_AI,'r--'); grid on; xlim([0 numel(Y_AI)]);

%% skill score come from https://www.jkscoe.or.kr/journal/view.php?number=81

for i = 1:numel(YTest)
ss(i) = (1 - rmse(i) ./ max(YTest{i})) .* 100;
ss_poly(i) = (1 - rmse_poly(i) ./ max(YTest{i})) .* 100;
end

skill_score = (1 - abs((Y_AI' - sj_trans_out_2)./(sj_trans_out_2))) .* 100;
skill_score_poly = (1 - abs((f - sj_trans_out_2)./(sj_trans_out_2))) .* 100;

plot(skill_score_poly,'k'); hold on; plot(skill_score,'r');
% ylim([-inf 100])
ylim([0 100])

figure;
plot(ss_poly,'k'); hold on; plot(ss,'r');
% ylim([-inf 100]);
xlabel('year'); ylabel('Skill Score ( % )')
ylim([70 100]); grid on; set(gca,'fontsize',13);
xticks(1:8); yticks(70:2:100);
xticklabels(2011:2018)


%% PRE & PTE come from https://www.readcube.com/library/a8d2c46c-1a54-49df-9418-96e0e68aea69:b83477eb-b5ec-4bd8-b55a-c10534e82a16

for i = 1:numel(YTest)
PRE(i) = ( abs(max(YPred{i}) - max(YTest{i})) / max(YTest{i}) ) .* 100;
PRE_poly(i) = ( abs(max(poly_fy{i}) - max(YTest{i})) / max(YTest{i}) ) .* 100;
PTE(i)= abs( find(YPred{i} == max(YPred{i})) - find(YTest{i} == max(YTest{i})) );
PTE_poly(i)= abs( find(poly_fy{i} == max(poly_fy{i})) - find(YTest{i} == max(YTest{i})) );
end

figure; hold on;
plot(PRE_poly,'k'); hold on; plot(PRE,'r');
ylabel('Peak discharge error ratio ( % )');
xlabel('year'); xticks(1:8); yticks(0:10:120);
xticklabels(2011:2018); set(gca,'fontsize',13);
ylim([0 120]); grid on;


figure; hold on;
plot(PTE_poly,'k'); hold on; plot(PTE,'r');
ylabel('Peak timing error ( days )');
xlabel('year'); xticks(1:8); yticks(0:5:100);
xticklabels(2011:2018); set(gca,'fontsize',13);
ylim([0 75]); grid on;


% b = a;
% b = [1;3;2];
% modelfun = @(b,x)(b(1)+b(2)*exp(-b(3)*x))
fun = @(x,xdata)x(1)*exp(x(2)*xdata);

opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';
beta0 = [2;2;2];
% beta = nlinfit(x,y,modelfun,beta0,opts)
beta = nlinfit(x',y',@hougen,beta0,opts)


% function f = fSSR(a, xm, ym)
% yp = a(1)*xm.^a(2);
% F = sum((ym-yp).^2);
% 
% fminsearch(@fSSR, [1, 1], [], sj_prec', sj_trans_out')

