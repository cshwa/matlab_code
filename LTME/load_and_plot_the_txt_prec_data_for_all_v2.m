%% load xls data to mat   
    close all; clear all; clc; 
    
     pick_st = [400403, 400404, 400501, 400502, 400503, 400504, 400601, 400602, ...
     400603, 400701, 400702, 400703, 400704, 400705, 400706, 400707, 400801, 400802, 400901, 400902, 400903, ...
     400904, 400101, 400102, 400103, 400104, 400105, 400106, 400107, 400108, 400109, 400201, 400202, 400203, ...
     400301, 400302, 400303, 400401, 400402];
 
 file_txt_path = 'D:\장기생태\Dynamic\06_river\data\유역평균강수량\표준유역 강수량\표준유역 강수량\';
 list_file_txt = dir([file_txt_path, 'climate_obs_sm_basin_map_data_*.txt']);
    
merg_prec = {}; mergmerg_prec_date_prec = {}; 
 for i = 1:length(list_file_txt)
      tempo_mat=[];
    tempo_mat=readmatrix([file_txt_path,list_file_txt(i).name],'Delimiter',',');      
     for j = 1:length(pick_st)
        tempo_idx=[]; tempo_mat_d=[];
        check_empty = 0;
        tempo_idx =find(tempo_mat(:,2) == pick_st(j));
        tempo_mat_d = tempo_mat(tempo_idx,3);
        if sum(isnan(tempo_mat_d)) > 0
            disp( [list_file_txt(i).name(end-7:end-4), ' point on ', num2str(j), ' code ', num2str(pick_st(j)), ' size ' num2str(sum(isnan(tempo_mat_d))) ]);
            if sum(isnan(tempo_mat_d)) < 20
                nandx = find(isnan(tempo_mat_d) == 1);
                nonandx = find(isnan(tempo_mat_d) == 0);
                tempo_mat_d(nandx) = interp1(nonandx, tempo_mat_d(nonandx), nandx);
                if sum(isnan(tempo_mat_d)) > 0
                    disp( [list_file_txt(i).name(end-7:end-4), ' point on ', num2str(j), ' code ', num2str(pick_st(j)), ' size ' num2str(sum(isnan(tempo_mat_d))) ,'!!! re-nan' ]);
                end
            end
            %% make it empty
            if sum(isnan(tempo_mat_d)) > 0
                check_empty = 1;
            end
        end
        if check_empty ~= 1
            merg_prec{i,j} = tempo_mat_d;
            merg_prec_date{i,j} = tempo_mat(tempo_idx,1);
            % warning for lack of data
            if length(tempo_mat_d) < 365
               disp([list_file_txt(i).name(end-7:end-4), ' point on ', num2str(j), ' code ', num2str(pick_st(j)), ' size ' num2str(sum(isnan(tempo_mat_d))) ,' has lack of length' ])
            end
        elseif check_empty == 1
            merg_prec{i,j} = { };
            merg_prec_date{i,j} = { };     
        end
     end
 end
 
 merg_prec(:,36) = [];
 merg_prec(:,26) = [];
 merg_prec(1:17,:) = [];
  merg_prec(1,:) = [];
 
 merg_prec_date(:,36) = [];
 merg_prec_date(:,26) = [];
 merg_prec_date(1:17,:) = [];
 merg_prec_date(1,:) = [];
 
 for i = 1:size(merg_prec,1)
     for j = 1:size(merg_prec,2)
     if size(merg_prec{i,j},1) < 365
        merg_prec(i,j)
        merg_prec_date{i,j}(1)
        i
        j
     end
     end
 end
 
 %% make 1989~present
k=0
for i = 1999:1999
            k=k+1;
    for j = 1:12
        eom_d_raw(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

%make t-axis for yearly 
temp_indx_pre=[];
for i = 1:size(eom_d_raw,1); temp_indx_pre = [temp_indx_pre eom_d_raw(i,:)]; end;

for i = 1:length(temp_indx_pre)
    if i ==1
         t_raw_indx(i) = temp_indx_pre(i);
    else
        t_raw_indx(i) = sum(temp_indx_pre(1:i));
    end
end

k=0; m=0;
for i = 1:length(1999:1999)
    l=0
    for n = 1:12
        m = m+1;
    for j = 1:eom_d_raw(i,n)
        l=l+1; % index for 1~31(or 30 or 29)
        k=k+1; % index for 1~366
        yymmdd_txt(k,:)=[num2str(i+1998) num2str(n,'%02d') num2str(j,'%02d')];
        yymmdd_txt_slash(k,:)=[num2str(i+1998) '/' num2str(n,'%02d') '/'  num2str(j,'%02d')];
    end
    end
end

tempo_new_prec = NaN(365,1);
yymmdd_num=str2num(yymmdd_txt);
tempo_prec=merg_prec{2,8};
for i = 1:length(yymmdd_num)
if length(find(yymmdd_num(i) == merg_prec_date{2,8})) > 0
   temp_idx = find(yymmdd_num(i) == merg_prec_date{2,8});
   tempo_new_prec(i) = tempo_prec(temp_idx);
end
end


nandx = find(isnan(tempo_new_prec) == 1);
nonandx = find(isnan(tempo_new_prec) == 0);
tempo_new_prec(nandx) = interp1(nonandx, tempo_new_prec(nonandx), nandx);

merg_prec{2,8} = tempo_new_prec;

%river discharge 1980~1989, calculated from water lv. & discharge relation
%from ykang
sj_dis=load('D:\장기생태\Dynamic\06_river\data\songjung_discharge_1980to2020.mat');  % pre_merg_dis is discharge

t_year = 1998:2018
clearvars sj_prec sj_trans_out
order_i = 0
for i = t_year(1):t_year(end)
    order_i = order_i+1
    clearvars tempo_temp tempo_trans
tempo_trans = sj_dis.dis_pre_total{t_year(order_i)-1979}; %dis_pre_total
if order_i == 1
    sj_trans_out=tempo_trans;
else
    sj_trans_out=cat(1,sj_trans_out,tempo_trans);
end
end


t_year = 2019:2020
clearvars sj_prec_2 sj_trans_out_2
order_i = 0
for i = t_year(1):t_year(end)
    order_i = order_i+1
    clearvars tempo_temp tempo_trans
tempo_trans = sj_dis.dis_pre_total{t_year(order_i)-1979}; %dis_pre_total
if order_i == 1
    sj_prec_2=tempo_prec;
else
    sj_prec_2=cat(1,sj_prec_2,tempo_prec);
end
end


%% LSTM-sequence-to-sequence exam (from https://kr.mathworks.com/help/deeplearning/ug/sequence-to-sequence-regression-using-deep-learning.html )

clearvars train_prec train_trans
order_i = 0
t_year = 1998:2012
% t_year = 1999:2018
for i =t_year(1):t_year(end)
     order_i = order_i+1
train_trans{order_i,1} = sj_dis.dis_pre_total{t_year(order_i)-1979}'; %dis_pre_total
end

clearvars test_prec test_trans
t_year = 2013:2020
% t_year = 1997:1998
order_i = 0
for i =t_year(1):t_year(end)
     order_i = order_i+1
test_trans{order_i,1} =  sj_dis.dis_pre_total{t_year(order_i)-1979}'; %dis_pre_total
end

clearvars XTrain YTrain XTest YTest

XTrain_pre = merg_prec(1:15,:);
YTrain = train_trans;

XTest_pre = merg_prec(16:23,:);
YTest = test_trans;


% numFeatures = size(XTrain{1},1)

%% make each cell row = feature (st. here), col = time
%% cell matrix's row = time
clearvars net YPred

for i = 1:size(XTrain_pre,1)
XTrain{i,1} = [XTrain_pre{i,:}]';
end

for i = 1:size(XTest_pre,1)
XTest{i,1} = [XTest_pre{i,:}]';
end

numFeatures = size(XTrain{1},1);

%% 훈련 예측 변수 정규화하기
mu = mean([XTrain{:}],2);
sig = std([XTrain{:}],0,2);

for i = 1:numel(XTrain)
    XTrain{i} = (XTrain{i} - mu) ./ sig;
end

for i=1:numel(XTrain)
    sequence = XTrain{i};
    sequenceLengths(i) = size(sequence,2);
%     sequenceLengths(i) = 2;
end

miniBatchSize = 3;
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

% every-epoch

maxEpochs = 2000;

options = trainingOptions('adam', ...
    'MaxEpochs',maxEpochs, ...
    'MiniBatchSize',miniBatchSize, ...
    'InitialLearnRate',0.01, ...
    'GradientThreshold',1, ...
    'Shuffle','every-epoch', ...
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

% idx = randperm(numel(YPred),8);
% idx = randperm(numel(YPred),2);
figure
for i = 1:numel(YPred)
    subplot(2,4,i)
% subplot(2,1,i)
    
    plot(YTest{i},'-','linew',2)
    hold on
    plot(YPred{i},'--','linew',2)
    grid on;
    hold off  
    ylim([0 4500])
    xlim([1 365])
    title("Test Observation " + (i+2012))
    xlabel("days")
    ylabel("m^3/s")
end
legend(["Test Data" "Predicted"],'Location','southeast')

clearvars YTestLast YPredLast
for i = 1:numel(YTest)
    if i == 1
    YTestLast = YTest{i};
    YPredLast = YPred{i};
    else
    YTestLast = cat(2,YTestLast, YTest{i});
    YPredLast = cat(2,YPredLast, YPred{i});
    end
end

figure
rmse = sqrt(mean((YPredLast - YTestLast).^2))
histogram(YPredLast - YTestLast)
title("RMSE = " + rmse)
ylabel("Frequency")
xlabel("Error")
xlim([-1000 1000])


figure; 
for i = 1:length(YPred)
rmse(i) = sqrt(mean((YPred{i} - YTest{i}).^2))
end
histogram(YPredLast - YTestLast)
title("RMSE = " + rmse)
ylabel("Frequency")
xlabel("Error")
xlim([-1000 1000])
return


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

