for sig=3:3  %sigma
clearvars b1 yCalc1 X b nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_recon yCalc3
% raw
w_no3_04=w_no3(1:tx_tick(25));
idx1 = find(isnan(w_no3_04) == 0);
w_no3_042=w_no3_04;
upper_bc(sig) = mean(w_no3_04(idx1)) + sig*std(w_no3_04(idx1));
lower_bc(sig) = mean(w_no3_04(idx1)) - sig*std(w_no3_04(idx1));
w_no3_042(find(w_no3_04 > mean(w_no3_04(idx1)) + sig*std(w_no3_04(idx1))))=NaN;
w_no3_042(find(w_no3_04 < mean(w_no3_04(idx1)) - sig*std(w_no3_04(idx1))))=NaN;
 
nanmean(w_no3_042(idx1))

%regression
%slope y = b1*x
nonan=w_no3_042(~isnan(w_no3_042));
idx = find(isnan(w_no3_042) == 0);
nanaxisx=find(isnan(w_no3_042) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
b0(sig) = b(1);  b1(sig) = b(2);
yCalc2{sig} = [1:length(w_no3_04)].*b1(sig) + b0(sig);
plot(yCalc2{sig},'-.','color','k','linew',2)
line([1:1:length(w_no3_04)], upper_bc(sig),'color', color_spec(sig),'linew',2)
line([1:1:length(w_no3_04)], lower_bc(sig),'color', color_spec(sig),'linew',2)
yCalc3 = yCalc2{sig};
plot([1:1:length(w_no3_04)],mean(w_no3_04(idx1)),'-.','color',color_spec(sig-1),'linew',2);

%     for i = 1:length(yCalc3)
%             w_no3_recon(i) = yCalc3(i)  + yp_w_no3(indx_366{i});
%     end
% plot([1:length(w_no3)], w_no3_recon,'-.','color',color_spec(sig),'linew',1)
end

sumjin
%% DO
9.9777          10.9930

%% CHL
3.7214          4.8574

%% NH4
0.2029          0.0608

%% NO3
1.4567          1.2198


namgang
%% DO
8.4265          9.5404

%% CHL
7.4269          6.3259

%% NH4
0.0345         0.0543

%% NO3
0.9861         0.7880
