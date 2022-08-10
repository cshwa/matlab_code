close all; clear; clc;   % -v3
% % station name tag
% % (for preventing the miss picking using tag from 1 col's name on excel) 
%coastal
j=[7;];
for i = 1:1
name_tag_7{i} = ['대한해협연안 ' num2str(j(i),'%02d')] 
end

% combining the tag and outter point excluding
name_tag = name_tag_7; 
size_tag = length(name_tag);

% %  when skip the outer point
% % for i = 1:length(name_tag_7)-18  % 연안 3, 6, 13:28 has to be remove (18)  [it's located at out of domain] 
% %  new_i = [1:28]; new_i(3)=[]; new_i(6)=[]; new_i(1)=[];
% %  name_tag{size_tag+i} = name_tag_7{new_i(i)};
% % end


%% pick the row on the excel which has same name with tag

[raw_co txt_co]=xlsread('해양환경측정망(연안해역환경측정망).xls','sheet','');
txt_matc_co = txt_co(3:end,1);
txt_date_co = txt_co(3:end,3);
temp_sur_co = txt_co(3:end,4); 
temp_bot_co = txt_co(3:end,5); 
salt_sur_co = txt_co(3:end,6); 
salt_bot_co = txt_co(3:end,7);
do_sur_co = txt_co(3:end,10); 
do_bot_co = txt_co(3:end,11);
nh4_sur_co = txt_co(3:end,14); 
nh4_bot_co = txt_co(3:end,15);
no3_sur_co = txt_co(3:end,18); 
no3_bot_co = txt_co(3:end,19);
chl_sur_co = txt_co(3:end,32); 
chl_bot_co = txt_co(3:end,33);

merge_txt = [txt_matc_co;]; % name list
merge_date = [txt_date_co;]; % date list
merge_temp_sur = [temp_sur_co;];
merge_temp_bot = [temp_bot_co;];
merge_salt_sur = [salt_sur_co;];
merge_salt_bot = [salt_bot_co;];
merge_do_sur = [do_sur_co;];
merge_do_bot = [do_bot_co;];
merge_nh4_sur = [nh4_sur_co;];
merge_nh4_bot = [nh4_bot_co;];
merge_no3_sur = [no3_sur_co;];
merge_no3_bot = [no3_bot_co;];
merge_chl_sur = [chl_sur_co;];
merge_chl_bot = [chl_bot_co;];

merge_data_txt = [txt_co(3:end,4:end);];


%merge
for i = 1:length(name_tag)
   if  sum(strcmp(name_tag{i}, merge_txt)) ~= 0
       indx{i} = find([strcmp(name_tag{i}, merge_txt)] == 1)     
   end
end

%% make date to be 'yymm' form
for i = 1:length(merge_date)
temp = char(merge_date{i});
if size(temp) ~= 7
    temp=temp(1,1:7);
end
merge_yymm{i,1} = temp;
end

%% make date to be 'mm' form
for i = 1:length(merge_date)
temp = char(merge_date{i});
temp = temp(1,6:7);
merge_mm{i,1} = temp;
end

%% make 1997 to 2018 'yymm' form
k=0
for i = 1997:2018
    for j = 1:12
        k=k+1;
        ref_date{k,1} = [num2str(i) '-' num2str(j,'%02d')];
    end
end

%% make 1997 to 2018 'yymm' form
for j = 1:12
 ref_date_mm{j,1} = [num2str(j,'%02d')];
end

% matched date 'yymm' form
for j = 1:length(indx) % st. axis
    for i = 1:length(ref_date) % date axis
       if  sum(strcmp(ref_date{i}, merge_yymm(indx{j}))) ~= 0
           indx_date{j,i} = find([strcmp(ref_date{i}, merge_yymm(indx{j}))] == 1);     
       end
    end
end

% matched date 'mm' form
for j = 1:length(indx) % st. axis
    for i = 1:length(ref_date_mm) % date axis
       if  sum(strcmp(ref_date_mm{i}, merge_mm(indx{j}))) ~= 0
           indx_date_mm{j,i} = find([strcmp(ref_date_mm{i}, merge_mm(indx{j}))] == 1);     
       end
    end
end


% make climate
%temp
clearvars temp
for i = 1:length(indx)
    temp = merge_temp_sur(indx{i});
    for j = 1:size(indx_date,2) %mth
        temp_sur_clim(i,j) = mean(str2num(char(temp(indx_date{i,j}))));
    end
end

clearvars temp
for i = 1:length(indx)
    temp = merge_temp_bot(indx{i});
    for j = 1:size(indx_date,2) %mth
        temp_bot_clim(i,j) = mean(str2num(char(temp(indx_date{i,j}))));
    end
end

%salt
clearvars temp
for i = 1:length(indx)
    temp = merge_salt_sur(indx{i});
    for j = 1:size(indx_date,2) %mth
        salt_sur_clim(i,j) = mean(str2num(char(temp(indx_date{i,j}))));
    end
end

clearvars temp
for i = 1:length(indx)
    temp = merge_salt_bot(indx{i});
    for j = 1:size(indx_date,2) %mth
        salt_bot_clim(i,j) = mean(str2num(char(temp(indx_date{i,j}))));
    end
end

%do
clearvars temp
for i = 1:length(indx)
    temp = merge_do_sur(indx{i});
    for j = 1:size(indx_date,2) %mth
        do_sur_clim(i,j) = mean(str2num(char(temp(indx_date{i,j}))));
    end
end

clearvars temp
for i = 1:length(indx)
    temp = merge_do_bot(indx{i});
    for j = 1:size(indx_date,2) %mth
        do_bot_clim(i,j) = mean(str2num(char(temp(indx_date{i,j}))));
    end
end

%nh4
clearvars temp
for i = 1:length(indx)
    temp = merge_nh4_sur(indx{i});
    for j = 1:size(indx_date,2) %mth
        nh4_sur_clim(i,j) = mean(str2num(char(temp(indx_date{i,j}))));
    end
end

clearvars temp
for i = 1:length(indx)
    temp = merge_nh4_bot(indx{i});
    for j = 1:size(indx_date,2) %mth
        nh4_bot_clim(i,j) = mean(str2num(char(temp(indx_date{i,j}))));
    end
end

%no3
clearvars temp
for i = 1:length(indx)
    temp = merge_no3_sur(indx{i});
    for j = 1:size(indx_date,2) %mth
        no3_sur_clim(i,j) = mean(str2num(char(temp(indx_date{i,j}))));
    end
end

clearvars temp
for i = 1:length(indx)
    temp = merge_no3_bot(indx{i});
    for j = 1:size(indx_date,2) %mth
        no3_bot_clim(i,j) = mean(str2num(char(temp(indx_date{i,j}))));
    end
end

%chl
clearvars temp
for i = 1:length(indx)
    temp = merge_chl_sur(indx{i});
    for j = 1:size(indx_date,2) %mth
        chl_sur_clim(i,j) = mean(str2num(char(temp(indx_date{i,j}))));
    end
end

clearvars temp
for i = 1:length(indx)
    temp = merge_chl_bot(indx{i});
    for j = 1:size(indx_date,2) %mth
        chl_bot_clim(i,j) = mean(str2num(char(temp(indx_date{i,j}))));
    end
end

return

t_tick=[1:12:length(1997:2018)*12]
tx_tick = t_tick;
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
gray = [128/255 128/255 128/255];
for i = 1:1
    clearvars temp_s temp_b
   temp_s= nh4_sur_clim(i,:);
   temp_b= nh4_bot_clim(i,:);
figure; hold on;
 plot(temp_s,'*','color','b','linew',2);
 plot(temp_b,'*','color','r','linew',2);   
   
t=1:length(1997:2018)*12;
%      temp_t_s(isnan(temp_t_s)) = interp1(t(~isnan(temp_t_s)), temp_t_s(~isnan(temp_t_s)), t(isnan(temp_t_s)) ); 
temp_s(isnan(temp_s)) = interp1( t(~isnan(temp_s)), temp_s(~isnan(temp_s)), t(isnan(temp_s)) ); 
temp_b(isnan(temp_b)) = interp1( t(~isnan(temp_b)), temp_b(~isnan(temp_b)), t(isnan(temp_b)) ); 

interp_nh4_b(:,i) = temp_b;
interp_nh4_s(:,i) = temp_s;

 plot(temp_s,'color','b');
 plot(temp_b,'color','r');  
 ylim([0 160])
xlabel('time(year)','fontsize',13)
ylabel('nh4 (ug/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:3:end)]);
set(gca,'xlim',[tx_tick(1) tx_tick(end)]);
set(gca,'xticklabel',1997:3:2018);
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',tx_tick(1:1:end));
title([name_tag{i}])
% legend('surface','bottom')
% print('-dpng',[name_tag{i},' nh4 st.', num2str(i),'.png'])
end
% save('KODC_st7_NH4.mat','nh4_sur_clim','nh4_bot_clim','t_tick');


sig = 3; %sigma
t_tick=[1:12:length(1997:2018)*12]
tx_tick = t_tick;
clearvars b1* X* b* yCalc* nonan* idx* idx1* *_bc_*
% color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
% gray = [128/255 128/255 128/255];
i = 1
    clearvars temp_ss temp_bb
   temp_ss= nh4_sur_clim(i,:);
   temp_bb= nh4_bot_clim(i,:);
figure; hold on;
%  plot(temp_ss,'*','color','b','linew',2);
%  plot(temp_bb,'*','color','r','linew',2);   
 
 t=1:length(1997:2018)*12;
%      temp_t_s(isnan(temp_t_s)) = interp1(t(~isnan(temp_t_s)), temp_t_s(~isnan(temp_t_s)), t(isnan(temp_t_s)) ); 
temp_ss(isnan(temp_ss)) = interp1( t(~isnan(temp_ss)), temp_ss(~isnan(temp_ss)), t(isnan(temp_ss)) ); 
temp_bb(isnan(temp_bb)) = interp1( t(~isnan(temp_bb)), temp_bb(~isnan(temp_bb)), t(isnan(temp_bb)) ); 

 
t_regime = length(1997:2008)*12 + 5;
temp_ss_front = temp_ss(1:t_regime);
temp_ss_back = temp_ss(t_regime+1:end);
temp_bb_front = temp_bb(1:t_regime);
temp_bb_back = temp_bb(t_regime+1:end);

idx_s_front = find(isnan(temp_ss_front) == 0);
idx_s_back = find(isnan(temp_ss_back) == 0);
idx_b_front = find(isnan(temp_bb_front) == 0);
idx_b_back = find(isnan(temp_bb_back) == 0);

upper_bc_s_front = mean(temp_ss_front(idx_s_front)) + sig*std(temp_ss_front(idx_s_front));
lower_bc_s_front = mean(temp_ss_front(idx_s_front)) - sig*std(temp_ss_front(idx_s_front));
upper_bc_b_front = mean(temp_bb_front(idx_s_front)) + sig*std(temp_bb_front(idx_b_front));
lower_bc_b_front = mean(temp_bb_front(idx_s_front)) - sig*std(temp_bb_front(idx_b_front));

upper_bc_s_back = mean(temp_ss_back(idx_s_back)) + sig*std(temp_ss_back(idx_s_back));
lower_bc_s_back = mean(temp_ss_back(idx_s_back)) - sig*std(temp_ss_back(idx_s_back));
upper_bc_b_back = mean(temp_bb_back(idx_s_back)) + sig*std(temp_bb_back(idx_b_back));
lower_bc_b_back = mean(temp_bb_back(idx_s_back)) - sig*std(temp_bb_back(idx_b_back));

temp_ss_front_r = temp_ss_front;
temp_ss_back_r = temp_ss_back;
temp_bb_front_r = temp_bb_front;
temp_bb_back_r = temp_bb_back;

temp_ss_front_r(find(temp_ss_front_r > mean(temp_ss_front_r(idx_s_front)) + sig*std(temp_ss_front_r(idx_s_front)))) = NaN;
temp_ss_back_r(find(temp_ss_back_r > mean(temp_ss_back_r(idx_s_back)) + sig*std(temp_ss_back_r(idx_s_back)))) = NaN; 
temp_bb_front_r(find(temp_bb_front_r > mean(temp_bb_front_r(idx_b_front)) + sig*std(temp_bb_front_r(idx_b_front)))) = NaN;
temp_bb_back_r(find(temp_bb_back_r > mean(temp_bb_back_r(idx_b_back)) + sig*std(temp_bb_back_r(idx_b_back)))) = NaN; 
temp_ss_front_r(find(temp_ss_front_r < mean(temp_ss_front_r(idx_s_front)) - sig*std(temp_ss_front_r(idx_s_front)))) = NaN;
temp_ss_back_r(find(temp_ss_back_r < mean(temp_ss_back_r(idx_s_back)) - sig*std(temp_ss_back_r(idx_s_back)))) = NaN; 
temp_bb_front_r(find(temp_bb_front_r < mean(temp_bb_front_r(idx_b_front)) - sig*std(temp_bb_front_r(idx_b_front)))) = NaN;
temp_bb_back_r(find(temp_bb_back_r < mean(temp_bb_back_r(idx_b_back)) - sig*std(temp_bb_back_r(idx_b_back)))) = NaN; 

line([1:0.1:t_regime],upper_bc_s_front,'color','b','linew',2);
line([t_regime+1:0.1:t_regime+length(temp_ss_back)],upper_bc_s_back,'color','b','linew',2);
line([1:0.1:t_regime],upper_bc_b_front,'color','r','linew',2);
line([t_regime+1:0.1:t_regime+length(temp_bb_back)],upper_bc_b_back,'color','r','linew',2);
line([1:0.1:t_regime],lower_bc_s_front,'color','b','linew',2);
line([t_regime+1:0.1:t_regime+length(temp_ss_back)],lower_bc_s_back,'color','b','linew',2);
line([1:0.1:t_regime],lower_bc_b_front,'color','r','linew',2);
line([t_regime+1:0.1:t_regime+length(temp_bb_back)],lower_bc_b_back,'color','r','linew',2);
   
%% reigme mean
 plot([1:.1:t_regime],nanmean(temp_ss_front_r),'color','c','linew',4);
 plot([t_regime+1:.1:t_regime+length(temp_ss_back)],nanmean(temp_ss_back_r),'color','c','linew',4);
 plot([1:.1:t_regime],nanmean(temp_bb_front_r),'color','m','linew',4);
 plot([t_regime+1:.1:t_regime+length(temp_bb_back)],nanmean(temp_bb_back_r),'color','m','linew',4); 
 
 plot([1:t_regime],temp_ss_front_r,'*','color','b','linew',2);
 plot([t_regime+1:t_regime+length(temp_ss_back)],temp_ss_back_r,'*','color','b','linew',2);
 plot([1:t_regime],temp_bb_front_r,'*','color','r','linew',2);
 plot([t_regime+1:t_regime+length(temp_bb_back)],temp_bb_back_r,'*','color','r','linew',2);
 
 ylim([0 160])
xlabel('time(year)','fontsize',13)
ylabel('nh4 (ug/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:3:end)]);
set(gca,'xlim',[tx_tick(1) tx_tick(end)]);
set(gca,'xticklabel',1997:3:2018);
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',tx_tick(1:1:end));
title([name_tag{i}])
% legend('surface','bottom')

%regression
%slope y = b1*x
nonan_s_front=temp_ss_front_r(~isnan(temp_ss_front_r));
nonan_s_back=temp_ss_back_r(~isnan(temp_ss_back_r));
nonan_b_front=temp_bb_front_r(~isnan(temp_bb_front_r));
nonan_b_back=temp_bb_back_r(~isnan(temp_bb_back_r));

idx1_s_front = find(isnan(temp_ss_front_r) == 0);
idx1_s_back = find(isnan(temp_ss_back_r) == 0);
idx1_b_front = find(isnan(temp_bb_front_r) == 0);
idx1_b_back = find(isnan(temp_bb_back_r) == 0);

%% surface
% Slope & Intercept y = b0 + b1*x
X_s_front = [ones(length(idx1_s_front'),1) idx1_s_front']; %b_0 b_1 
b_s_front = X_s_front\nonan_s_front';
b0_s_front = b_s_front(1);  b1_s_front = b_s_front(2);
yCalc_s_front = [1:length(temp_ss_front)].*b1_s_front + b0_s_front;

% Slope & Intercept y = b0 + b1*x
X_s_back = [ones(length(idx1_s_back'),1) idx1_s_back']; %b_0 b_1 
b_s_back = X_s_back\nonan_s_back';
b0_s_back = b_s_back(1);  b1_s_back = b_s_back(2);
yCalc_s_back = [length(temp_ss_front)+1:length(temp_ss_front)+length(temp_ss_back)].*b1_s_back + b0_s_back;

%% bottom
% Slope & Intercept y = b0 + b1*x
X_b_front = [ones(length(idx1_b_front'),1) idx1_b_front']; %b_0 b_1 
b_b_front = X_b_front\nonan_b_front';
b0_b_front = b_b_front(1);  b1_b_front = b_b_front(2);
yCalc_b_front = [1:length(temp_bb_front)].*b1_b_front + b0_b_front;

% Slope & Intercept y = b0 + b1*x
X_b_back = [ones(length(idx1_b_back'),1) idx1_b_back']; %b_0 b_1 
b_b_back = X_b_back\nonan_b_back';
b0_b_back = b_b_back(1);  b1_b_back = b_b_back(2);
yCalc_b_back = [length(temp_bb_front)+1:length(temp_bb_front)+length(temp_bb_back)].*b1_b_back + b0_b_back;

plot([1:length(temp_ss_front)], yCalc_s_front,'--' ,'color','b','linew',2)
plot([length(temp_ss_front)+1:length(temp_ss_front)+length(temp_ss_back)], yCalc_s_back,'--' ,'color','b','linew',2)
plot([1:length(temp_bb_front)], yCalc_b_front,'--' ,'color','r','linew',2)
plot([length(temp_bb_front)+1:length(temp_bb_front)+length(temp_bb_back)], yCalc_b_back,'--' ,'color','r','linew',2)

% print('-dpng',[name_tag{i},' regime nh4 st.', num2str(i),'.png'])

%% outlier extract (3sig)
%combine
temp_ss_3sig = [temp_ss_front_r'; temp_ss_back_r';];
temp_bb_3sig = [temp_bb_front_r'; temp_bb_back_r';];
temp_ss_3sig_int = temp_ss_3sig;
temp_bb_3sig_int = temp_bb_3sig;
%interp
t=1:length(1997:2018)*12;
%      temp_t_s(isnan(temp_t_s)) = interp1(t(~isnan(temp_t_s)), temp_t_s(~isnan(temp_t_s)), t(isnan(temp_t_s)) ); 
temp_ss_3sig_int(isnan(temp_ss_3sig_int)) = interp1( t(~isnan(temp_ss_3sig_int)), temp_ss_3sig_int(~isnan(temp_ss_3sig_int)), t(isnan(temp_ss_3sig_int)) ); 
temp_bb_3sig_int(isnan(temp_bb_3sig_int)) = interp1( t(~isnan(temp_bb_3sig_int)), temp_bb_3sig_int(~isnan(temp_bb_3sig_int)), t(isnan(temp_bb_3sig_int)) ); 

figure; hold on;
 plot([1:t_regime],temp_ss_front_r,'*','color','b','linew',2);
 plot([t_regime+1:t_regime+length(temp_ss_back)],temp_ss_back_r,'*','color','b','linew',2);
 plot([1:t_regime],temp_bb_front_r,'*','color','r','linew',2);
 plot([t_regime+1:t_regime+length(temp_bb_back)],temp_bb_back_r,'*','color','r','linew',2);
 
 plot([1:length(temp_bb_3sig_int)],temp_ss_3sig_int,'--','color','b','linew',2);
 plot([1:length(temp_bb_3sig_int)],temp_bb_3sig_int,'--','color','r','linew',2);
 
 % whole climate
 for i = 1:12
 nh4_mth_full(i) = nanmean(temp_ss_3sig_int(i:12:end));
 nh4_mth_full_c{i} = temp_ss_3sig_int(i:12:end);
 nh4_mth_full_c_b{i} = temp_bb_3sig_int(i:12:end);
 if i ~= 12  % 12mth is 21 sample
    nh4_mth_full_m(:,i) = temp_ss_3sig_int(i:12:end);
    nh4_mth_full_m_b(:,i) = temp_bb_3sig_int(i:12:end);
 elseif i ==12 % 12mth is 21 sample
    nh4_mth_full_m(:,i) = [temp_ss_3sig_int(i:12:end); NaN;];
    nh4_mth_full_m_b(:,i) = [temp_bb_3sig_int(i:12:end); NaN;];
 end
       
 nh4_mth_1st(i) = nanmean(temp_ss_3sig_int(i:12:t_regime));
 nh4_mth_1st_c{i} = temp_ss_3sig_int(i:12:t_regime);
 nh4_mth_1st_b_c{i} = temp_bb_3sig_int(i:12:t_regime);
 if i == 6 || 7 || 8 || 9 || 10 || 11 || 12
     nh4_mth_2nd(i) = nanmean(temp_ss_3sig_int(t_regime+7+i:12:end)) + temp_ss_3sig_int(144+i-1) /2;
     nh4_mth_2nd_c{i} = [temp_ss_3sig_int(t_regime+7+i:12:end); temp_ss_3sig_int(144+i);];
     nh4_mth_2nd_b_c{i} = [temp_bb_3sig_int(t_regime+7+i:12:end); temp_bb_3sig_int(144+i);];
 elseif i == 1 || 2 || 3 || 4 || 5
    nh4_mth_2nd(i) = nanmean(temp_ss_3sig_int(t_regime+7+i:12:end));
    nh4_mth_2nd_c{i} = temp_ss_3sig_int(t_regime+7+i:12:end);
    nh4_mth_2nd_b_c{i} = temp_bb_3sig_int(t_regime+7+i:12:end);
 end
 end
 
 % 1st regime climate
 for i =1:12
 if length(nh4_mth_1st_c{i}) == 13
     nh4_mth_1st_m(:,i) = [nh4_mth_1st_c{i};];
     nh4_mth_1st_b_m(:,i) = [nh4_mth_1st_b_c{i}];
 elseif length(nh4_mth_1st_c{i}) == 12
     nh4_mth_1st_m(:,i) = [nh4_mth_1st_c{i}; NaN;];
     nh4_mth_1st_b_m(:,i) = [nh4_mth_1st_b_c{i}; NaN;];
 end
 end
 
  % 2nd regime climate
 for i =1:12
 if length(nh4_mth_2nd_c{i}) == 10
     nh4_mth_2nd_m(:,i) = [nh4_mth_2nd_c{i};];
     nh4_mth_2nd_b_m(:,i) = [nh4_mth_2nd_b_c{i}];
 elseif length(nh4_mth_2nd_c{i}) == 9
     nh4_mth_2nd_m(:,i) = [nh4_mth_2nd_c{i}; NaN;];
     nh4_mth_2nd_b_m(:,i) = [nh4_mth_2nd_b_c{i}; NaN;];
 end
 end
  
 
figure; hold on;
 plot(nh4_mth_full,'--','color','b','linew',2);
 plot(nh4_mth_1st,'--','color','r','linew',2);
 plot(nh4_mth_2nd,'--','color','g','linew',2);
  for i = 1:12
     plot(i,nh4_mth_full_c{i},'*','color','b','linew',2);
 end
xlim([1 12])
xlabel('time(month)','fontsize',13)
ylabel('nh4 (ug/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[1:1:12]);
set(gca,'xticklabel',1:1:12);
legend('full','1st','2nd');


nh4_mm_s = [nh4_mth_full_m(:,1); nh4_mth_full_m(:,2); nh4_mth_full_m(:,3); ...
     nh4_mth_full_m(:,4); nh4_mth_full_m(:,5); nh4_mth_full_m(:,6); nh4_mth_full_m(:,7); ...
     nh4_mth_full_m(:,8); nh4_mth_full_m(:,9); nh4_mth_full_m(:,10); nh4_mth_full_m(:,11); nh4_mth_full_m(:,12);];
 
nh4_mm_b = [nh4_mth_full_m_b(:,1); nh4_mth_full_m_b(:,2); nh4_mth_full_m_b(:,3); ...
     nh4_mth_full_m_b(:,4); nh4_mth_full_m_b(:,5); nh4_mth_full_m_b(:,6); nh4_mth_full_m_b(:,7); ...
     nh4_mth_full_m_b(:,8); nh4_mth_full_m_b(:,9); nh4_mth_full_m_b(:,10); nh4_mth_full_m_b(:,11); nh4_mth_full_m_b(:,12);];

nh4_mm_s_1st = [nh4_mth_1st_m(:,1); nh4_mth_1st_m(:,2); nh4_mth_1st_m(:,3); ...
     nh4_mth_1st_m(:,4); nh4_mth_1st_m(:,5); nh4_mth_1st_m(:,6); nh4_mth_1st_m(:,7); ...
     nh4_mth_1st_m(:,8); nh4_mth_1st_m(:,9); nh4_mth_1st_m(:,10); nh4_mth_1st_m(:,11); nh4_mth_1st_m(:,12);];

 nh4_mm_b_1st = [nh4_mth_1st_b_m(:,1); nh4_mth_1st_b_m(:,2); nh4_mth_1st_b_m(:,3); ...
     nh4_mth_1st_b_m(:,4); nh4_mth_1st_b_m(:,5); nh4_mth_1st_b_m(:,6); nh4_mth_1st_b_m(:,7); ...
     nh4_mth_1st_b_m(:,8); nh4_mth_1st_b_m(:,9); nh4_mth_1st_b_m(:,10); nh4_mth_1st_b_m(:,11); nh4_mth_1st_b_m(:,12);];
 
nh4_mm_s_2nd = [nh4_mth_2nd_m(:,1); nh4_mth_2nd_m(:,2); nh4_mth_2nd_m(:,3); ...
     nh4_mth_2nd_m(:,4); nh4_mth_2nd_m(:,5); nh4_mth_2nd_m(:,6); nh4_mth_2nd_m(:,7); ...
     nh4_mth_2nd_m(:,8); nh4_mth_2nd_m(:,9); nh4_mth_2nd_m(:,10); nh4_mth_2nd_m(:,11); nh4_mth_2nd_m(:,12);];

 nh4_mm_b_2nd = [nh4_mth_2nd_b_m(:,1); nh4_mth_2nd_b_m(:,2); nh4_mth_2nd_b_m(:,3); ...
     nh4_mth_2nd_b_m(:,4); nh4_mth_2nd_b_m(:,5); nh4_mth_2nd_b_m(:,6); nh4_mth_2nd_b_m(:,7); ...
     nh4_mth_2nd_b_m(:,8); nh4_mth_2nd_b_m(:,9); nh4_mth_2nd_b_m(:,10); nh4_mth_2nd_b_m(:,11); nh4_mth_2nd_b_m(:,12);];
 
 
xx = repmat([1:12],22,1); 
xx_mm = [xx(:,1); xx(:,2); xx(:,3); ...
     xx(:,4); xx(:,5); xx(:,6); xx(:,7); ...
     xx(:,8); xx(:,9); xx(:,10); xx(:,11); xx(:,12);];

xx_1st = repmat([1:12],13,1); 
xx_mm_1st = [xx_1st(:,1); xx_1st(:,2); xx_1st(:,3); ...
     xx_1st(:,4); xx_1st(:,5); xx_1st(:,6); xx_1st(:,7); ...
     xx_1st(:,8); xx_1st(:,9); xx_1st(:,10); xx_1st(:,11); xx_1st(:,12);];
 xx_mm_1st_raw = [xx_1st(:,1); xx_1st(:,2); xx_1st(:,3); ...
     xx_1st(:,4); xx_1st(:,5); xx_1st(:,6); xx_1st(:,7); ...
     xx_1st(:,8); xx_1st(:,9); xx_1st(:,10); xx_1st(:,11); xx_1st(:,12);];
 
 xx_2nd = repmat([1:12],10,1); 
xx_mm_2nd = [xx_2nd(:,1); xx_2nd(:,2); xx_2nd(:,3); ...
     xx_2nd(:,4); xx_2nd(:,5); xx_2nd(:,6); xx_2nd(:,7); ...
     xx_2nd(:,8); xx_2nd(:,9); xx_2nd(:,10); xx_2nd(:,11); xx_2nd(:,12);];
 xx_mm_2nd_raw = [xx_2nd(:,1); xx_2nd(:,2); xx_2nd(:,3); ...
     xx_2nd(:,4); xx_2nd(:,5); xx_2nd(:,6); xx_2nd(:,7); ...
     xx_2nd(:,8); xx_2nd(:,9); xx_2nd(:,10); xx_2nd(:,11); xx_2nd(:,12);];
 
%  xx_1st = repmat([1:12],13,1); 
% xx_mm_1st = [xx_1st(:,1); xx_1st(:,2); xx_1st(:,3); ...
%      xx_1st(:,4); xx_1st(:,5); xx_1st(:,6); xx_1st(:,7); ...
%      xx_1st(:,8); xx_1st(:,9); xx_1st(:,10); xx_1st(:,11); xx_1st(:,12);];

% eliminate regime mean
nh4_mm_s_raw = nh4_mm_s;
nh4_mm_b_raw = nh4_mm_b;
nh4_mm_s_raw_1st = nh4_mm_s_1st;
nh4_mm_b_raw_1st = nh4_mm_b_1st;
nh4_mm_s_raw_2nd = nh4_mm_s_2nd;
nh4_mm_b_raw_2nd = nh4_mm_b_2nd;

nh4_mm_s =  nh4_mm_s - nanmean(nh4_mm_s_raw);
nh4_mm_b =  nh4_mm_b - nanmean(nh4_mm_b_raw);
nh4_mm_s_1st =  nh4_mm_s_1st - nanmean(nh4_mm_s_raw_1st);
nh4_mm_b_1st =  nh4_mm_b_1st - nanmean(nh4_mm_b_raw_1st);
nh4_mm_s_2nd =  nh4_mm_s_2nd - nanmean(nh4_mm_s_raw_2nd);
nh4_mm_b_2nd =  nh4_mm_b_2nd - nanmean(nh4_mm_b_raw_2nd);

nh4_mm_s_2nd_adv =  nh4_mm_s_2nd - (nanmean(nh4_mm_s_raw_2nd) - nanmean(nh4_mm_s_raw_1st)) - nanmean(nh4_mm_s_raw_1st);
nh4_mm_b_2nd_adv =  nh4_mm_b_2nd - (nanmean(nh4_mm_b_raw_2nd) - nanmean(nh4_mm_b_raw_1st)) - nanmean(nh4_mm_b_raw_1st);


%removing nan
nan_s = find(isnan(nh4_mm_s) ==1);
nan_b = find(isnan(nh4_mm_b) ==1);
nan_s_1st = find(isnan(nh4_mm_s_1st) ==1);
nan_b_1st = find(isnan(nh4_mm_b_1st) ==1);
nan_s_2nd = find(isnan(nh4_mm_s_2nd) ==1);
nan_b_2nd = find(isnan(nh4_mm_b_2nd) ==1);

nh4_mm_s(nan_s)=[];
nh4_mm_b(nan_b)=[];
xx_mm(nan_s)=[];

nh4_mm_s_1st(nan_s_1st)=[];
nh4_mm_b_1st(nan_b_1st)=[];
xx_mm_1st(nan_s_1st)=[];

nh4_mm_s_2nd(nan_s_2nd)=[];
nh4_mm_b_2nd(nan_b_2nd)=[];
xx_mm_2nd(nan_b_2nd)=[];


y_s_pre = polyfit(xx_mm, nh4_mm_s,3)
y_s = polyval(y_s_pre,1:12);

y_b_pre = polyfit(xx_mm, nh4_mm_b,3)
y_b = polyval(y_b_pre,1:12);

y_s_pre_1st = polyfit(xx_mm_1st, nh4_mm_s_1st,3)
y_s_1st = polyval(y_s_pre_1st,1:12);

y_b_pre_1st = polyfit(xx_mm_1st, nh4_mm_b_1st,3)
y_b_1st = polyval(y_b_pre_1st,1:12);

y_s_pre_2nd = polyfit(xx_mm_2nd, nh4_mm_s_2nd,3)
y_s_2nd = polyval(y_s_pre_2nd,1:12);

y_b_pre_2nd = polyfit(xx_mm_2nd, nh4_mm_b_2nd,3)
y_b_2nd = polyval(y_b_pre_2nd,1:12);

% plot(y_s) 
% plot(y_b,'r')

figure; hold on;
%  plot(nh4_mth_full,'--','color','b','linew',2);
%  plot(nh4_mth_1st,'--','color','r','linew',2);
%  plot(nh4_mth_2nd,'--','color','g','linew',2);
 plot(y_s + nanmean(nh4_mm_s_raw),'--','color','b','linew',2)
 plot(y_b + nanmean(nh4_mm_b_raw),'--','color','r','linew',2)
 plot(y_s_1st + nanmean(nh4_mm_s_raw_1st),'*-','color','b','linew',2)
 plot(y_b_1st + nanmean(nh4_mm_b_raw_1st),'*-','color','r','linew',2)
 plot(y_s_2nd + nanmean(nh4_mm_s_raw_2nd),'+-','color','b','linew',2)
 plot(y_b_2nd + nanmean(nh4_mm_b_raw_2nd),'+-','color','r','linew',2)  
  for i = 1:12
     plot(i,nh4_mth_full_c{i},'*','color','b','linew',2);
     plot(i,nh4_mth_full_c_b{i},'+','color','r','linew',2);
  end
xlim([1 12])
xlabel('time(month)','fontsize',13)
ylabel('nh4 (ug/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[1:1:12]);
set(gca,'xticklabel',1:1:12);
legend('full-sur','full-bot');


% plot(xx_mm_1st_raw,nh4_mm_s_raw_1st,'*','color','b','linew',2);
% plot(xx_mm_1st_raw,nh4_mm_b_raw_1st,'+','color','r','linew',2);
plot(xx_mm_2nd_raw,nh4_mm_s_raw_2nd,'*','color','b','linew',2);
plot(xx_mm_2nd_raw,nh4_mm_b_raw_2nd,'+','color','r','linew',2);


%%% compare regime polynomial

figure; hold on;
 plot(y_s + nanmean(nh4_mm_s_raw),'+-','color','b','linew',2)
%  plot(y_b + nanmean(nh4_mm_b_raw),'--','color','r','linew',2)
 plot(y_s_1st + nanmean(nh4_mm_s_raw_1st),'+-','color','r','linew',2)
%  plot(y_b_1st + nanmean(nh4_mm_b_raw_1st),'+','color','r','linew',2)
 plot(y_s_2nd + nanmean(nh4_mm_s_raw_2nd),'+-','color','g','linew',2)
%  plot(y_b_2nd + nanmean(nh4_mm_b_raw_2nd),'+-','color','r','linew',2)  
  
xlim([1 12])
xlabel('time(month)','fontsize',13)
ylabel('nh4 (ug/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[1:1:12]);
set(gca,'xticklabel',1:1:12);
legend('full-sur','1st-sur','2nd-sur');


figure; hold on;
%  plot(y_s + nanmean(nh4_mm_s_raw),'+-','color','b','linew',2)
 plot(y_b + nanmean(nh4_mm_b_raw),'+-','color','b','linew',2)
%  plot(y_s_1st + nanmean(nh4_mm_s_raw_1st),'+-','color','r','linew',2)
 plot(y_b_1st + nanmean(nh4_mm_b_raw_1st),'+-','color','r','linew',2)
%  plot(y_s_2nd + nanmean(nh4_mm_s_raw_2nd),'+-','color','g','linew',2)
 plot(y_b_2nd + nanmean(nh4_mm_b_raw_2nd),'+-','color','g','linew',2)  
  
xlim([1 12])
xlabel('time(month)','fontsize',13)
ylabel('nh4 (ug/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[1:1:12]);
set(gca,'xticklabel',1:1:12);
legend('full-bot','1st-bot','2nd-bot');


figure; hold on;
 plot(y_s + nanmean(nh4_mm_s_raw),'+-','color','b','linew',2)
 plot(y_b + nanmean(nh4_mm_b_raw),'<-','color','b','linew',2)
 plot(y_s_1st + nanmean(nh4_mm_s_raw_1st),'+-','color','r','linew',2)
 plot(y_b_1st + nanmean(nh4_mm_b_raw_1st),'<-','color','r','linew',2)
 plot(y_s_2nd + nanmean(nh4_mm_s_raw_2nd),'+-','color','g','linew',2)
 plot(y_b_2nd + nanmean(nh4_mm_b_raw_2nd),'<-','color','g','linew',2)  
  
xlim([1 12])
xlabel('time(month)','fontsize',13)
ylabel('nh4 (ug/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[1:1:12]);
set(gca,'xticklabel',1:1:12);
legend('surface','bottom');




