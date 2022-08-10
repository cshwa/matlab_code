close all; clear; clc;   % -v3

% % station name tag
% % (for preventing the miss picking using tag from 1 col's name on excel) 

%coastal
j=[3; 5; 6; 7; 13; 14; 16; 17;];
for i = 1:8
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

%% pick matched name with tag
% %port
% for i = 1:length(name_tag)
%    if  sum(strcmp(name_tag{i}, txt_matc_p)) ~= 0
%        indx{i} = find([strcmp(name_tag{i}, txt_matc_p)] == 1)     
%    end
% end

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
for i = 1:8
    clearvars temp_s temp_b
   temp_s= do_sur_clim(i,:);
   temp_b= do_bot_clim(i,:);
figure; hold on;
 plot(temp_s,'*','color','b','linew',2);
 plot(temp_b,'*','color','r','linew',2);   
   
t=1:length(1997:2018)*12;
%      temp_t_s(isnan(temp_t_s)) = interp1(t(~isnan(temp_t_s)), temp_t_s(~isnan(temp_t_s)), t(isnan(temp_t_s)) ); 
temp_s(isnan(temp_s)) = interp1( t(~isnan(temp_s)), temp_s(~isnan(temp_s)), t(isnan(temp_s)) ); 
temp_b(isnan(temp_b)) = interp1( t(~isnan(temp_b)), temp_b(~isnan(temp_b)), t(isnan(temp_b)) ); 

interp_do_b(:,i) = temp_b;
interp_do_s(:,i) = temp_s;

 plot(temp_s,'color','b');
 plot(temp_b,'color','r');  
 ylim([2 16])
xlabel('time(year)','fontsize',13)
ylabel('do (mg/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:3:end)]);
set(gca,'xlim',[tx_tick(1) tx_tick(end)]);
set(gca,'xticklabel',1997:3:2018);
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',tx_tick(1:1:end));
title([name_tag{i}])
% legend('surface','bottom')
% print('-dpng',[name_tag{i},' DO st.', num2str(i),'.png'])
end

t_tick=[1:12:length(1997:2018)*12]
tx_tick = t_tick;
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
gray = [128/255 128/255 128/255];                                                                                                                                                                                                                                                                                                                                                                                                                                                       
for i = 1:8
    clearvars temp_s temp_b
   temp_s= chl_sur_clim(i,:);
   temp_b= chl_bot_clim(i,:);
figure; hold on;
 plot(temp_s,'*','color','b','linew',2);
 plot(temp_b,'*','color','r','linew',2);   
   
t=1:length(1997:2018)*12;
%      temp_t_s(isnan(temp_t_s)) = interp1(t(~isnan(temp_t_s)), temp_t_s(~isnan(temp_t_s)), t(isnan(temp_t_s)) ); 
temp_s(isnan(temp_s)) = interp1( t(~isnan(temp_s)), temp_s(~isnan(temp_s)), t(isnan(temp_s)) ); 
temp_b(isnan(temp_b)) = interp1( t(~isnan(temp_b)), temp_b(~isnan(temp_b)), t(isnan(temp_b)) ); 

interp_chl_b(:,i) = temp_b;
interp_chl_s(:,i) = temp_s;

 plot(temp_s,'color','b');
 plot(temp_b,'color','r');  
 ylim([0 20])
xlabel('time(year)','fontsize',13)
ylabel('chl (ug/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:3:end)]);
set(gca,'xlim',[tx_tick(1) tx_tick(end)]);
set(gca,'xticklabel',1997:3:2018);
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',tx_tick(1:1:end));
title([name_tag{i}])
% legend('surface','bottom')
% print('-dpng',[name_tag{i},' chl st.', num2str(i),'.png'])
end

t_tick=[1:12:length(1997:2018)*12]
tx_tick = t_tick;
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
gray = [128/255 128/255 128/255];
for i = 1:8
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

t_tick=[1:12:length(1997:2018)*12]
tx_tick = t_tick;
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
gray = [128/255 128/255 128/255];
for i = 1:8
    clearvars temp_s temp_b
   temp_s= no3_sur_clim(i,:);
   temp_b= no3_bot_clim(i,:);
figure; hold on;
 plot(temp_s,'*','color','b','linew',2);
 plot(temp_b,'*','color','r','linew',2);   
   
t=1:length(1997:2018)*12;
%      temp_t_s(isnan(temp_t_s)) = interp1(t(~isnan(temp_t_s)), temp_t_s(~isnan(temp_t_s)), t(isnan(temp_t_s)) ); 
temp_s(isnan(temp_s)) = interp1( t(~isnan(temp_s)), temp_s(~isnan(temp_s)), t(isnan(temp_s)) ); 
temp_b(isnan(temp_b)) = interp1( t(~isnan(temp_b)), temp_b(~isnan(temp_b)), t(isnan(temp_b)) ); 

interp_no3_b(:,i) = temp_b;
interp_no3_s(:,i) = temp_s;

 plot(temp_s,'color','b');
 plot(temp_b,'color','r');  
 ylim([0 250])
xlabel('time(year)','fontsize',13)
ylabel('no3 (ug/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:3:end)]);
set(gca,'xlim',[tx_tick(1) tx_tick(end)]);
set(gca,'xticklabel',1997:3:2018);
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',tx_tick(1:1:end));
title([name_tag{i}])
% legend('surface','bottom')
% print('-dpng',[name_tag{i},' no3 st.', num2str(i),'.png'])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sig = 3; %sigma
t_tick=[1:12:length(1997:2018)*12]
tx_tick = t_tick;
clearvars b1* X* b* yCalc* nonan* idx* idx1* *_bc_*
% color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
% gray = [128/255 128/255 128/255];
for i = 1:8
    clearvars temp_s temp_b
   temp_s= do_sur_clim(i,:);
   temp_b= do_bot_clim(i,:);
figure; hold on;
 plot(temp_s,'*','color','b','linew',2);
 plot(temp_b,'*','color','r','linew',2);   
 
 t=1:length(1997:2018)*12;
%      temp_t_s(isnan(temp_t_s)) = interp1(t(~isnan(temp_t_s)), temp_t_s(~isnan(temp_t_s)), t(isnan(temp_t_s)) ); 
temp_s(isnan(temp_s)) = interp1( t(~isnan(temp_s)), temp_s(~isnan(temp_s)), t(isnan(temp_s)) ); 
temp_b(isnan(temp_b)) = interp1( t(~isnan(temp_b)), temp_b(~isnan(temp_b)), t(isnan(temp_b)) ); 

 
t_regime = length(1997:2004)*12;
temp_s_front = temp_s(1:t_regime);
temp_s_back = temp_s(t_regime+1:end);
temp_b_front = temp_b(1:t_regime);
temp_b_back = temp_b(t_regime+1:end);

idx_s_front = find(isnan(temp_s_front) == 0);
idx_s_back = find(isnan(temp_s_back) == 0);
idx_b_front = find(isnan(temp_b_front) == 0);
idx_b_back = find(isnan(temp_b_back) == 0);

upper_bc_s_front = mean(temp_s_front(idx_s_front)) + sig*std(temp_s_front(idx_s_front));
lower_bc_s_front = mean(temp_s_front(idx_s_front)) - sig*std(temp_s_front(idx_s_front));
upper_bc_b_front = mean(temp_b_front(idx_s_front)) + sig*std(temp_b_front(idx_b_front));
lower_bc_b_front = mean(temp_b_front(idx_s_front)) - sig*std(temp_b_front(idx_b_front));

upper_bc_s_back = mean(temp_s_back(idx_s_back)) + sig*std(temp_s_back(idx_s_back));
lower_bc_s_back = mean(temp_s_back(idx_s_back)) - sig*std(temp_s_back(idx_s_back));
upper_bc_b_back = mean(temp_b_back(idx_s_back)) + sig*std(temp_b_back(idx_b_back));
lower_bc_b_back = mean(temp_b_back(idx_s_back)) - sig*std(temp_b_back(idx_b_back));

temp_s_front_r = temp_s_front;
temp_s_back_r = temp_s_back;
temp_b_front_r = temp_b_front;
temp_b_back_r = temp_b_back;


temp_s_front_r(find(temp_s_front_r > mean(temp_s_front_r(idx_s_front)) + sig*std(temp_s_front_r(idx_s_front)))) = NaN;
temp_s_back_r(find(temp_s_back_r > mean(temp_s_back_r(idx_s_back)) + sig*std(temp_s_back_r(idx_s_back)))) = NaN; 
temp_b_front_r(find(temp_b_front_r > mean(temp_b_front_r(idx_b_front)) + sig*std(temp_b_front_r(idx_b_front)))) = NaN;
temp_b_back_r(find(temp_b_back_r > mean(temp_b_back_r(idx_b_back)) + sig*std(temp_b_back_r(idx_b_back)))) = NaN; 
temp_s_front_r(find(temp_s_front_r < mean(temp_s_front_r(idx_s_front)) - sig*std(temp_s_front_r(idx_s_front)))) = NaN;
temp_s_back_r(find(temp_s_back_r < mean(temp_s_back_r(idx_s_back)) - sig*std(temp_s_back_r(idx_s_back)))) = NaN; 
temp_b_front_r(find(temp_b_front_r < mean(temp_b_front_r(idx_b_front)) - sig*std(temp_b_front_r(idx_b_front)))) = NaN;
temp_b_back_r(find(temp_b_back_r < mean(temp_b_back_r(idx_b_back)) - sig*std(temp_b_back_r(idx_b_back)))) = NaN; 

% line([1:0.1:t_regime],upper_bc_s_front,'color','b','linew',2);
% line([t_regime+1:0.1:t_regime+length(temp_s_back)],upper_bc_s_back,'color','b','linew',2);
% line([1:0.1:t_regime],upper_bc_b_front,'color','r','linew',2);
% line([t_regime+1:0.1:t_regime+length(temp_b_back)],upper_bc_b_back,'color','r','linew',2);
% line([1:0.1:t_regime],lower_bc_s_front,'color','b','linew',2);
% line([t_regime+1:0.1:t_regime+length(temp_s_back)],lower_bc_s_back,'color','b','linew',2);
% line([1:0.1:t_regime],lower_bc_b_front,'color','r','linew',2);
% line([t_regime+1:0.1:t_regime+length(temp_b_back)],lower_bc_b_back,'color','r','linew',2);
   
%% reigme mean
 plot([1:.1:t_regime],nanmean(temp_s_front_r),'color','b','linew',2);
 plot([t_regime+1:.1:t_regime+length(temp_s_back)],nanmean(temp_s_back_r),'color','b','linew',2);
 plot([1:.1:t_regime],nanmean(temp_b_front_r),'color','r','linew',2);
 plot([t_regime+1:.1:t_regime+length(temp_b_back)],nanmean(temp_b_back_r),'color','r','linew',2); 
 
 ylim([2 16])
xlabel('time(year)','fontsize',13)
ylabel('do (mg/L)','fontsize',13)
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
nonan_s_front=temp_s_front_r(~isnan(temp_s_front_r));
nonan_s_back=temp_s_back_r(~isnan(temp_s_back_r));
nonan_b_front=temp_b_front_r(~isnan(temp_b_front_r));
nonan_b_back=temp_b_back_r(~isnan(temp_b_back_r));

idx1_s_front = find(isnan(temp_s_front_r) == 0);
idx1_s_back = find(isnan(temp_s_back_r) == 0);
idx1_b_front = find(isnan(temp_b_front_r) == 0);
idx1_b_back = find(isnan(temp_b_back_r) == 0);

%% surface
% Slope & Intercept y = b0 + b1*x
X_s_front = [ones(length(idx1_s_front'),1) idx1_s_front']; %b_0 b_1 
b_s_front = X_s_front\nonan_s_front';
b0_s_front = b_s_front(1);  b1_s_front = b_s_front(2);
yCalc_s_front = [1:length(temp_s_front)].*b1_s_front + b0_s_front;

% Slope & Intercept y = b0 + b1*x
X_s_back = [ones(length(idx1_s_back'),1) idx1_s_back']; %b_0 b_1 
b_s_back = X_s_back\nonan_s_back';
b0_s_back = b_s_back(1);  b1_s_back = b_s_back(2);
yCalc_s_back = [length(temp_s_front)+1:length(temp_s_front)+length(temp_s_back)].*b1_s_back + b0_s_back;

%% bottom
% Slope & Intercept y = b0 + b1*x
X_b_front = [ones(length(idx1_b_front'),1) idx1_b_front']; %b_0 b_1 
b_b_front = X_b_front\nonan_b_front';
b0_b_front = b_b_front(1);  b1_b_front = b_b_front(2);
yCalc_b_front = [1:length(temp_b_front)].*b1_b_front + b0_b_front;

% Slope & Intercept y = b0 + b1*x
X_b_back = [ones(length(idx1_b_back'),1) idx1_b_back']; %b_0 b_1 
b_b_back = X_b_back\nonan_b_back';
b0_b_back = b_b_back(1);  b1_b_back = b_b_back(2);
yCalc_b_back = [length(temp_b_front)+1:length(temp_b_front)+length(temp_b_back)].*b1_b_back + b0_b_back;

% plot([1:length(temp_s_front)], yCalc_s_front,'--' ,'color','b','linew',2)
% plot([length(temp_s_front)+1:length(temp_s_front)+length(temp_s_back)], yCalc_s_back,'--' ,'color','b','linew',2)
% plot([1:length(temp_b_front)], yCalc_b_front,'--' ,'color','r','linew',2)
% plot([length(temp_b_front)+1:length(temp_b_front)+length(temp_b_back)], yCalc_b_back,'--' ,'color','r','linew',2)

print('-dpng',[name_tag{i},' regime DO silver st.', num2str(i),'.png'])
end


sig = 3; %sigma
t_tick=[1:12:length(1997:2018)*12]
tx_tick = t_tick;
clearvars b1* X* b* yCalc* nonan* idx* idx1* *_bc_*
% color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
% gray = [128/255 128/255 128/255];
for i = 1:8
    clearvars temp_s temp_b
   temp_s= chl_sur_clim(i,:);
   temp_b= chl_bot_clim(i,:);
figure; hold on;
 plot(temp_s,'*','color','b','linew',2);
 plot(temp_b,'*','color','r','linew',2);   
 
 t=1:length(1997:2018)*12;
%      temp_t_s(isnan(temp_t_s)) = interp1(t(~isnan(temp_t_s)), temp_t_s(~isnan(temp_t_s)), t(isnan(temp_t_s)) ); 
temp_s(isnan(temp_s)) = interp1( t(~isnan(temp_s)), temp_s(~isnan(temp_s)), t(isnan(temp_s)) ); 
temp_b(isnan(temp_b)) = interp1( t(~isnan(temp_b)), temp_b(~isnan(temp_b)), t(isnan(temp_b)) ); 

 
t_regime = length(1997:2004)*12;
temp_s_front = temp_s(1:t_regime);
temp_s_back = temp_s(t_regime+1:end);
temp_b_front = temp_b(1:t_regime);
temp_b_back = temp_b(t_regime+1:end);

idx_s_front = find(isnan(temp_s_front) == 0);
idx_s_back = find(isnan(temp_s_back) == 0);
idx_b_front = find(isnan(temp_b_front) == 0);
idx_b_back = find(isnan(temp_b_back) == 0);

upper_bc_s_front = mean(temp_s_front(idx_s_front)) + sig*std(temp_s_front(idx_s_front));
lower_bc_s_front = mean(temp_s_front(idx_s_front)) - sig*std(temp_s_front(idx_s_front));
upper_bc_b_front = mean(temp_b_front(idx_s_front)) + sig*std(temp_b_front(idx_b_front));
lower_bc_b_front = mean(temp_b_front(idx_s_front)) - sig*std(temp_b_front(idx_b_front));

upper_bc_s_back = mean(temp_s_back(idx_s_back)) + sig*std(temp_s_back(idx_s_back));
lower_bc_s_back = mean(temp_s_back(idx_s_back)) - sig*std(temp_s_back(idx_s_back));
upper_bc_b_back = mean(temp_b_back(idx_s_back)) + sig*std(temp_b_back(idx_b_back));
lower_bc_b_back = mean(temp_b_back(idx_s_back)) - sig*std(temp_b_back(idx_b_back));

temp_s_front_r = temp_s_front;
temp_s_back_r = temp_s_back;
temp_b_front_r = temp_b_front;
temp_b_back_r = temp_b_back;


temp_s_front_r(find(temp_s_front_r > mean(temp_s_front_r(idx_s_front)) + sig*std(temp_s_front_r(idx_s_front)))) = NaN;
temp_s_back_r(find(temp_s_back_r > mean(temp_s_back_r(idx_s_back)) + sig*std(temp_s_back_r(idx_s_back)))) = NaN; 
temp_b_front_r(find(temp_b_front_r > mean(temp_b_front_r(idx_b_front)) + sig*std(temp_b_front_r(idx_b_front)))) = NaN;
temp_b_back_r(find(temp_b_back_r > mean(temp_b_back_r(idx_b_back)) + sig*std(temp_b_back_r(idx_b_back)))) = NaN; 
temp_s_front_r(find(temp_s_front_r < mean(temp_s_front_r(idx_s_front)) - sig*std(temp_s_front_r(idx_s_front)))) = NaN;
temp_s_back_r(find(temp_s_back_r < mean(temp_s_back_r(idx_s_back)) - sig*std(temp_s_back_r(idx_s_back)))) = NaN; 
temp_b_front_r(find(temp_b_front_r < mean(temp_b_front_r(idx_b_front)) - sig*std(temp_b_front_r(idx_b_front)))) = NaN;
temp_b_back_r(find(temp_b_back_r < mean(temp_b_back_r(idx_b_back)) - sig*std(temp_b_back_r(idx_b_back)))) = NaN; 

% line([1:0.1:t_regime],upper_bc_s_front,'color','b','linew',2);
% line([t_regime+1:0.1:t_regime+length(temp_s_back)],upper_bc_s_back,'color','b','linew',2);
% line([1:0.1:t_regime],upper_bc_b_front,'color','r','linew',2);
% line([t_regime+1:0.1:t_regime+length(temp_b_back)],upper_bc_b_back,'color','r','linew',2);
% line([1:0.1:t_regime],lower_bc_s_front,'color','b','linew',2);
% line([t_regime+1:0.1:t_regime+length(temp_s_back)],lower_bc_s_back,'color','b','linew',2);
% line([1:0.1:t_regime],lower_bc_b_front,'color','r','linew',2);
% line([t_regime+1:0.1:t_regime+length(temp_b_back)],lower_bc_b_back,'color','r','linew',2);
   
%% reigme mean
 plot([1:.1:t_regime],nanmean(temp_s_front_r),'color','b','linew',4);
 plot([t_regime+1:.1:t_regime+length(temp_s_back)],nanmean(temp_s_back_r),'color','b','linew',4);
 plot([1:.1:t_regime],nanmean(temp_b_front_r),'color','r','linew',4);
 plot([t_regime+1:.1:t_regime+length(temp_b_back)],nanmean(temp_b_back_r),'color','r','linew',4); 
 
 ylim([0 20])
xlabel('time(year)','fontsize',13)
ylabel('chl (ug/L)','fontsize',13)
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
nonan_s_front=temp_s_front_r(~isnan(temp_s_front_r));
nonan_s_back=temp_s_back_r(~isnan(temp_s_back_r));
nonan_b_front=temp_b_front_r(~isnan(temp_b_front_r));
nonan_b_back=temp_b_back_r(~isnan(temp_b_back_r));

idx1_s_front = find(isnan(temp_s_front_r) == 0);
idx1_s_back = find(isnan(temp_s_back_r) == 0);
idx1_b_front = find(isnan(temp_b_front_r) == 0);
idx1_b_back = find(isnan(temp_b_back_r) == 0);

%% surface
% Slope & Intercept y = b0 + b1*x
X_s_front = [ones(length(idx1_s_front'),1) idx1_s_front']; %b_0 b_1 
b_s_front = X_s_front\nonan_s_front';
b0_s_front = b_s_front(1);  b1_s_front = b_s_front(2);
yCalc_s_front = [1:length(temp_s_front)].*b1_s_front + b0_s_front;

% Slope & Intercept y = b0 + b1*x
X_s_back = [ones(length(idx1_s_back'),1) idx1_s_back']; %b_0 b_1 
b_s_back = X_s_back\nonan_s_back';
b0_s_back = b_s_back(1);  b1_s_back = b_s_back(2);
yCalc_s_back = [length(temp_s_front)+1:length(temp_s_front)+length(temp_s_back)].*b1_s_back + b0_s_back;

%% bottom
% Slope & Intercept y = b0 + b1*x
X_b_front = [ones(length(idx1_b_front'),1) idx1_b_front']; %b_0 b_1 
b_b_front = X_b_front\nonan_b_front';
b0_b_front = b_b_front(1);  b1_b_front = b_b_front(2);
yCalc_b_front = [1:length(temp_b_front)].*b1_b_front + b0_b_front;

% Slope & Intercept y = b0 + b1*x
X_b_back = [ones(length(idx1_b_back'),1) idx1_b_back']; %b_0 b_1 
b_b_back = X_b_back\nonan_b_back';
b0_b_back = b_b_back(1);  b1_b_back = b_b_back(2);
yCalc_b_back = [length(temp_b_front)+1:length(temp_b_front)+length(temp_b_back)].*b1_b_back + b0_b_back;

% plot([1:length(temp_s_front)], yCalc_s_front,'--' ,'color','b','linew',2)
% plot([length(temp_s_front)+1:length(temp_s_front)+length(temp_s_back)], yCalc_s_back,'--' ,'color','b','linew',2)
% plot([1:length(temp_b_front)], yCalc_b_front,'--' ,'color','r','linew',2)
% plot([length(temp_b_front)+1:length(temp_b_front)+length(temp_b_back)], yCalc_b_back,'--' ,'color','r','linew',2)

print('-dpng',[name_tag{i},' regime chl silver st.', num2str(i),'.png'])
end


sig = 3; %sigma
t_tick=[1:12:length(1997:2018)*12]
tx_tick = t_tick;
clearvars b1* X* b* yCalc* nonan* idx* idx1* *_bc_*
% color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
% gray = [128/255 128/255 128/255];
for i = 1:8
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

 
t_regime = length(1997:2004)*12;
temp_s_front = temp_s(1:t_regime);
temp_s_back = temp_s(t_regime+1:end);
temp_b_front = temp_b(1:t_regime);
temp_b_back = temp_b(t_regime+1:end);

idx_s_front = find(isnan(temp_s_front) == 0);
idx_s_back = find(isnan(temp_s_back) == 0);
idx_b_front = find(isnan(temp_b_front) == 0);
idx_b_back = find(isnan(temp_b_back) == 0);

upper_bc_s_front = mean(temp_s_front(idx_s_front)) + sig*std(temp_s_front(idx_s_front));
lower_bc_s_front = mean(temp_s_front(idx_s_front)) - sig*std(temp_s_front(idx_s_front));
upper_bc_b_front = mean(temp_b_front(idx_s_front)) + sig*std(temp_b_front(idx_b_front));
lower_bc_b_front = mean(temp_b_front(idx_s_front)) - sig*std(temp_b_front(idx_b_front));

upper_bc_s_back = mean(temp_s_back(idx_s_back)) + sig*std(temp_s_back(idx_s_back));
lower_bc_s_back = mean(temp_s_back(idx_s_back)) - sig*std(temp_s_back(idx_s_back));
upper_bc_b_back = mean(temp_b_back(idx_s_back)) + sig*std(temp_b_back(idx_b_back));
lower_bc_b_back = mean(temp_b_back(idx_s_back)) - sig*std(temp_b_back(idx_b_back));

temp_s_front_r = temp_s_front;
temp_s_back_r = temp_s_back;
temp_b_front_r = temp_b_front;
temp_b_back_r = temp_b_back;


temp_s_front_r(find(temp_s_front_r > mean(temp_s_front_r(idx_s_front)) + sig*std(temp_s_front_r(idx_s_front)))) = NaN;
temp_s_back_r(find(temp_s_back_r > mean(temp_s_back_r(idx_s_back)) + sig*std(temp_s_back_r(idx_s_back)))) = NaN; 
temp_b_front_r(find(temp_b_front_r > mean(temp_b_front_r(idx_b_front)) + sig*std(temp_b_front_r(idx_b_front)))) = NaN;
temp_b_back_r(find(temp_b_back_r > mean(temp_b_back_r(idx_b_back)) + sig*std(temp_b_back_r(idx_b_back)))) = NaN; 
temp_s_front_r(find(temp_s_front_r < mean(temp_s_front_r(idx_s_front)) - sig*std(temp_s_front_r(idx_s_front)))) = NaN;
temp_s_back_r(find(temp_s_back_r < mean(temp_s_back_r(idx_s_back)) - sig*std(temp_s_back_r(idx_s_back)))) = NaN; 
temp_b_front_r(find(temp_b_front_r < mean(temp_b_front_r(idx_b_front)) - sig*std(temp_b_front_r(idx_b_front)))) = NaN;
temp_b_back_r(find(temp_b_back_r < mean(temp_b_back_r(idx_b_back)) - sig*std(temp_b_back_r(idx_b_back)))) = NaN; 

% line([1:0.1:t_regime],upper_bc_s_front,'color','b','linew',2);
% line([t_regime+1:0.1:t_regime+length(temp_s_back)],upper_bc_s_back,'color','b','linew',2);
% line([1:0.1:t_regime],upper_bc_b_front,'color','r','linew',2);
% line([t_regime+1:0.1:t_regime+length(temp_b_back)],upper_bc_b_back,'color','r','linew',2);
% line([1:0.1:t_regime],lower_bc_s_front,'color','b','linew',2);
% line([t_regime+1:0.1:t_regime+length(temp_s_back)],lower_bc_s_back,'color','b','linew',2);
% line([1:0.1:t_regime],lower_bc_b_front,'color','r','linew',2);
% line([t_regime+1:0.1:t_regime+length(temp_b_back)],lower_bc_b_back,'color','r','linew',2);
   
%% reigme mean
 plot([1:.1:t_regime],nanmean(temp_s_front_r),'color','b','linew',2);
 plot([t_regime+1:.1:t_regime+length(temp_s_back)],nanmean(temp_s_back_r),'color','b','linew',2);
 plot([1:.1:t_regime],nanmean(temp_b_front_r),'color','r','linew',2);
 plot([t_regime+1:.1:t_regime+length(temp_b_back)],nanmean(temp_b_back_r),'color','r','linew',2); 
 
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
nonan_s_front=temp_s_front_r(~isnan(temp_s_front_r));
nonan_s_back=temp_s_back_r(~isnan(temp_s_back_r));
nonan_b_front=temp_b_front_r(~isnan(temp_b_front_r));
nonan_b_back=temp_b_back_r(~isnan(temp_b_back_r));

idx1_s_front = find(isnan(temp_s_front_r) == 0);
idx1_s_back = find(isnan(temp_s_back_r) == 0);
idx1_b_front = find(isnan(temp_b_front_r) == 0);
idx1_b_back = find(isnan(temp_b_back_r) == 0);

%% surface
% Slope & Intercept y = b0 + b1*x
X_s_front = [ones(length(idx1_s_front'),1) idx1_s_front']; %b_0 b_1 
b_s_front = X_s_front\nonan_s_front';
b0_s_front = b_s_front(1);  b1_s_front = b_s_front(2);
yCalc_s_front = [1:length(temp_s_front)].*b1_s_front + b0_s_front;

% Slope & Intercept y = b0 + b1*x
X_s_back = [ones(length(idx1_s_back'),1) idx1_s_back']; %b_0 b_1 
b_s_back = X_s_back\nonan_s_back';
b0_s_back = b_s_back(1);  b1_s_back = b_s_back(2);
yCalc_s_back = [length(temp_s_front)+1:length(temp_s_front)+length(temp_s_back)].*b1_s_back + b0_s_back;

%% bottom
% Slope & Intercept y = b0 + b1*x
X_b_front = [ones(length(idx1_b_front'),1) idx1_b_front']; %b_0 b_1 
b_b_front = X_b_front\nonan_b_front';
b0_b_front = b_b_front(1);  b1_b_front = b_b_front(2);
yCalc_b_front = [1:length(temp_b_front)].*b1_b_front + b0_b_front;

% Slope & Intercept y = b0 + b1*x
X_b_back = [ones(length(idx1_b_back'),1) idx1_b_back']; %b_0 b_1 
b_b_back = X_b_back\nonan_b_back';
b0_b_back = b_b_back(1);  b1_b_back = b_b_back(2);
yCalc_b_back = [length(temp_b_front)+1:length(temp_b_front)+length(temp_b_back)].*b1_b_back + b0_b_back;

% plot([1:length(temp_s_front)], yCalc_s_front,'--' ,'color','b','linew',2)
% plot([length(temp_s_front)+1:length(temp_s_front)+length(temp_s_back)], yCalc_s_back,'--' ,'color','b','linew',2)
% plot([1:length(temp_b_front)], yCalc_b_front,'--' ,'color','r','linew',2)
% plot([length(temp_b_front)+1:length(temp_b_front)+length(temp_b_back)], yCalc_b_back,'--' ,'color','r','linew',2)

print('-dpng',[name_tag{i},' regime nh4 silver st.', num2str(i),'.png'])
end


sig = 3; %sigma
t_tick=[1:12:length(1997:2018)*12]
tx_tick = t_tick;
clearvars b1* X* b* yCalc* nonan* idx* idx1* *_bc_*
% color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
% gray = [128/255 128/255 128/255];
for i = 1:8
    clearvars temp_s temp_b
   temp_s= no3_sur_clim(i,:);
   temp_b= no3_bot_clim(i,:);
figure; hold on;
 plot(temp_s,'*','color','b','linew',2);
 plot(temp_b,'*','color','r','linew',2);   
 
 t=1:length(1997:2018)*12;
%      temp_t_s(isnan(temp_t_s)) = interp1(t(~isnan(temp_t_s)), temp_t_s(~isnan(temp_t_s)), t(isnan(temp_t_s)) ); 
temp_s(isnan(temp_s)) = interp1( t(~isnan(temp_s)), temp_s(~isnan(temp_s)), t(isnan(temp_s)) ); 
temp_b(isnan(temp_b)) = interp1( t(~isnan(temp_b)), temp_b(~isnan(temp_b)), t(isnan(temp_b)) ); 

 
t_regime = length(1997:2004)*12;
temp_s_front = temp_s(1:t_regime);
temp_s_back = temp_s(t_regime+1:end);
temp_b_front = temp_b(1:t_regime);
temp_b_back = temp_b(t_regime+1:end);

idx_s_front = find(isnan(temp_s_front) == 0);
idx_s_back = find(isnan(temp_s_back) == 0);
idx_b_front = find(isnan(temp_b_front) == 0);
idx_b_back = find(isnan(temp_b_back) == 0);

upper_bc_s_front = mean(temp_s_front(idx_s_front)) + sig*std(temp_s_front(idx_s_front));
lower_bc_s_front = mean(temp_s_front(idx_s_front)) - sig*std(temp_s_front(idx_s_front));
upper_bc_b_front = mean(temp_b_front(idx_s_front)) + sig*std(temp_b_front(idx_b_front));
lower_bc_b_front = mean(temp_b_front(idx_s_front)) - sig*std(temp_b_front(idx_b_front));

upper_bc_s_back = mean(temp_s_back(idx_s_back)) + sig*std(temp_s_back(idx_s_back));
lower_bc_s_back = mean(temp_s_back(idx_s_back)) - sig*std(temp_s_back(idx_s_back));
upper_bc_b_back = mean(temp_b_back(idx_s_back)) + sig*std(temp_b_back(idx_b_back));
lower_bc_b_back = mean(temp_b_back(idx_s_back)) - sig*std(temp_b_back(idx_b_back));

temp_s_front_r = temp_s_front;
temp_s_back_r = temp_s_back;
temp_b_front_r = temp_b_front;
temp_b_back_r = temp_b_back;


temp_s_front_r(find(temp_s_front_r > mean(temp_s_front_r(idx_s_front)) + sig*std(temp_s_front_r(idx_s_front)))) = NaN;
temp_s_back_r(find(temp_s_back_r > mean(temp_s_back_r(idx_s_back)) + sig*std(temp_s_back_r(idx_s_back)))) = NaN; 
temp_b_front_r(find(temp_b_front_r > mean(temp_b_front_r(idx_b_front)) + sig*std(temp_b_front_r(idx_b_front)))) = NaN;
temp_b_back_r(find(temp_b_back_r > mean(temp_b_back_r(idx_b_back)) + sig*std(temp_b_back_r(idx_b_back)))) = NaN; 
temp_s_front_r(find(temp_s_front_r < mean(temp_s_front_r(idx_s_front)) - sig*std(temp_s_front_r(idx_s_front)))) = NaN;
temp_s_back_r(find(temp_s_back_r < mean(temp_s_back_r(idx_s_back)) - sig*std(temp_s_back_r(idx_s_back)))) = NaN; 
temp_b_front_r(find(temp_b_front_r < mean(temp_b_front_r(idx_b_front)) - sig*std(temp_b_front_r(idx_b_front)))) = NaN;
temp_b_back_r(find(temp_b_back_r < mean(temp_b_back_r(idx_b_back)) - sig*std(temp_b_back_r(idx_b_back)))) = NaN; 

% line([1:0.1:t_regime],upper_bc_s_front,'color','b','linew',2);
% line([t_regime+1:0.1:t_regime+length(temp_s_back)],upper_bc_s_back,'color','b','linew',2);
% line([1:0.1:t_regime],upper_bc_b_front,'color','r','linew',2);
% line([t_regime+1:0.1:t_regime+length(temp_b_back)],upper_bc_b_back,'color','r','linew',2);
% line([1:0.1:t_regime],lower_bc_s_front,'color','b','linew',2);
% line([t_regime+1:0.1:t_regime+length(temp_s_back)],lower_bc_s_back,'color','b','linew',2);
% line([1:0.1:t_regime],lower_bc_b_front,'color','r','linew',2);
% line([t_regime+1:0.1:t_regime+length(temp_b_back)],lower_bc_b_back,'color','r','linew',2);
%    
%% reigme mean
 plot([1:.1:t_regime],nanmean(temp_s_front_r),'color','b','linew',4);
 plot([t_regime+1:.1:t_regime+length(temp_s_back)],nanmean(temp_s_back_r),'color','b','linew',4);
 plot([1:.1:t_regime],nanmean(temp_b_front_r),'color','r','linew',4);
 plot([t_regime+1:.1:t_regime+length(temp_b_back)],nanmean(temp_b_back_r),'color','r','linew',4); 
 
 ylim([0 250])
xlabel('time(year)','fontsize',13)
ylabel('no3 (ug/L)','fontsize',13)
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
nonan_s_front=temp_s_front_r(~isnan(temp_s_front_r));
nonan_s_back=temp_s_back_r(~isnan(temp_s_back_r));
nonan_b_front=temp_b_front_r(~isnan(temp_b_front_r));
nonan_b_back=temp_b_back_r(~isnan(temp_b_back_r));

idx1_s_front = find(isnan(temp_s_front_r) == 0);
idx1_s_back = find(isnan(temp_s_back_r) == 0);
idx1_b_front = find(isnan(temp_b_front_r) == 0);
idx1_b_back = find(isnan(temp_b_back_r) == 0);

%% surface
% Slope & Intercept y = b0 + b1*x
X_s_front = [ones(length(idx1_s_front'),1) idx1_s_front']; %b_0 b_1 
b_s_front = X_s_front\nonan_s_front';
b0_s_front = b_s_front(1);  b1_s_front = b_s_front(2);
yCalc_s_front = [1:length(temp_s_front)].*b1_s_front + b0_s_front;

% Slope & Intercept y = b0 + b1*x
X_s_back = [ones(length(idx1_s_back'),1) idx1_s_back']; %b_0 b_1 
b_s_back = X_s_back\nonan_s_back';
b0_s_back = b_s_back(1);  b1_s_back = b_s_back(2);
yCalc_s_back = [length(temp_s_front)+1:length(temp_s_front)+length(temp_s_back)].*b1_s_back + b0_s_back;

%% bottom
% Slope & Intercept y = b0 + b1*x
X_b_front = [ones(length(idx1_b_front'),1) idx1_b_front']; %b_0 b_1 
b_b_front = X_b_front\nonan_b_front';
b0_b_front = b_b_front(1);  b1_b_front = b_b_front(2);
yCalc_b_front = [1:length(temp_b_front)].*b1_b_front + b0_b_front;

% Slope & Intercept y = b0 + b1*x
X_b_back = [ones(length(idx1_b_back'),1) idx1_b_back']; %b_0 b_1 
b_b_back = X_b_back\nonan_b_back';
b0_b_back = b_b_back(1);  b1_b_back = b_b_back(2);
yCalc_b_back = [length(temp_b_front)+1:length(temp_b_front)+length(temp_b_back)].*b1_b_back + b0_b_back;

% plot([1:length(temp_s_front)], yCalc_s_front,'--' ,'color','b','linew',2)
% plot([length(temp_s_front)+1:length(temp_s_front)+length(temp_s_back)], yCalc_s_back,'--' ,'color','b','linew',2)
% plot([1:length(temp_b_front)], yCalc_b_front,'--' ,'color','r','linew',2)
% plot([length(temp_b_front)+1:length(temp_b_front)+length(temp_b_back)], yCalc_b_back,'--' ,'color','r','linew',2)

print('-dpng',[name_tag{i},' regime no3 silver st.', num2str(i),'.png'])
end



