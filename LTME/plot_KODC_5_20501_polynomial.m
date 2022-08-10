close all; clear; clc;   % -v3

% % station name tag
% % (for preventing the miss picking using tag from 1 col's name on excel) 

% port
name_tag_1{1} = [num2str(01,'%02d'),'번'] 
 
% combining the tag and outter point excluding
% name_tag = name_tag_1'; 
% name_tag{end+1:end+length(name_tag_2)} = name_tag_2; 
% name_tag{end+1:end+length(name_tag_3)} = name_tag_3; 
% size_tag = length(name_tag);


%% pick the row on the excel which has same name with tag
[raw_p txt_p]=xlsread('정선해양조사_205line_from80.xls','sheet','');
txt_matc_p = txt_p(3:end,2); % name list
txt_date_p = txt_p(3:end,3); % date list
txt_cut=txt_p(3:end,:);

depth_nan_p = strcmp(txt_cut(:,4),''); % detect NaN
for i = 1:length(txt_cut)
   if depth_nan_p(i) == 1 % when it's nan
       depth_p(i) = NaN;
   else
       depth_p(i) = str2num(char(txt_cut(i,4)));
   end
end

temp_nan_p = strcmp(txt_cut(:,5),''); % detect NaN
for i = 1:length(txt_cut)
   if temp_nan_p(i) == 1 % when it's nan
       temp_p(i) = NaN;
   else
       temp_p(i) = str2num(char(txt_cut(i,5)));
   end
end

salt_nan_p = strcmp(txt_cut(:,7),''); % detect NaN
for i = 1:length(txt_cut)
   if salt_nan_p(i) == 1 % when it's nan
       salt_p(i) = NaN;
   else
       salt_p(i) = str2num(char(txt_cut(i,7)));
   end
end

do_nan_p = strcmp(txt_cut(:,9),''); % detect NaN
for i = 1:length(txt_cut)
   if do_nan_p(i) == 1 % when it's nan
       do_p(i) = NaN;
   else
       do_p(i) = str2num(char(txt_cut(i,9)));
   end
end

no3_nan_p = strcmp(txt_cut(:,14),''); % detect NaN
for i = 1:length(txt_cut)
   if no3_nan_p(i) == 1 % when it's nan
       no3_p(i) = NaN;
   else
       no3_p(i) = str2num(char(txt_cut(i,14)));
   end
end

%% pick matched name with tag
% station & date
for i = 1:length(name_tag_1)
   if  sum(strcmp(name_tag_1{i}, txt_matc_p)) ~= 0
       indx_st{i} = find([strcmp(name_tag_1{i}, txt_matc_p)] == 1)
       date_st{i} = txt_date_p(indx_st{i});
   end
end

for i = 1:length(name_tag_1)
    depth_kodc{i} = depth_p(indx_st{i});
    temp_kodc{i} = temp_p(indx_st{i});
    salt_kodc{i} = salt_p(indx_st{i});
    do_kodc{i} = do_p(indx_st{i});
    no3_kodc{i} = no3_p(indx_st{i});
end


%% make date to be 'yymm' form
for i = 1:length(name_tag_1)
clearvars temp temp_ymd temp_ymd_c temp_yymmdd temp_yymmdd_c
temp = char(date_st{i});
temp_ymd=temp(:,1:end-12);
temp_yymmdd=temp(:,1:end-9);

for j = 1:length(temp_ymd)
    temp_ymd_c{j} = temp_ymd(j,:);
    temp_yymmdd_c{j} = temp_yymmdd(j,:);
                    
end
date_ymd{i,1} = temp_ymd_c;
date_yymmdd{i,1} = temp_yymmdd_c;
end

% make 1980~present
k=0
for i = 1980:2019
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

k=0; m=0;
for i = 1:40
    l=0
        ref_yy(i,:)=[num2str(i+1979)];
    for n = 1:12
        m = m+1;
        ref_yymm{m}=[num2str(i+1979) '-' num2str(n,'%02d')];
    for j = 1:eom_d(i,n)
        l=l+1; % index for 1~31(or 30 or 29)
        k=k+1; % index for 1~366
        ref_yymmdd{k}=[num2str(i+1979) '-' num2str(n,'%02d') '-'  num2str(j,'%02d')];
        ref_mmdd{k}=[num2str(n,'%02d') '-'  num2str(j,'%02d')];
    end
    end
end


% matched date 'yymm' form on surf & bottom
for j = 1:length(name_tag_1) % st. axis
    for i = 1:length(ref_yymmdd) % date axis
       if  sum(strcmp(ref_yymmdd{i}, date_ymd{j})) ~= 0
           clearvars indx_date
           indx_date = find([strcmp(ref_yymmdd{i}, date_yymmdd{j})] == 1); 
           indx_date_s{j,i} = indx_date(find(depth_kodc{j}(1,indx_date) == 0));
           indx_date_b{j,i} = indx_date(find(depth_kodc{j}(1,indx_date) == max(depth_kodc{j}(1,indx_date))));
       end
    end 
end

indx_date_s = cell(1,length(ref_yymmdd));
indx_date_b = cell(1,length(ref_yymmdd));
for j = 1:length(name_tag_1) % st. axis
    for i = 1:length(ref_yymmdd) % date axis
       if  sum(strcmp(ref_yymmdd{i}, date_yymmdd{j})) ~= 0
           clearvars indx_date
           indx_date = find([strcmp(ref_yymmdd{i}, date_yymmdd{j})] == 1); 
           indx_date_s{j,i} = indx_date(find(depth_kodc{j}(1,indx_date) == 0));
           indx_date_b{j,i} = indx_date(find(depth_kodc{j}(1,indx_date) == max(depth_kodc{j}(1,indx_date))));
%        elseif sum(strcmp(ref_yymmdd{i}, date_yymmdd{j})) == 0
%            indx_date_s{j,i} = NaN;
%            indx_date_b{j,i} = NaN;
       end
    end
end


% clearvars indx_date_*
% for j = 1:3 % st. axis
%     for i = 1:length(ref_yymm) % date axis
%        if  sum(strcmp(ref_yymm{i}, date_ymd{j})) ~= 0
%            clearvars indx_date
%            indx_date = find([strcmp(ref_yymm{i}, date_ymd{j})] == 1);
%            if sum(indx_date(find(depth_kodc{j}(1,indx_date)==0))) ~= 0
%             indx_date_s(j,i) = indx_date(find(depth_kodc{j}(1,indx_date)==0));
%            elseif  sum(indx_date(find(depth_kodc{j}(1,indx_date)==0))) == 0
%             indx_date_s(j,i) = NaN;
%            end
%            
% %            indx_date_b(j,i) = indx_date(find(depth_kodc{j}(1,indx_date)==max(depth_kodc{j}(1,indx_date))));
%        end
%     end
% end


% make climate

%temp
for i = 1:size(indx_date_s,1)
    clearvars temp_kodc_k
    temp_kodc_k = temp_kodc{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_kodc_k(indx_date_s{i,j}))) ~= 0     
        temp_sur(i,j) = nanmean(temp_kodc_k(indx_date_s{i,j}));
    elseif sum(size(temp_kodc_k(indx_date_s{i,j}))) == 0 
        temp_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_kodc_k
    temp_kodc_k = temp_kodc{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_kodc_k(indx_date_b{i,j}))) ~= 0     
        temp_bot(i,j) = nanmean(temp_kodc_k(indx_date_b{i,j}));
    elseif sum(size(temp_kodc_k(indx_date_b{i,j}))) == 0 
        temp_bot(i,j) = NaN;
    end
    end
end

%salt
for i = 1:size(indx_date_s,1)
    clearvars temp_kodc_k
    temp_kodc_k = salt_kodc{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_kodc_k(indx_date_s{i,j}))) ~= 0     
        salt_sur(i,j) = nanmean(temp_kodc_k(indx_date_s{i,j}));
    elseif sum(size(temp_kodc_k(indx_date_s{i,j}))) == 0 
        salt_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_kodc_k
    temp_kodc_k = salt_kodc{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_kodc_k(indx_date_b{i,j}))) ~= 0     
        salt_bot(i,j) = nanmean(temp_kodc_k(indx_date_b{i,j}));
    elseif sum(size(temp_kodc_k(indx_date_b{i,j}))) == 0 
        salt_bot(i,j) = NaN;
    end
    end
end

%do
for i = 1:size(indx_date_s,1)
    clearvars temp_kodc_k
    temp_kodc_k = do_kodc{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_kodc_k(indx_date_s{i,j}))) ~= 0     
        do_sur(i,j) = nanmean(temp_kodc_k(indx_date_s{i,j}));
    elseif sum(size(temp_kodc_k(indx_date_s{i,j}))) == 0 
        do_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_kodc_k
    temp_kodc_k = do_kodc{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_kodc_k(indx_date_b{i,j}))) ~= 0     
        do_bot(i,j) = nanmean(temp_kodc_k(indx_date_b{i,j}));
    elseif sum(size(temp_kodc_k(indx_date_b{i,j}))) == 0 
        do_bot(i,j) = NaN;
    end
    end
end


%no3
for i = 1:size(indx_date_s,1)
    clearvars temp_kodc_k
    temp_kodc_k = no3_kodc{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_kodc_k(indx_date_s{i,j}))) ~= 0     
        no3_sur(i,j) = nanmean(temp_kodc_k(indx_date_s{i,j}));
    elseif sum(size(temp_kodc_k(indx_date_s{i,j}))) == 0 
        no3_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_kodc_k
    temp_kodc_k = no3_kodc{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_kodc_k(indx_date_b{i,j}))) ~= 0     
        no3_bot(i,j) = nanmean(temp_kodc_k(indx_date_b{i,j}));
    elseif sum(size(temp_kodc_k(indx_date_b{i,j}))) == 0 
        no3_bot(i,j) = NaN;
    end
    end
end

% check when they were overapped
%surface
for i = 1:size(indx_date_s,1)
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_p(indx_date_s{i,j}))) > 2     
       disp(['i val = ', num2str(i)]); disp(['j val = ', num2str(j)]);
    end
    end
end

%bottom
for i = 1:size(indx_date_b,1)
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_p(indx_date_b{i,j}))) > 2     
       disp(['i val = ', num2str(i)]); disp(['j val = ', num2str(j)]);
    end
    end
end

% save('KODC_data_monthly_20501.mat','*_bot','*_sur','ref_yymm', 'date_ymd');


%% interp >> to do regime shift test

do_sur = do_sur(1,:)';
do_bot = do_bot(1,:)';
no3_sur =  no3_sur(1,:)'.* 14;  %umol/L -> ug/L
no3_bot = no3_bot(1,:)' .* 14;
temp_sur = temp_sur(1,:)';
temp_bot = temp_bot(1,:)';
salt_sur = salt_sur(1,:)';
salt_bot = salt_bot(1,:)';

%% extract over 3sig
sig = 3; %% sigma
% sur
clearvars idx1
idx1 = find(isnan(do_sur) == 0);
regime_do=do_sur;
regime_do(find(do_sur > mean(do_sur(idx1)) + sig*std(do_sur(idx1))))=NaN;
regime_do(find(do_sur < mean(do_sur(idx1)) - sig*std(do_sur(idx1))))=NaN;
regm_do = nanmean(regime_do)

clearvars idx1
idx1 = find(isnan(no3_sur) == 0);
regime_no3=no3_sur;
regime_no3(find(no3_sur > mean(no3_sur(idx1)) + sig*std(no3_sur(idx1))))=NaN;
regime_no3(find(no3_sur < mean(no3_sur(idx1)) - sig*std(no3_sur(idx1))))=NaN;
regm_no3 = nanmean(regime_no3)


clearvars idx1
idx1 = find(isnan(temp_sur) == 0);
regime_temp=temp_sur;
regime_temp(find(temp_sur > mean(temp_sur(idx1)) + sig*std(temp_sur(idx1))))=NaN;
regime_temp(find(temp_sur < mean(temp_sur(idx1)) - sig*std(temp_sur(idx1))))=NaN;
regm_temp = nanmean(regime_temp)


clearvars idx1
idx1 = find(isnan(salt_sur) == 0);
regime_salt=salt_sur;
regime_salt(find(salt_sur > mean(salt_sur(idx1)) + sig*std(salt_sur(idx1))))=NaN;
regime_salt(find(salt_sur < mean(salt_sur(idx1)) - sig*std(salt_sur(idx1))))=NaN;
regm_salt = nanmean(regime_salt)

%bot
clearvars idx1
idx1 = find(isnan(do_bot) == 0);
regime_do_b=do_bot;
regime_do_b(find(do_bot > mean(do_bot(idx1)) + sig*std(do_bot(idx1))))=NaN;
regime_do_b(find(do_bot < mean(do_bot(idx1)) - sig*std(do_bot(idx1))))=NaN;
regm_do_b = nanmean(regime_do_b)

clearvars idx1
idx1 = find(isnan(no3_bot) == 0);
regime_no3_b=no3_bot;
regime_no3_b(find(no3_bot > mean(no3_bot(idx1)) + sig*std(no3_bot(idx1))))=NaN;
regime_no3_b(find(no3_bot < mean(no3_bot(idx1)) - sig*std(no3_bot(idx1))))=NaN;
regm_no3_b = nanmean(regime_no3_b)


clearvars idx1
idx1 = find(isnan(temp_bot) == 0);
regime_temp_b=temp_bot;
regime_temp_b(find(temp_bot > mean(temp_bot(idx1)) + sig*std(temp_bot(idx1))))=NaN;
regime_temp_b(find(temp_bot < mean(temp_bot(idx1)) - sig*std(temp_bot(idx1))))=NaN;
regm_temp_b = nanmean(regime_temp)


clearvars idx1
idx1 = find(isnan(salt_bot) == 0);
regime_salt_b=salt_bot;
regime_salt_b(find(salt_bot > mean(salt_bot(idx1)) + sig*std(salt_bot(idx1))))=NaN;
regime_salt_b(find(salt_bot < mean(salt_bot(idx1)) - sig*std(salt_bot(idx1))))=NaN;
regm_salt_b = nanmean(regime_salt)

%% extract regime mean
regime_do = regime_do - regm_do;
regime_no3 = regime_no3 - regm_no3;
regime_temp = regime_temp - regm_temp;
regime_salt = regime_salt - regm_salt;

regime_do_b = regime_do_b - regm_do_b;
regime_no3_b = regime_no3_b - regm_no3_b;
regime_temp_b = regime_temp_b - regm_temp_b;
regime_salt_b = regime_salt_b - regm_salt_b;

%% matched 366

% make 366 mm-dd
for i = 1:12
  eom_d(i) = eomday(1980,i); % 1980 is leap-yr
end

k=0
for i = 1:12
    l=0
    for j = 1:eom_d(i)
        l=l+1; % index for 1~31(or 30 or 29)
        k=k+1; % index for 1~366
        mmdd(k,:)=[num2str(i,'%02d') '-'  num2str(l,'%02d')]
    end
end

% make it to cell-array
for i = 1:length(mmdd)
    mmdd_c{i,1} = mmdd(i,:); % delete year
end

% pick matched date from water temp date
for i = 1:length(mmdd_c)
       mmdd_indx{i} = find([strcmp(mmdd_c{i}, ref_mmdd)] == 1)
end


%confirm how may data on each the days
size_c = [];
for i = 1:366; size_c(i)=length(mmdd_indx{i}); end
bar(size_c)

% make quasi_climate 1980~2019
for i = 1:length(mmdd_indx)
    if size(mmdd_indx{i},1) == 0     
        reg_clim_do(i) = NaN;
        reg_clim_no3(i) = NaN;
        reg_clim_temp(i) = NaN; 
        reg_clim_salt(i) = NaN;
        reg_clim_do_b(i) = NaN;
        reg_clim_no3_b(i) = NaN;
        reg_clim_temp_b(i) = NaN; 
        reg_clim_salt_b(i) = NaN;
    else
        reg_clim_do(i) = nanmean(regime_do(mmdd_indx{i}));
        reg_clim_no3(i) = nanmean(regime_no3(mmdd_indx{i}));
        reg_clim_temp(i) = nanmean(regime_temp(mmdd_indx{i}));
        reg_clim_salt(i) = nanmean(regime_salt(mmdd_indx{i}));
        reg_clim_do_b(i) = nanmean(regime_do(mmdd_indx{i}));
        reg_clim_no3_b(i) = nanmean(regime_no3(mmdd_indx{i}));
        reg_clim_temp_b(i) = nanmean(regime_temp(mmdd_indx{i}));
        reg_clim_salt_b(i) = nanmean(regime_salt(mmdd_indx{i}));
    end
end


%surface
%DO
clearvars xp_w_do pf_w_do
xp = 1:366;
xp_w_do = find(isnan(reg_clim_do)==0);
pf_w_do = polyfit(xp_w_do,reg_clim_do(xp_w_do),3);
yp_w_do_04 = polyval(pf_w_do,xp);
yp_w_do_04 = yp_w_do_04 + regm_do;

%NO3
clearvars xp_w_no3 pf_w_no3
xp = 1:366;
xp_w_no3 = find(isnan(reg_clim_no3)==0);
pf_w_no3 = polyfit(xp_w_no3,reg_clim_no3(xp_w_no3),3);
yp_w_no3_04 = polyval(pf_w_no3,xp);
yp_w_no3_04 = yp_w_no3_04 + regm_no3;

%salt
clearvars xp_w_salt pf_w_salt
xp = 1:366;
xp_w_salt = find(isnan(reg_clim_salt)==0);
pf_w_salt = polyfit(xp_w_salt,reg_clim_salt(xp_w_salt),3);
yp_w_salt_04 = polyval(pf_w_salt,xp);
yp_w_salt_04 = yp_w_salt_04 + regm_salt;

%temp
clearvars xp_w_temp pf_w_temp
xp = 1:366;
xp_w_temp = find(isnan(reg_clim_temp)==0);
pf_w_temp = polyfit(xp_w_temp,reg_clim_temp(xp_w_temp),3);
yp_w_temp_04 = polyval(pf_w_temp,xp);
yp_w_temp_04 = yp_w_temp_04 + regm_temp;


%bottom
%DO
clearvars xp_w_do_b pf_w_do_b
xp = 1:366;
xp_w_do_b = find(isnan(reg_clim_do_b)==0);
pf_w_do_b = polyfit(xp_w_do_b,reg_clim_do_b(xp_w_do_b),3);
yp_w_do_04_b = polyval(pf_w_do_b,xp);
yp_w_do_04_b = yp_w_do_04_b + regm_do_b;

%NO3
clearvars xp_w_no3_b pf_w_no3_b
xp = 1:366;
xp_w_no3_b = find(isnan(reg_clim_no3_b)==0);
pf_w_no3_b = polyfit(xp_w_no3_b,reg_clim_no3_b(xp_w_no3_b),3);
yp_w_no3_04_b = polyval(pf_w_no3_b,xp);
yp_w_no3_04_b = yp_w_no3_04_b + regm_no3_b;

%salt
clearvars xp_w_salt_b pf_w_salt_b
xp = 1:366;
xp_w_salt_b = find(isnan(reg_clim_salt_b)==0);
pf_w_salt_b = polyfit(xp_w_salt_b,reg_clim_salt_b(xp_w_salt_b),3);
yp_w_salt_04_b = polyval(pf_w_salt_b,xp);
yp_w_salt_04_b = yp_w_salt_04_b + regm_salt_b;

%temp
clearvars xp_w_temp_b pf_w_temp_b
xp = 1:366;
xp_w_temp_b = find(isnan(reg_clim_temp_b)==0);
pf_w_temp_b = polyfit(xp_w_temp_b,reg_clim_temp_b(xp_w_temp_b),3);
yp_w_temp_04_b = polyval(pf_w_temp_b,xp);
yp_w_temp_04_b = yp_w_temp_04_b + regm_temp_b;

return

%% DO
figure;
scatter(1:366,reg_clim_do + regm_do);
hold on
plot(1:366, yp_w_do_04,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('DO (mg/L)','fontsize',13)
title('정선관측(표층)-polynomial DO.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([2 18]) 
xlim([1 366])

%% salt
figure;
scatter(1:366,reg_clim_salt + regm_salt);
hold on
plot(1:366, yp_w_salt_04,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('salta (mg/m^3)','fontsize',13)
title('sumjin(songjung)-polynomial salt.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([32 35])
xlim([1 366])


%% no3
figure;
scatter(1:366,reg_clim_no3 + regm_no3);
hold on
plot(1:366, yp_w_no3_04,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('NO3 (ug/L)','fontsize',13)
title('정선관측(표층)-polynomial no3.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([0 inf])
xlim([1 366])

%% temp
figure;
scatter(1:366,reg_clim_temp + regm_temp);
hold on
plot(1:366, yp_w_temp_04,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('temp (mg/L)','fontsize',13)
title('sumjin(songjung)-polynomial temp.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([0 inf])
xlim([1 366])



%bottom
%% DO
figure;
scatter(1:366,reg_clim_do_b + regm_do_b);
hold on
plot(1:366, yp_w_do_04_b,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('DO (mg/L)','fontsize',13)
title('정선관측(저층)-polynomial DO.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([2 18]) 
xlim([1 366])

%% salt
figure;
scatter(1:366,reg_clim_salt_b + regm_salt_b);
hold on
plot(1:366, yp_w_salt_04_b,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('salta (mg/m^3)','fontsize',13)
title('sumjin(songjung)-polynomial salt.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([32 35])
xlim([1 366])


%% no3
figure;
scatter(1:366,reg_clim_no3_b + regm_no3_b);
hold on
plot(1:366, yp_w_no3_04_b,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('NO3 (ug/L)','fontsize',13)
title('정선관측(저층)-polynomial no3.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([0 inf])
xlim([1 366])

%% temp
figure;
scatter(1:366,reg_clim_temp_b + regm_temp_b);
hold on
plot(1:366, yp_w_temp_04_b,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('temp (mg/L)','fontsize',13)
title('sumjin(songjung)-polynomial temp.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([0 inf])
xlim([1 366])

%% ALL
figure;
hold on;
% plot(1:366, yp_w_do_04_b,'--','color','b','linew',2)
% plot(1:366, yp_w_salt_04_b./1000,'--','color','g','linew',2)
plot(1:366, yp_w_no3_04,'--','color','b','linew',2)
plot(1:366, yp_w_no3_04_b,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('Concentration (ug/L)','fontsize',13)
title('정선관측(표층&저층)-polynomial.','fontsize',13)
grid on;
legend('no3-sur','no3-bot');
set(gca,'fontsize',13)
xlim([1 366])

%% plot raw
% make 1980~present
k=0
for i = 1980:2019
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end
t_tick_pre=sum(eom_d,2);

clearvars t_tick
for i = 1:length(1980:2019)-1
    t_tick(i)=sum(t_tick_pre(1:i))+1;
end
t_tick=[1; t_tick';];
tx_tick = t_tick;

%plot-NO3
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
gray = [128/255 128/255 128/255];
    clearvars temp_s temp_b
   temp_s= no3_sur;
   temp_b= no3_bot;
figure; hold on;
%  plot(temp_s,'*','color','b','linew',2);
 plot(temp_b,'*','color','r','linew',2);   
   
t=1:length(no3_sur);
%      temp_t_s(isnan(temp_t_s)) = interp1(t(~isnan(temp_t_s)), temp_t_s(~isnan(temp_t_s)), t(isnan(temp_t_s)) ); 
temp_s(isnan(temp_s)) = interp1( t(~isnan(temp_s)), temp_s(~isnan(temp_s)), t(isnan(temp_s)) ); 
temp_b(isnan(temp_b)) = interp1( t(~isnan(temp_b)), temp_b(~isnan(temp_b)), t(isnan(temp_b)) ); 

interp_nh4_b = temp_b;
interp_nh4_s = temp_s;

%  plot(temp_s,'color','b');
 plot(temp_b,'color','r');  
 ylim([0 inf])
xlabel('time(year)','fontsize',13)
ylabel('no3 (ug/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:3:end)]);
set(gca,'xlim',[tx_tick(8) tx_tick(end)]);
set(gca,'xticklabel',1980:3:2019);
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',tx_tick(1:1:end));
% legend('surface','bottom')


%plot-DO
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
gray = [128/255 128/255 128/255];
    clearvars temp_s temp_b
   temp_s= do_sur;
   temp_b= do_bot;
figure; hold on;
 plot(temp_s,'*','color','b','linew',2);
 plot(temp_b,'*','color','r','linew',2);
%  plot(temp_sur,'*','color','c','linew',2);
%  plot(temp_bot,'*','color','m','linew',2);   
   
t=1:length(no3_sur);
%      temp_t_s(isnan(temp_t_s)) = interp1(t(~isnan(temp_t_s)), temp_t_s(~isnan(temp_t_s)), t(isnan(temp_t_s)) ); 
temp_s(isnan(temp_s)) = interp1( t(~isnan(temp_s)), temp_s(~isnan(temp_s)), t(isnan(temp_s)) ); 
temp_b(isnan(temp_b)) = interp1( t(~isnan(temp_b)), temp_b(~isnan(temp_b)), t(isnan(temp_b)) ); 

interp_nh4_b = temp_b;
interp_nh4_s = temp_s;

 plot(temp_s,'color','b');
 plot(temp_b,'color','r');  
 ylim([0 inf])
xlabel('time(year)','fontsize',13)
ylabel('do (mg/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:1:11)]);
set(gca,'xlim',[tx_tick(1) tx_tick(11)]);
set(gca,'xticklabel',1980:1:1990);
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',tx_tick(1:1:end));
legend('surface','bottom')

%plot-temp
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
gray = [128/255 128/255 128/255];
    clearvars temp_s temp_b
   temp_s= temp_sur;
   temp_b= temp_bot;
figure; hold on;
 plot(temp_s,'*','color','b','linew',2);
 plot(temp_b,'*','color','r','linew',2);
%  plot(temp_sur,'*','color','c','linew',2);
%  plot(temp_bot,'*','color','m','linew',2);   
   
t=1:length(no3_sur);
%      temp_t_s(isnan(temp_t_s)) = interp1(t(~isnan(temp_t_s)), temp_t_s(~isnan(temp_t_s)), t(isnan(temp_t_s)) ); 
temp_s(isnan(temp_s)) = interp1( t(~isnan(temp_s)), temp_s(~isnan(temp_s)), t(isnan(temp_s)) ); 
temp_b(isnan(temp_b)) = interp1( t(~isnan(temp_b)), temp_b(~isnan(temp_b)), t(isnan(temp_b)) ); 

interp_nh4_b = temp_b;
interp_nh4_s = temp_s;

 plot(temp_s,'color','b');
 plot(temp_b,'color','r');  
 ylim([0 inf])
xlabel('time(year)','fontsize',13)
ylabel('temp (^oC)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:1:11)]);
set(gca,'xlim',[tx_tick(1) tx_tick(11)]);
set(gca,'xticklabel',1980:1:1990);
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',tx_tick(1:1:end));
legend('surface','bottom')


%plot- DO vs. temp
figure; hold on;
plot(temp_sur,do_sur,'.','color','b');
ylabel('DO (mg/L)','fontsize',13)
xlabel('Temp(^oC)','fontsize',13)
grid on; set(gca,'fontsize',13)

%regression
%slope y = b1*x
do_s_nonan=find(isnan(do_sur)==0);
do_b_nonan=find(isnan(do_bot)==0);

y_s_do_1 = polyfit(temp_sur(do_s_nonan), do_sur(do_s_nonan),1)
% re_do_sur = polyval(y_s_do_1, 5:0.01:30);
% plot(5:0.01:30,re_do_sur,'color','r');

yCalc = [5:0.01:30].* y_s_do_1(1) + y_s_do_1(2);
plot(5:0.01:30,yCalc,'-.','color','k','linew',2);
text(6,4.5,['y = ' num2str(y_s_do_1(1),'%.3f') 'x + ' num2str(y_s_do_1(2),'%.3f')],'fontsize',13);


%plot- DO vs. temp  bot
figure; hold on;
plot(temp_bot,do_bot,'.','color','b');
ylabel('DO (mg/L)','fontsize',13)
xlabel('Temp(^oC)','fontsize',13)
grid on; set(gca,'fontsize',13)

%regression
%slope y = b1*x
do_s_nonan=find(isnan(do_bot)==0);
do_b_nonan=find(isnan(do_bot)==0);
do_b_nonan(find(isnan(temp_bot(do_b_nonan)) == 1)) = [];

y_b_do_1 = polyfit(temp_bot(do_b_nonan), do_bot(do_b_nonan),1)
% re_do_bot = polyval(y_s_do_1, 5:0.01:30);
% plot(5:0.01:30,re_do_bot,'color','r');

yCalc_b = [0:0.01:30].* y_b_do_1(1) + y_b_do_1(2);
plot(0:0.01:30,yCalc_b,'-.','color','k','linew',2);
text(2,2,['y = ' num2str(y_b_do_1(1),'%.3f') 'x + ' num2str(y_b_do_1(2),'%.3f')],'fontsize',13);
plot(repmat([13],length(1:0.01:9),1),1:0.01:9,'-.','color','r','linew',2);

%plot-DO & recon.DO
% temp_s_nonan=find(isnan(temp_sur)==0);
% yCalc = temp_sur(temp_s_nonan) .* y_s_do_1(1) + y_s_do_1(2);

temp_b_nonan=find(isnan(temp_bot)==0);
yCalc_b = temp_bot(temp_b_nonan) .* y_b_do_1(1) + y_b_do_1(2);

clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
gray = [128/255 128/255 128/255];
    clearvars temp_s temp_b
   temp_s= do_sur;
   temp_b= do_bot;
   
figure;     
hold on;
%  plot(temp_s,'*','color','b','linew',2); 
% plot(temp_s_nonan, yCalc,'*','color','r','linew',2);
 plot(temp_b,'*','color','r','linew',2);
 plot(temp_b_nonan, yCalc_b,'*','color','b','linew',2);
    
t=1:length(no3_sur);
%      temp_t_s(isnan(temp_t_s)) = interp1(t(~isnan(temp_t_s)), temp_t_s(~isnan(temp_t_s)), t(isnan(temp_t_s)) ); 
temp_s(isnan(temp_s)) = interp1( t(~isnan(temp_s)), temp_s(~isnan(temp_s)), t(isnan(temp_s)) ); 
temp_b(isnan(temp_b)) = interp1( t(~isnan(temp_b)), temp_b(~isnan(temp_b)), t(isnan(temp_b)) ); 

interp_nh4_b = temp_b;
interp_nh4_s = temp_s;

%  plot(temp_s,'color','b');
 plot(temp_b,'color','r');  
 ylim([-inf inf])
xlabel('time(year)','fontsize',13)
ylabel('do (mg/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:1:11)]);
set(gca,'xlim',[tx_tick(1) tx_tick(11)]);
set(gca,'xticklabel',1980:1:1990);

% legend('surface','surface_r_e')
legend('bottom','bottom_r_e')

set(gca,'xtick',[tx_tick(1:5:end)]);
set(gca,'xlim',[tx_tick(1) tx_tick(end)]);
set(gca,'xticklabel',1980:5:2019);
set(gcf, 'Units', 'Inches', 'Position', [1, 1, 9, 7], ...
    'PaperUnits', 'Inches', 'PaperSize', [9, 7])

%plot- DO vs. DO  bot
figure; hold on;
plot(do_sur,do_bot,'.','color','b');
ylabel('DO (mg/L)','fontsize',13)
xlabel('DO-sur (mg/L)','fontsize',13)
grid on; set(gca,'fontsize',13)

%regression
%slope y = b1*x
do_s_nonan=find(isnan(do_bot)==0);
do_b_nonan=find(isnan(do_bot)==0);
do_b_nonan(find(isnan(do_sur(do_b_nonan)) == 1)) = [];

y_b_do_2 = polyfit(do_sur(do_b_nonan), do_bot(do_b_nonan),1)
% re_do_bot = polyval(y_s_do_1, 5:0.01:30);
% plot(5:0.01:30,re_do_bot,'color','r');

yCalc_b_2 = [4:0.01:8].* y_b_do_2(1) + y_b_do_2(2);
plot(4:0.01:8,yCalc_b_2,'-.','color','k','linew',2);
text(4.5,7.3,['y = ' num2str(y_b_do_2(1),'%.3f') 'x + ' num2str(y_b_do_2(2),'%.3f')],'fontsize',13);



%plot-DO & recon.DO  2
% temp_s_nonan, yCalc
do_re_2 = yCalc .* y_b_do_2(1) + y_b_do_2(2);

clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
gray = [128/255 128/255 128/255];
    clearvars temp_b
   temp_b= do_bot;
figure; hold on;
 plot(temp_b,'*','color','r','linew',2);
%  plot(temp_s_nonan, do_re_2,'*','color','b','linew',2);
    
t=1:length(no3_sur);
temp_b(isnan(temp_b)) = interp1( t(~isnan(temp_b)), temp_b(~isnan(temp_b)), t(isnan(temp_b)) ); 

 plot(temp_b,'color','r');  
 ylim([-inf inf])
xlabel('time(year)','fontsize',13)
ylabel('do (mg/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:1:11)]);
set(gca,'xlim',[tx_tick(1) tx_tick(11)]);
set(gca,'xticklabel',1980:1:1990);
legend('bottom','bottom_r_e')

set(gca,'xtick',[tx_tick(1:5:end)]);
set(gca,'xlim',[tx_tick(1) tx_tick(end)]);
set(gca,'xticklabel',1980:5:2019);

set(gcf, 'Units', 'Inches', 'Position', [1, 1, 9, 7], ...
    'PaperUnits', 'Inches', 'PaperSize', [9, 7])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%plot- DO vs. temps diff(surf-bot)
clearvars temp_diff*

below13=find(temp_bot < 13);
bt13=find(temp_bot >= 13);
temp_diff = temp_sur(bt13) - temp_bot(bt13);

% % make quasi_climate 1980~2019
% for i = 1:length(mmdd_indx)
%     if size(mmdd_indx{i},1) == 0     
%         reg_temp_diff{i} = {};
%     else
%         reg_temp_diff{i} = temp_diff(mmdd_indx{i})
%     end
% end
% 
% figure; hold on;
% for i = 1:366
%     plot(i,reg_temp_diff{i},'.','color','b');
% end
% xlabel('time(days)','fontsize',13)
% ylabel('temp diff(surf - bot) (^oC)','fontsize',13)
% title('temp diff(surf - bot)','fontsize',13)
% grid on
% set(gca,'fontsize',13)
% ylim([-inf inf])
% xlim([1 366])
% 
% xx = repmat([1:366],40,1);
% clearvars xx_mm
% for i = 1:366
% if i == 1
%     xx_mm = xx(:,i);  
% elseif i~= 1
%     xx_mm = [xx_mm(:,1); xx(:,i);];
% end
% end
% 
% 
% size_c = [];
% for i = 1:366; size_c(i)=length(reg_temp_diff{i}); end
% bar(size_c)
% 
% clearvars temp_diff_mtx
% for i = 1:366
%     if length(reg_temp_diff{i}) == 40
%        temp_diff_mtx(:,i) = [reg_temp_diff{i}];  
%     elseif length(reg_temp_diff{i}) ~= 40
%        temp_diff_mtx(:,i) = [reg_temp_diff{i}; NaN(40 - length(reg_temp_diff{i}),1);];
%     end
% end
% 
% clearvars temp_diff_fix
% for i = 1:366
% if i == 1
%     temp_diff_fix = temp_diff_mtx(:,i);  
% elseif i~= 1
%     temp_diff_fix = [temp_diff_fix(:,1); temp_diff_mtx(:,i);];
% end
% end
% 
% diff_f_nan = find(isnan(temp_diff_fix)==1);
% xx_mm(diff_f_nan) = [];
% temp_diff_fix(diff_f_nan) = [];
% 
% plot(xx_mm,temp_diff_fix,'.','color','b');
% xlabel('time(days)','fontsize',13)
% ylabel('temp diff(surf - bot) (^oC)','fontsize',13)
% title('temp diff(surf - bot)','fontsize',13)
% grid on
% set(gca,'fontsize',13)
% ylim([-inf inf])
% xlim([1 366])

%% bt13
figure; hold on;
plot(temp_diff, do_bot(bt13),'.','color','b');
ylabel('DO (mg/L)','fontsize',13)
xlabel('Temp(^oC)','fontsize',13)
grid on; set(gca,'fontsize',13)

%regression
%slope y = b1*x
bt13_re=bt13;
bt13_re(find(isnan(do_bot(bt13)) == 1)) = []; %  do_bot(temp_diff_nonan) 's NaN removing
temp_diff_re =  temp_sur(bt13_re) - temp_bot(bt13_re);

y_b_do_diff = polyfit(temp_diff_re, do_bot(bt13_re),1)
% re_do_sur = polyval(y_s_do_1, 5:0.01:30);
% plot(5:0.01:30,re_do_sur,'color','r');

yCalc = [-2:0.01:16].* y_b_do_diff(1) + y_b_do_diff(2);
plot(-2:0.01:16,yCalc,'-.','color','k','linew',2);
text(6,6,['y = ' num2str(y_b_do_diff(1),'%.3f') 'x + ' num2str(y_b_do_diff(2),'%.3f')],'fontsize',13);

%%below 13
figure; hold on;
plot(temp_bot(below13), do_bot(below13),'.','color','b');
ylabel('DO (mg/L)','fontsize',13)
xlabel('Temp(^oC)','fontsize',13)
grid on; set(gca,'fontsize',13)

%regression
%slope y = b1*x
below13_re=below13;
below13_re(find(isnan(do_bot(below13)) == 1)) = []; %  do_bot(temp_diff_nonan) 's NaN removing
temp_b_below = temp_bot(below13_re);

y_b_do_below = polyfit(temp_b_below, do_bot(below13_re),1)
% re_do_sur = polyval(y_s_do_1, 5:0.01:30);
% plot(5:0.01:30,re_do_sur,'color','r');

yCalc_below = [2:0.01:14].* y_b_do_below(1) + y_b_do_below(2);
plot(2:0.01:14,yCalc_below,'-.','color','k','linew',2);
text(9,4.5,['y = ' num2str(temp_b_below(1),'%.3f') 'x + ' num2str(temp_b_below(2),'%.3f')],'fontsize',13);


%% recon. diff T13
%plot-DO & recon.DO
clearvars temp_b yCalc_b_re yCal_con_int t t_in yCalc_below yCal_con_mer
% temp_b_nonan=find(isnan(temp_bot)==0);
yCalc_b_re = temp_diff .* y_b_do_diff(1) + y_b_do_diff(2);
yCalc_below = temp_b_below .* y_b_do_below(1) + y_b_do_below(2);

yCal_con_mer = [yCalc_b_re; yCalc_below;]; yCal_con_int = yCal_con_mer;
xx_con_mer = [bt13; below13_re;];

clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
gray = [128/255 128/255 128/255];
    clearvars temp_s temp_b
   temp_b= do_bot;
   
figure;     
hold on;
 plot(temp_b,'*','color','r','linew',2);
 plot(bt13, yCalc_b_re,'*','color','b','linew',2);
 plot(below13_re, yCalc_below, '*', 'color','c','linew',2);

 t=1:length(no3_sur);
%      temp_t_s(isnan(temp_t_s)) = interp1(t(~isnan(temp_t_s)), temp_t_s(~isnan(temp_t_s)), t(isnan(temp_t_s)) ); 
temp_b(isnan(temp_b)) = interp1( t(~isnan(temp_b)), temp_b(~isnan(temp_b)), t(isnan(temp_b)) ); 

[t t_in]=sort(xx_con_mer);
yCal_con_int = yCal_con_int(t_in);
yCal_con_int(isnan(yCal_con_int)) = interp1( t(~isnan(yCal_con_int)), yCal_con_int(~isnan(yCal_con_int)), t(isnan(yCal_con_int)) ); 

interp_nh4_b = temp_b;

 plot(temp_b,'color','r');
 plot(t,yCal_con_int,'color','g');  
 ylim([-inf inf])
xlabel('time(year)','fontsize',13)
ylabel('do (mg/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:1:11)]);
set(gca,'xlim',[tx_tick(1) tx_tick(11)]);
set(gca,'xticklabel',1980:1:1990);

% legend('surface','surface_r_e')
legend('bottom_O_B_S','bottom_r_e_-_u_p_1_3','bottom_r_e_-_d_o_w_n_1_3')

set(gca,'xtick',[tx_tick(1:5:end)]);
set(gca,'xlim',[tx_tick(1) tx_tick(end)]);
set(gca,'xticklabel',1980:5:2019);
set(gcf, 'Units', 'Inches', 'Position', [1, 1, 9, 7], ...
    'PaperUnits', 'Inches', 'PaperSize', [9, 7])







