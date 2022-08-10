close all; clear; clc;   % -v3

% % station name tag
% % (for preventing the miss picking using tag from 1 col's name on excel) 

% port
% st_selet = [5,13,14,15,16,17,22,24,25,26,27];
st_selet = [13,14,15,16,17,22,24,25,26,27]; %no num.5

for i = 1:length(st_selet)
name_tag_1{i} = [num2str(st_selet(i),'%02d'),'¹ø'] 
end

% combining the tag and outter point excluding
% name_tag = name_tag_1'; 
% name_tag{end+1:end+length(name_tag_2)} = name_tag_2; 
% name_tag{end+1:end+length(name_tag_3)} = name_tag_3; 
% size_tag = length(name_tag);

%% pick the row on the excel which has same name with tag
[raw_p txt_p]=xlsread('Á¤¼±ÇØ¾çÁ¶»ç_400line_from80.xls','sheet','');
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

clearvars ref_yy
k=0; m=0;
for i = 1:40
    l=0       
    for n = 1:12
        m = m+1;
        ref_yymm{m}=[num2str(i+1979) '-' num2str(n,'%02d')];
    for j = 1:eom_d(i,n)
        l=l+1; % index for 1~31(or 30 or 29)
        k=k+1; % index for 1~366
        ref_yymmdd{k}=[num2str(i+1979) '-' num2str(n,'%02d') '-'  num2str(j,'%02d')];
        ref_mmdd{k}=[num2str(n,'%02d') '-'  num2str(j,'%02d')];
        ref_yy{k}=[num2str(i+1979)];
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

% do_sur = do_sur(1,:)';
% do_bot = do_bot(1,:)';
% no3_sur =  no3_sur(1,:)'.* 14;  %umol/L -> ug/L
% no3_bot = no3_bot(1,:)' .* 14;
% temp_sur = temp_sur(1,:)';
% temp_bot = temp_bot(1,:)';
% salt_sur = salt_sur(1,:)';
% salt_bot = salt_bot(1,:)';

do_sur = do_sur';
do_bot = do_bot';
no3_sur =  no3_sur'.* 14;  %umol/L -> ug/L
no3_bot = no3_bot' .* 14;
temp_sur = temp_sur';
temp_bot = temp_bot';
salt_sur = salt_sur';
salt_bot = salt_bot';


% 2001~2010 data only no-nan to compare with yj-tak's output
%% matched 366
% make 366 mm-dd
% for i = 1:12
%   eom_d(i) = eomday(1980,i); % 1980 is leap-yr
% end

% k=0
% for i = 1:12
%     l=0
%     for j = 1:eom_d(i)
%         l=l+1; % index for 1~31(or 30 or 29)
%         k=k+1; % index for 1~366
%         mmdd(k,:)=[num2str(i,'%02d') '-'  num2str(l,'%02d')]
%     end
% end

% make it to cell-array
clearvars yy_c
for i = 1:40
    yy_c{i,1} = [num2str(1979+i,'%04d')]; % delete year
end

% pick matched date from water temp date
clearvars mmdd_indx yy_indx
for i = 1:length(yy_c)
       yy_indx{i} = find([strcmp(yy_c{i}, ref_yy)] == 1)
end

%confirm how may data on each the days
size_c = [];
for i = 1:40; size_c(i)=length(yy_indx{i}); end
bar(size_c)

% % make 1980~present
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

for i = 1:size(eom_d,1)
    for j = 1:size(eom_d,2)
        eom_d_each(i,j) = sum(eom_d(i,1:j));
    end
end

% make format to daily
for i = 1:22
    tx_tick(17+i) % if i=1, 2001 1/1 
    last_d_mon(i,:) = tx_tick(17+i)-1 + eom_d_each(17+i,:); % if i=1, 2001 : last day on each mon
    first_d_mon(i,1) = tx_tick(17+i);
    first_d_mon(i,2:12) = last_d_mon(i,1:11)+1;
end

% 1997:2018 has to be NaN; also 2011.01.01 ~ has to be NaN;
% % 1:tx_tick(18)-1  : ~1996.12.31 
% % tx_tick(39)+1:end  : 2019.01.01

do_sur(1:tx_tick(18)-1,:) = NaN; do_sur(tx_tick(39)+1:end,:) = NaN;
do_bot(1:tx_tick(18)-1,:) = NaN; do_bot(tx_tick(39)+1:end,:) = NaN;
no3_sur(1:tx_tick(18)-1,:) = NaN; no3_sur(tx_tick(39)+1:end,:) = NaN;
no3_bot(1:tx_tick(18)-1,:) = NaN; no3_bot(tx_tick(39)+1:end,:) = NaN;
temp_sur(1:tx_tick(18)-1,:) = NaN; temp_sur(tx_tick(39)+1:end,:) = NaN;
temp_bot(1:tx_tick(18)-1,:) = NaN; temp_bot(tx_tick(39)+1:end,:) = NaN;
salt_sur(1:tx_tick(18)-1,:) = NaN; salt_sur(tx_tick(39)+1:end,:) = NaN;
salt_bot(1:tx_tick(18)-1,:) = NaN; salt_bot(tx_tick(39)+1:end,:) = NaN;


%% extract over 3sig
for i = 1:length(name_tag_1)
    clearvars regime_*
    sig = 3; %% sigma
    
    % sur
    clearvars idx1 tempo_data
    tempo_data= do_sur(:,i); 
    idx1 = find(isnan(tempo_data) == 0);
    regime_do=tempo_data;
    regime_do(find(tempo_data > mean(tempo_data(idx1)) + sig*std(tempo_data(idx1))))=NaN;
    regime_do(find(tempo_data < mean(tempo_data(idx1)) - sig*std(tempo_data(idx1))))=NaN;
    mtx_regime_do(:,i) = regime_do;
    regm_do(:,i) = nanmean(regime_do)

    clearvars idx1 tempo_data
    tempo_data= no3_sur(:,i); 
    idx1 = find(isnan(tempo_data) == 0);
    regime_no3=tempo_data;
    regime_no3(find(tempo_data > mean(tempo_data(idx1)) + sig*std(tempo_data(idx1))))=NaN;
    regime_no3(find(tempo_data < mean(tempo_data(idx1)) - sig*std(tempo_data(idx1))))=NaN;
    mtx_regime_no3(:,i) = regime_no3;
    regm_no3(:,i) = nanmean(regime_no3)


    clearvars idx1 tempo_data
    tempo_data= temp_sur(:,i); 
    idx1 = find(isnan(tempo_data) == 0);
    regime_temp=tempo_data;
    regime_temp(find(tempo_data > mean(tempo_data(idx1)) + sig*std(tempo_data(idx1))))=NaN;
    regime_temp(find(tempo_data < mean(tempo_data(idx1)) - sig*std(tempo_data(idx1))))=NaN;
    mtx_regime_temp(:,i) = regime_temp;
    regm_temp(:,i) = nanmean(regime_temp)


    clearvars idx1 tempo_data
    tempo_data= salt_sur(:,i);
    idx1 = find(isnan(tempo_data) == 0);
    regime_salt=tempo_data;
    regime_salt(find(tempo_data > mean(tempo_data(idx1)) + sig*std(tempo_data(idx1))))=NaN;
    regime_salt(find(tempo_data < mean(tempo_data(idx1)) - sig*std(tempo_data(idx1))))=NaN;
    mtx_regime_salt(:,i) = regime_salt;
    regm_salt(:,i) = nanmean(regime_salt)

    %bot
    clearvars idx1 tempo_data
    tempo_data= do_bot(:,i);
    idx1 = find(isnan(tempo_data) == 0);
    regime_do_b=tempo_data;
    regime_do_b(find(tempo_data > mean(tempo_data(idx1)) + sig*std(tempo_data(idx1))))=NaN;
    regime_do_b(find(tempo_data < mean(tempo_data(idx1)) - sig*std(tempo_data(idx1))))=NaN;
    mtx_regime_do_b(:,i) = regime_do_b;
    regm_do_b(:,i) = nanmean(regime_do_b)

    clearvars idx1 tempo_data
    tempo_data= no3_bot(:,i);
    idx1 = find(isnan(tempo_data) == 0);
    regime_no3_b=tempo_data;
    regime_no3_b(find(tempo_data > mean(tempo_data(idx1)) + sig*std(tempo_data(idx1))))=NaN;
    regime_no3_b(find(tempo_data < mean(tempo_data(idx1)) - sig*std(tempo_data(idx1))))=NaN;
    mtx_regime_no3_b(:,i) = regime_no3_b;
    regm_no3_b(:,i) = nanmean(regime_no3_b)


    clearvars idx1 tempo_data
    tempo_data= temp_bot(:,i);
    idx1 = find(isnan(tempo_data) == 0);
    regime_temp_b=tempo_data;
    regime_temp_b(find(tempo_data > mean(tempo_data(idx1)) + sig*std(tempo_data(idx1))))=NaN;
    regime_temp_b(find(tempo_data < mean(tempo_data(idx1)) - sig*std(tempo_data(idx1))))=NaN;
    mtx_regime_temp_b(:,i) = regime_temp_b;
    regm_temp_b(:,i) = nanmean(regime_temp)


    clearvars idx1 tempo_data
    tempo_data= salt_bot(:,i);
    idx1 = find(isnan(tempo_data) == 0);
    regime_salt_b=tempo_data;
    regime_salt_b(find(tempo_data > mean(tempo_data(idx1)) + sig*std(tempo_data(idx1))))=NaN;
    regime_salt_b(find(tempo_data < mean(tempo_data(idx1)) - sig*std(tempo_data(idx1))))=NaN;
    mtx_regime_salt_b(:,i) = regime_salt_b;
    regm_salt_b(:,i) = nanmean(regime_salt)

%% extract regime mean
% regime_do = regime_do - regm_do;
% regime_no3 = regime_no3 - regm_no3;
% regime_temp = regime_temp - regm_temp;
% regime_salt = regime_salt - regm_salt;
% 
% regime_do_b = regime_do_b - regm_do_b;
% regime_no3_b = regime_no3_b - regm_no3_b;
% regime_temp_b = regime_temp_b - regm_temp_b;
% regime_salt_b = regime_salt_b - regm_salt_b;
end

clearvars regime_*
regime_do = mtx_regime_do;
regime_no3 = mtx_regime_no3;
regime_temp = mtx_regime_temp;
regime_salt = mtx_regime_salt;

regime_do_b = mtx_regime_do_b;
regime_no3_b = mtx_regime_no3_b;
regime_temp_b = mtx_regime_temp_b;
regime_salt_b = mtx_regime_salt_b;


%% matched 366

% date also *5 dim

% % make 366 mm-dd
% for i = 1:12
%   eom_d(i) = eomday(1980,i); % 1980 is leap-yr
% end
% 
% k=0
% for i = 1:12
%     l=0
%     for j = 1:eom_d(i)
%         l=l+1; % index for 1~31(or 30 or 29)
%         k=k+1; % index for 1~366
%         mmdd(k,:)=[num2str(i,'%02d') '-'  num2str(l,'%02d')]
%     end
% end
% 
% % make it to cell-array
% for i = 1:length(mmdd)
%     mmdd_c{i,1} = mmdd(i,:); % delete year
% end
% 
% % pick matched date from water temp date
% for i = 1:length(mmdd_c)
%        mmdd_indx{i} = find([strcmp(mmdd_c{i}, ref_mmdd)] == 1)
% end
% 
% 
% %confirm how may data on each the days
% size_c = [];
% for i = 1:366; size_c(i)=length(mmdd_indx{i}); end
% bar(size_c)

% make quasi_climate 1980~2019 

% diminish dim.
clearvars reg_clim_*
for j = 1:length(name_tag_1)
for i = 1:length(yy_indx)
clearvars regime_*_cut
    regime_do_cut = regime_do(:,j);
    regime_no3_cut = regime_no3(:,j);
    regime_temp_cut = regime_temp(:,j);
    regime_salt_cut = regime_salt(:,j);
    regime_do_b_cut = regime_do_b(:,j);
    regime_no3_b_cut = regime_no3_b(:,j);
    regime_temp_b_cut = regime_temp_b(:,j);
    regime_salt_b_cut = regime_salt_b(:,j);
    
    if size(yy_indx{i},1) == 0     
        reg_clim_do(j,i) = NaN;
        reg_clim_no3(j,i) = NaN;
        reg_clim_temp(j,i) = NaN; 
        reg_clim_salt(j,i) = NaN;
        reg_clim_do_b(j,i) = NaN;
        reg_clim_no3_b(j,i) = NaN;
        reg_clim_temp_b(j,i) = NaN; 
        reg_clim_salt_b(j,i) = NaN;
    else
        reg_clim_do(j,i) = nanmean(regime_do_cut(yy_indx{i}));
        reg_clim_no3(j,i) = nanmean(regime_no3_cut(yy_indx{i}));
        reg_clim_temp(j,i) = nanmean(regime_temp_cut(yy_indx{i}));
        reg_clim_salt(j,i) = nanmean(regime_salt_cut(yy_indx{i}));
        reg_clim_do_b(j,i) = nanmean(regime_do_b_cut(yy_indx{i}));
        reg_clim_no3_b(j,i) = nanmean(regime_no3_b_cut(yy_indx{i}));
        reg_clim_temp_b(j,i) = nanmean(regime_temp_b_cut(yy_indx{i}));
        reg_clim_salt_b(j,i) = nanmean(regime_salt_b_cut(yy_indx{i}));
    end
end
end

%% DO
% figure;
% scatter(1:366,reg_clim_do);
% hold on
% xlabel('time(days)','fontsize',13)
% ylabel('DO (mg/L)','fontsize',13)
% title('Á¤¼±°üÃø(Ç¥Ãþ)-polynomial DO.','fontsize',13)
% grid on
% set(gca,'fontsize',13)
% ylim([2 10]) 
% xlim([1 366])

%% salt
color_pic = lines(size(reg_clim_salt,1));
marker_sty = {'o','+','x','^','>','h','p','s','d','.','*','v','<'};
figure;
for i = 1:size(reg_clim_salt,1)
% scatter(1:40,reg_clim_salt(i,:),color_pic{i});
scatter(1:40,reg_clim_salt(i,:),marker_sty{i},'MarkerEdgeColor',color_pic(i,:));
hold on
end
xlabel('time(year)','fontsize',13)
ylabel('salt (mg/m^3)','fontsize',13)
set(gca,'xtick',[1:5:40]);
set(gca,'xlim',[1 40]);
set(gca,'xticklabel',1980:5:2019);
title('Á¤¼±°üÃø-Ç¥Ãþ¿°ºÐ ¿¬Æò±Õ','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([32 35])
legend('400-13','400-14','400-15','400-16','400-17','400-22','400-24','400-25','400-26','400-27')
% legend('205-01','205-02','205-03','205-04','205-05')

return

%% no3
figure;
for i = 1:size(reg_clim_salt,1)
scatter(1:40,reg_clim_no3(i,:),color_pic{i});
hold on
end
xlabel('time(year)','fontsize',13)
ylabel('NO3 (ug/L)','fontsize',13)
set(gca,'xtick',[1:5:40]);
set(gca,'xlim',[1 40]);
set(gca,'xticklabel',1980:5:2019);
title('Á¤¼±°üÃø-Ç¥ÃþÁú»ê¿° ¿¬Æò±Õ','fontsize',13)
grid on
set(gca,'fontsize',13)
% ylim([32 35])
% legend('205-01','205-02','205-03','205-04','205-05')

%% temp
figure;
for i = 1:size(reg_clim_salt,1)
scatter(1:40,reg_clim_temp(i,:),color_pic{i});
hold on
end
xlabel('time(year)','fontsize',13)
ylabel('temperature (^oC)','fontsize',13)
set(gca,'xtick',[1:5:40]);
set(gca,'xlim',[1 40]);
set(gca,'xticklabel',1980:5:2019);
title('Á¤¼±°üÃø-Ç¥Ãþ¼ö¿Â ¿¬Æò±Õ','fontsize',13)
grid on
set(gca,'fontsize',13)



%bottom
%% DO
figure;
scatter(1:366,reg_clim_do_b);
hold on
xlabel('time(days)','fontsize',13)
ylabel('DO (mg/L)','fontsize',13)
title('Á¤¼±°üÃø(ÀúÃþ)-polynomial DO.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([2 10]) 
xlim([1 366])

%% salt
figure;
for i = 1:size(reg_clim_salt,1)
scatter(1:40,reg_clim_salt_b(i,:),color_pic{i});
hold on
end
xlabel('time(year)','fontsize',13)
ylabel('salt (mg/m^3)','fontsize',13)
set(gca,'xtick',[1:5:40]);
set(gca,'xlim',[1 40]);
set(gca,'xticklabel',1980:5:2019);
title('Á¤¼±°üÃø-¹Ù´ÚÃþ¿°ºÐ ¿¬Æò±Õ','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([32 35])
% legend('205-01','205-02','205-03','205-04','205-05')

%% no3
figure;
for i = 1:size(reg_clim_salt,1)
scatter(1:40,reg_clim_no3_b(i,:),color_pic{i});
hold on
end
xlabel('time(year)','fontsize',13)
ylabel('NO3 (ug/L)','fontsize',13)
set(gca,'xtick',[1:5:40]);
set(gca,'xlim',[1 40]);
set(gca,'xticklabel',1980:5:2019);
title('Á¤¼±°üÃø-¹Ù´ÚÃþÁú»ê¿° ¿¬Æò±Õ','fontsize',13)
grid on
set(gca,'fontsize',13)
% ylim([32 35])
% legend('205-01','205-02','205-03','205-04','205-05')

%% temp
figure;
for i = 1:size(reg_clim_salt,1)
scatter(1:40,reg_clim_temp_b(i,:),color_pic{i});
hold on
end
xlabel('time(year)','fontsize',13)
ylabel('temperature (^oC)','fontsize',13)
set(gca,'xtick',[1:5:40]);
set(gca,'xlim',[1 40]);
set(gca,'xticklabel',1980:5:2019);
title('Á¤¼±°üÃø-¹Ù´ÚÃþ¼ö¿Â ¿¬Æò±Õ','fontsize',13)
grid on
set(gca,'fontsize',13)


% right_here
%% linear trend coeff. and save file

%% salt
color_pic = lines(size(reg_clim_salt,1));
marker_sty = {'o','+','x','^','>','h','p','s','d','.','*','v','<'};
xp = 1:40;
figure; hold on;
for i = 1:size(reg_clim_salt,1)
  clearvars reg_data_salt xp_w_salt pf_w_salt
reg_data_salt = reg_clim_salt(i,:);
xp_w_salt = find(isnan(reg_data_salt)==0);
pf_w_salt = polyfit(xp_w_salt,reg_data_salt(xp_w_salt),1);
yp_w_salt(i,:) = polyval(pf_w_salt,xp);
scatter(1:40,reg_clim_salt(i,:),marker_sty{i},'MarkerEdgeColor',color_pic(i,:));
plot(1:40, yp_w_salt(i,:),'color',color_pic(i,:));
coeff_salt(i,:) = pf_w_salt;
hold on
end
xlabel('time(year)','fontsize',13)
ylabel('salt (mg/m^3)','fontsize',13)
set(gca,'xtick',[1:5:40]);
set(gca,'xlim',[1 40]);
set(gca,'xticklabel',1980:5:2019);
title('Á¤¼±°üÃø-Ç¥Ãþ¿°ºÐ ¿¬Æò±Õ','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([32 35])


%% no3
figure; hold on;
for i = 1:size(reg_clim_salt,1)
clearvars reg_data_no3 xp_w_no3 pf_w_no3
if i ~= 1 && i ~= 3 && i ~=5 && i ~= 6 && i ~= 7 && i ~= 9
reg_data_no3 = reg_clim_no3(i,:);
xp_w_no3 = find(isnan(reg_data_no3)==0);
pf_w_no3 = polyfit(xp_w_no3,reg_data_no3(xp_w_no3),1);
yp_w_no3(i,:) = polyval(pf_w_no3,xp);
scatter(1:40,reg_clim_no3(i,:),marker_sty{i},'MarkerEdgeColor',color_pic(i,:));
plot(1:40, yp_w_no3(i,:),'color',color_pic(i,:));
coeff_no3(i,:) = pf_w_no3;
hold on
end
end
xlabel('time(year)','fontsize',13)
ylabel('NO3 (ug/L)','fontsize',13)
set(gca,'xtick',[1:5:40]);
set(gca,'xlim',[1 40]);
set(gca,'xticklabel',1980:5:2019);
title('Á¤¼±°üÃø-Ç¥ÃþÁú»ê¿° ¿¬Æò±Õ','fontsize',13)
grid on
set(gca,'fontsize',13)
% ylim([32 35])
% legend('205-01','205-02','205-03','205-04','205-05')

%temp
xp = 1:40;
figure; hold on;
for i = 1:size(reg_clim_salt,1)
  clearvars reg_data_temp xp_w_temp pf_w_temp
reg_data_temp = reg_clim_temp(i,:);
xp_w_temp = find(isnan(reg_data_temp)==0);
pf_w_temp = polyfit(xp_w_temp,reg_data_temp(xp_w_temp),1);
yp_w_temp(i,:) = polyval(pf_w_temp,xp);
scatter(1:40,reg_clim_temp(i,:),marker_sty{i},'MarkerEdgeColor',color_pic(i,:));
plot(1:40, yp_w_temp(i,:),'color',color_pic(i,:));
coeff_temp(i,:) = pf_w_temp;
hold on
end
xlabel('time(year)','fontsize',13)
ylabel('temperature (^oC)','fontsize',13)
set(gca,'xtick',[1:5:40]);
set(gca,'xlim',[1 40]);
set(gca,'xticklabel',1980:5:2019);
title('Á¤¼±°üÃø-Ç¥Ãþ¼ö¿Â ¿¬Æò±Õ','fontsize',13)
grid on
set(gca,'fontsize',13)


lon_400(1) = 128.5;
lat_400(1) = 34.0767;
lon_400(2) = 128.4;
lat_400(2) = 34.225;
lon_400(3) = 128.3067;
lat_400(3) = 34.37;
lon_400(4) = 128.2033;
lat_400(4) = 34.5167;
lon_400(5) = 128.0083;
lat_400(5) = 34.45;
lon_400(6) = 128.3983;
lat_400(6) = 33.7383;
lon_400(7) = 127.865;
lat_400(7) = 33.5933;
lon_400(8) =127.5667;
lat_400(8) = 33.56;
lon_400(9) = 127.28;
lat_400(9) = 33.5283;
lon_400(10) = 127.075;
lat_400(10) = 33.505;

save('kodc_400_1997_2018_just_mean_yearly.mat','-v7.3');







