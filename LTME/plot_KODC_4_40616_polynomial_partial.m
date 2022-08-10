close all; clear; clc;   % -v3

% % station name tag
% % (for preventing the miss picking using tag from 1 col's name on excel) 

% port
name_tag_1{1} = [num2str(16,'%02d'),'¹ø'] 
 
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

po4_nan_p = strcmp(txt_cut(:,12),''); % detect NaN
for i = 1:length(txt_cut)
   if po4_nan_p(i) == 1 % when it's nan
       po4_p(i) = NaN;
   else
       po4_p(i) = str2num(char(txt_cut(i,12)));
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
indx_st={};
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
    po4_kodc{i} = po4_p(indx_st{i});
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


%po4
for i = 1:size(indx_date_s,1)
    clearvars temp_kodc_k
    temp_kodc_k = po4_kodc{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_kodc_k(indx_date_s{i,j}))) ~= 0     
        po4_sur(i,j) = nanmean(temp_kodc_k(indx_date_s{i,j}));
    elseif sum(size(temp_kodc_k(indx_date_s{i,j}))) == 0 
        po4_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_kodc_k
    temp_kodc_k = po4_kodc{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_kodc_k(indx_date_b{i,j}))) ~= 0     
        po4_bot(i,j) = nanmean(temp_kodc_k(indx_date_b{i,j}));
    elseif sum(size(temp_kodc_k(indx_date_b{i,j}))) == 0 
        po4_bot(i,j) = NaN;
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
no3_sur =  no3_sur(1,:)'.* 14;   %umol N /L -> ug N /L (wrong unit on the excel) on the web it is umol/L
no3_bot = no3_bot(1,:)' .* 14;
temp_sur = temp_sur(1,:)';
temp_bot = temp_bot(1,:)';
salt_sur = salt_sur(1,:)';
salt_bot = salt_bot(1,:)';
po4_sur = po4_sur(1,:)'.* 30.973762; %umol P /L -> ug P /L
po4_bot = po4_bot(1,:)'.* 30.973762; %umol P /L -> ug P /L



%% divide the regime (in case, 2004)

% no3-kodc 40016   s: 2002-05-01, 2010-04-01 | b : 2000-11-01, 2010-04-01

%divide the data
% 1980 ~ 2000.10.31
t_regime_1 = 8187; % 2002-05-31
t_regime_b_1 = 7640; %  2000-11-30,

regime_no3_s_pre_1 = no3_sur(1:t_regime_1);
regime_no3_b_pre_1 = no3_bot(1:t_regime_b_1);

% 2000.11.01 ~ 2010.04.30
t_regime_2 = 11078; % 2010.04.30
t_regime_b_2 = t_regime_2; % 2010.04.30

regime_no3_s_pre_2 = no3_sur(t_regime_1+1:t_regime_2);
regime_no3_b_pre_2 = no3_bot(t_regime_b_1+1:t_regime_b_2);

% 2010.05.01 ~ 2019.12.31
regime_no3_s_pre_3 = no3_sur(t_regime_2+1:end);
regime_no3_b_pre_3 = no3_bot(t_regime_b_2+1:end);


% 1980 ~ 2011-03 (po4 bot)
% ref_yymmdd{11413} is 2011-03
t_regime_po4_b_1 = 11413;
t_regime_po4_b_2 = 11414;
regime_po4_b_pre_1 = po4_bot(1:t_regime_po4_b_1);
regime_po4_b_pre_2 = po4_bot(t_regime_po4_b_2:end);



%% extract over 3sig
sig = 3; %% sigma
 
%% 1st
% sur
clearvars idx1
idx1 = find(isnan(regime_no3_s_pre_1) == 0);
regime_no3_1=regime_no3_s_pre_1;
regime_no3_1(find(regime_no3_s_pre_1 > mean(regime_no3_s_pre_1(idx1)) + sig*std(regime_no3_s_pre_1(idx1))))=NaN;
regime_no3_1(find(regime_no3_s_pre_1 < mean(regime_no3_s_pre_1(idx1)) - sig*std(regime_no3_s_pre_1(idx1))))=NaN;
regm_no3_1 = nanmean(regime_no3_1)

%bot
clearvars idx1
idx1 = find(isnan(regime_no3_b_pre_1) == 0);
regime_no3_b_1=regime_no3_b_pre_1;
regime_no3_b_1(find(regime_no3_b_pre_1 > mean(regime_no3_b_pre_1(idx1)) + sig*std(regime_no3_b_pre_1(idx1))))=NaN;
regime_no3_b_1(find(regime_no3_b_pre_1 < mean(regime_no3_b_pre_1(idx1)) - sig*std(regime_no3_b_pre_1(idx1))))=NaN;
regm_no3_b_1 = nanmean(regime_no3_b_1)

%% 2nd
clearvars idx1
idx1 = find(isnan(regime_no3_s_pre_2) == 0);
regime_no3_2=regime_no3_s_pre_2;
regime_no3_2(find(regime_no3_s_pre_2 > mean(regime_no3_s_pre_2(idx1)) + sig*std(regime_no3_s_pre_2(idx1))))=NaN;
regime_no3_2(find(regime_no3_s_pre_2 < mean(regime_no3_s_pre_2(idx1)) - sig*std(regime_no3_s_pre_2(idx1))))=NaN;
regm_no3_2 = nanmean(regime_no3_2)

%bot
clearvars idx1
idx1 = find(isnan(regime_no3_b_pre_2) == 0);
regime_no3_b_2=regime_no3_b_pre_2;
regime_no3_b_2(find(regime_no3_b_pre_2 > mean(regime_no3_b_pre_2(idx1)) + sig*std(regime_no3_b_pre_2(idx1))))=NaN;
regime_no3_b_2(find(regime_no3_b_pre_2 < mean(regime_no3_b_pre_2(idx1)) - sig*std(regime_no3_b_pre_2(idx1))))=NaN;
regm_no3_b_2 = nanmean(regime_no3_b_2)


%% 3rd
clearvars idx1
idx1 = find(isnan(regime_no3_s_pre_3) == 0);
regime_no3_3=regime_no3_s_pre_3;
regime_no3_3(find(regime_no3_s_pre_3 > mean(regime_no3_s_pre_3(idx1)) + sig*std(regime_no3_s_pre_3(idx1))))=NaN;
regime_no3_3(find(regime_no3_s_pre_3 < mean(regime_no3_s_pre_3(idx1)) - sig*std(regime_no3_s_pre_3(idx1))))=NaN;
regm_no3_3 = nanmean(regime_no3_3)

%bot
clearvars idx1
idx1 = find(isnan(regime_no3_b_pre_3) == 0);
regime_no3_b_3=regime_no3_b_pre_3;
regime_no3_b_3(find(regime_no3_b_pre_3 > mean(regime_no3_b_pre_3(idx1)) + sig*std(regime_no3_b_pre_3(idx1))))=NaN;
regime_no3_b_3(find(regime_no3_b_pre_3 < mean(regime_no3_b_pre_3(idx1)) - sig*std(regime_no3_b_pre_3(idx1))))=NaN;
regm_no3_b_3 = nanmean(regime_no3_b_3)

regime_no3_1_raw = regime_no3_1;
regime_no3_b_1_raw = regime_no3_b_1;
regime_no3_2_raw = regime_no3_2;
regime_no3_b_2_raw = regime_no3_b_2;
regime_no3_3_raw = regime_no3_3;
regime_no3_b_3_raw = regime_no3_b_3;

%% po4
% surface has no-regime
clearvars idx1
idx1 = find(isnan(po4_sur) == 0);
regime_po4=po4_sur;
regime_po4(find(po4_sur > mean(po4_sur(idx1)) + sig*std(po4_sur(idx1))))=NaN;
regime_po4(find(po4_sur < mean(po4_sur(idx1)) - sig*std(po4_sur(idx1))))=NaN;
regm_po4 = nanmean(regime_po4)

%bot
clearvars idx1
idx1 = find(isnan(regime_po4_b_pre_1) == 0);
regime_po4_b_1=regime_po4_b_pre_1;
regime_po4_b_1(find(regime_po4_b_pre_1 > mean(regime_po4_b_pre_1(idx1)) + sig*std(regime_po4_b_pre_1(idx1))))=NaN;
regime_po4_b_1(find(regime_po4_b_pre_1 < mean(regime_po4_b_pre_1(idx1)) - sig*std(regime_po4_b_pre_1(idx1))))=NaN;
regm_po4_b_1 = nanmean(regime_po4_b_1)

%% 2nd

%bot
clearvars idx1
idx1 = find(isnan(regime_po4_b_pre_2) == 0);
regime_po4_b_2=regime_po4_b_pre_2;
regime_po4_b_2(find(regime_po4_b_pre_2 > mean(regime_po4_b_pre_2(idx1)) + sig*std(regime_po4_b_pre_2(idx1))))=NaN;
regime_po4_b_2(find(regime_po4_b_pre_2 < mean(regime_po4_b_pre_2(idx1)) - sig*std(regime_po4_b_pre_2(idx1))))=NaN;
regm_po4_b_2 = nanmean(regime_po4_b_2)


%% extract regime mean
% 1st
regime_no3_1 = regime_no3_1 - regm_no3_1;
regime_no3_b_1 = regime_no3_b_1 - regm_no3_b_1;
% 2nd
regime_no3_2 = regime_no3_2 - (regm_no3_2 - regm_no3_1) - regm_no3_1;
regime_no3_b_2 = regime_no3_b_2 - (regm_no3_b_2 - regm_no3_b_1) - regm_no3_b_1;

regime_no3_only_2 = regime_no3_2_raw - regm_no3_2 ;
regime_no3_b_only_2 = regime_no3_b_2_raw - regm_no3_b_2;

% 3rd
regime_no3_3 = regime_no3_3 - (regm_no3_3 - regm_no3_1) - regm_no3_1;
regime_no3_b_3 = regime_no3_b_3 - (regm_no3_b_3 - regm_no3_b_1) - regm_no3_b_1;

regime_no3_only_3 = regime_no3_3_raw - regm_no3_3 ;
regime_no3_b_only_3 = regime_no3_b_3_raw - regm_no3_b_3;

%__ po4 ___________________________________________________________________
%sur 
regime_po4 = regime_po4 - regm_po4;
%bot
regime_po4_b_1 = regime_po4_b_1 - regm_po4_b_1;
regime_po4_b_only_2 = regime_po4_b_2 - regm_po4_b_2;

regime_po4_b_2 = regime_po4_b_2 - (regm_po4_b_2 - regm_po4_b_1) - regm_po4_b_1;

%__________________________________________________________________________

% merge matrix

regime_no3 = [regime_no3_1; regime_no3_2; regime_no3_3;];
regime_no3_b = [regime_no3_b_1; regime_no3_b_2; regime_no3_b_3;];
regime_po4_b = [regime_po4_b_1; regime_po4_b_2;];

regime_no3_on_1 = [regime_no3_1; NaN(length(regime_no3_2),1); NaN(length(regime_no3_3),1);];
% only 2, 3 's effect can be shown 
regime_no3_on_2 = [NaN(length(regime_no3_1),1); regime_no3_only_2; NaN(length(regime_no3_3),1);];
regime_no3_on_3 = [NaN(length(regime_no3_1),1); NaN(length(regime_no3_2),1); regime_no3_only_3;];



regime_no3_b_on_1 = [regime_no3_b_1; NaN(length(regime_no3_b_2),1); NaN(length(regime_no3_b_3),1);];
regime_no3_b_on_2 = [NaN(length(regime_no3_b_1),1); regime_no3_b_only_2; NaN(length(regime_no3_b_3),1);];
regime_no3_b_on_3 = [NaN(length(regime_no3_b_1),1); NaN(length(regime_no3_b_2),1); regime_no3_b_only_3;];

regime_po4_b_on_1 = [regime_po4_b_1; NaN(length(regime_po4_b_2),1);];
regime_po4_b_on_2 = [NaN(length(regime_po4_b_1),1); regime_po4_b_2;];
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
clearvars reg_clim_no3*
for i = 1:length(mmdd_indx)
    if size(mmdd_indx{i},1) == 0     
        reg_clim_no3(i) = NaN;
        reg_clim_no3_b(i) = NaN;
        reg_clim_no3_1(i) = NaN;
        reg_clim_no3_b_1(i) = NaN;
        reg_clim_no3_2(i) = NaN;
        reg_clim_no3_b_2(i) = NaN;
        reg_clim_no3_3(i) = NaN;
        reg_clim_no3_b_3(i) = NaN;
    else
        reg_clim_no3(i) = nanmean(regime_no3(mmdd_indx{i}));
        reg_clim_no3_b(i) = nanmean(regime_no3_b(mmdd_indx{i}));
        reg_clim_no3_1(i) = nanmean(regime_no3_on_1(mmdd_indx{i}));
        reg_clim_no3_b_1(i) = nanmean(regime_no3_b_on_1(mmdd_indx{i}));
        reg_clim_no3_2(i) = nanmean(regime_no3_on_2(mmdd_indx{i}));
        reg_clim_no3_b_2(i) = nanmean(regime_no3_b_on_2(mmdd_indx{i}));
        reg_clim_no3_3(i) = nanmean(regime_no3_on_3(mmdd_indx{i}));
        reg_clim_no3_b_3(i) = nanmean(regime_no3_b_on_3(mmdd_indx{i}));
        reg_clim_po4(i) = nanmean(regime_po4(mmdd_indx{i}));
        reg_clim_po4_b(i) = nanmean(regime_po4_b(mmdd_indx{i}));
        reg_clim_po4_b_1(i) = nanmean(regime_po4_b_on_1(mmdd_indx{i}));
        reg_clim_po4_b_2(i) = nanmean(regime_po4_b_on_2(mmdd_indx{i}));
    end
end

%% full
%surf
%NO3
clearvars xp_w_no3 pf_w_no3
xp = 1:366;
xp_w_no3 = find(isnan(reg_clim_no3)==0);
pf_w_no3 = polyfit(xp_w_no3,reg_clim_no3(xp_w_no3),3);
yp_w_no3_04 = polyval(pf_w_no3,xp);
yp_w_no3_04 = yp_w_no3_04 + regm_no3_1;

%bottom
%NO3
clearvars xp_w_no3_b pf_w_no3_b
xp = 1:366;
xp_w_no3_b = find(isnan(reg_clim_no3_b)==0);
pf_w_no3_b = polyfit(xp_w_no3_b,reg_clim_no3_b(xp_w_no3_b),3);
yp_w_no3_04_b = polyval(pf_w_no3_b,xp);
yp_w_no3_04_b = yp_w_no3_04_b + regm_no3_b_1;

%% 1st
%surf
%NO3
clearvars xp_w_no3 pf_w_no3
xp = 1:366;
xp_w_no3 = find(isnan(reg_clim_no3_1)==0);
pf_w_no3 = polyfit(xp_w_no3,reg_clim_no3_1(xp_w_no3),3);
yp_w_no3_1 = polyval(pf_w_no3,xp);
yp_w_no3_1 = yp_w_no3_1 + regm_no3_1;

%bottom
%NO3
clearvars xp_w_no3_b pf_w_no3_b
xp = 1:366;
xp_w_no3_b = find(isnan(reg_clim_no3_b_1)==0);
pf_w_no3_b = polyfit(xp_w_no3_b,reg_clim_no3_b_1(xp_w_no3_b),3);
yp_w_no3_b_1 = polyval(pf_w_no3_b,xp);
yp_w_no3_b_1 = yp_w_no3_b_1 + regm_no3_b_1;

%% 2nd
%surf
%NO3
clearvars xp_w_no3 pf_w_no3
xp = 1:366;
xp_w_no3 = find(isnan(reg_clim_no3_2)==0);
pf_w_no3 = polyfit(xp_w_no3,reg_clim_no3_2(xp_w_no3),3);
yp_w_no3_2 = polyval(pf_w_no3,xp);
yp_w_no3_2 = yp_w_no3_2 + regm_no3_2;

%bottom
%NO3
clearvars xp_w_no3_b pf_w_no3_b
xp = 1:366;
xp_w_no3_b = find(isnan(reg_clim_no3_b_2)==0);
pf_w_no3_b = polyfit(xp_w_no3_b,reg_clim_no3_b_2(xp_w_no3_b),3);
yp_w_no3_b_2 = polyval(pf_w_no3_b,xp);
yp_w_no3_b_2 = yp_w_no3_b_2 + regm_no3_b_2;

%% 3rd
%surf
%NO3
clearvars xp_w_no3 pf_w_no3
xp = 1:366;
xp_w_no3 = find(isnan(reg_clim_no3_3)==0);
pf_w_no3 = polyfit(xp_w_no3,reg_clim_no3_3(xp_w_no3),3);
yp_w_no3_3 = polyval(pf_w_no3,xp);
yp_w_no3_3 = yp_w_no3_3 + regm_no3_3;

%bottom
%NO3
clearvars xp_w_no3_b pf_w_no3_b
xp = 1:366;
xp_w_no3_b = find(isnan(reg_clim_no3_b_3)==0);
pf_w_no3_b = polyfit(xp_w_no3_b,reg_clim_no3_b_3(xp_w_no3_b),3);
yp_w_no3_b_3 = polyval(pf_w_no3_b,xp);
yp_w_no3_b_3 = yp_w_no3_b_3 + regm_no3_b_3;

%surf
%po4
clearvars xp_w_po4 pf_w_po4
xp = 1:366;
xp_w_po4 = find(isnan(reg_clim_po4)==0);
pf_w_po4 = polyfit(xp_w_po4,reg_clim_po4(xp_w_po4),3);
yp_w_po4_04 = polyval(pf_w_po4,xp);
yp_w_po4= yp_w_po4_04 + regm_po4;

%bottom
%po4
clearvars xp_w_po4_b pf_w_po4_b
xp = 1:366;
xp_w_po4_b = find(isnan(reg_clim_po4_b)==0);
pf_w_po4_b = polyfit(xp_w_po4_b,reg_clim_po4_b(xp_w_po4_b),3);
yp_w_po4_04_b = polyval(pf_w_po4_b,xp);
yp_w_po4_b = yp_w_po4_04_b + regm_po4_b_1;

%bottom
%po4 regime 1
clearvars xp_w_po4_b pf_w_po4_b
xp = 1:366;
xp_w_po4_b = find(isnan(reg_clim_po4_b_1)==0);
pf_w_po4_b = polyfit(xp_w_po4_b,reg_clim_po4_b(xp_w_po4_b),3);
yp_w_po4_04_b_1 = polyval(pf_w_po4_b,xp);
yp_w_po4_b_1 = yp_w_po4_04_b_1 + regm_po4_b_1;

%bottom
%po4 regime 2
clearvars xp_w_po4_b pf_w_po4_b
xp = 1:366;
xp_w_po4_b = find(isnan(reg_clim_po4_b_2)==0);
pf_w_po4_b = polyfit(xp_w_po4_b,reg_clim_po4_b(xp_w_po4_b),3);
yp_w_po4_04_b_2 = polyval(pf_w_po4_b,xp);
yp_w_po4_b_2 = yp_w_po4_04_b_2 + regm_po4_b_2;

return

%% po4
%full_adv
figure;
scatter(1:366,reg_clim_po4 + regm_po4);
hold on
plot(1:366, yp_w_po4,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('po4 (ug/L)','fontsize',13)
title('Á¤¼±°üÃø(Ç¥Ãþ)-polynomial po4.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([-inf inf])
xlim([1 366])

% whole_adv
figure;
scatter(1:366,reg_clim_po4_b + regm_po4_b_1);
hold on
plot(1:366, yp_w_po4_b,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('po4 (ug/L)','fontsize',13)
title('Á¤¼±°üÃø(ÀúÃþ)-polynomial po4.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([-inf inf])
xlim([1 366])

% 1st
figure;
scatter(1:366,reg_clim_po4_b_1 + regm_po4_b_1);
hold on
plot(1:366, yp_w_po4_b_1,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('po4 (ug/L)','fontsize',13)
title('Á¤¼±°üÃø(ÀúÃþ)-polynomial po4.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([-inf inf])
xlim([1 366])

% 2nd
figure;
scatter(1:366,reg_clim_po4_b_2 + regm_po4_b_2);
hold on
plot(1:366, yp_w_po4_b_2,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('po4 (ug/L)','fontsize',13)
title('Á¤¼±°üÃø(ÀúÃþ)-polynomial po4.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([-inf inf])
xlim([1 366])


%% no3
%full_adv
figure;
scatter(1:366,reg_clim_no3 + regm_no3_1);
hold on
plot(1:366, yp_w_no3_04,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('NO3 (ug/L)','fontsize',13)
title('Á¤¼±°üÃø(Ç¥Ãþ)-polynomial no3.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([-inf inf])
xlim([1 366])

figure;
scatter(1:366,reg_clim_no3_b + regm_no3_b_1);
hold on
plot(1:366, yp_w_no3_04_b,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('NO3 (ug/L)','fontsize',13)
title('Á¤¼±°üÃø(ÀúÃþ)-polynomial no3.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([-inf inf])
xlim([1 366])


% 1st
figure;
scatter(1:366,reg_clim_no3_1 + regm_no3_1);
hold on
plot(1:366, yp_w_no3_1,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('NO3 (ug/L)','fontsize',13)
title('Á¤¼±°üÃø(Ç¥Ãþ)-polynomial no3.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([-inf inf])
xlim([1 366])

figure;
scatter(1:366,reg_clim_no3_b_1 + regm_no3_b_1);
hold on
plot(1:366, yp_w_no3_b_1,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('NO3 (ug/L)','fontsize',13)
title('Á¤¼±°üÃø(ÀúÃþ)-polynomial no3.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([-inf inf])
xlim([1 366])

% 2nd
figure;
scatter(1:366,reg_clim_no3_2 + regm_no3_2);
hold on
plot(1:366, yp_w_no3_2,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('NO3 (ug/L)','fontsize',13)
title('Á¤¼±°üÃø(Ç¥Ãþ)-polynomial no3.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([-inf inf])
xlim([1 366])

figure;
scatter(1:366,reg_clim_no3_b_2 + regm_no3_b_2);
hold on
plot(1:366, yp_w_no3_b_2,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('NO3 (ug/L)','fontsize',13)
title('Á¤¼±°üÃø(ÀúÃþ)-polynomial no3.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([-inf inf])
xlim([1 366])

% 3rd
figure;
scatter(1:366,reg_clim_no3_3 + regm_no3_3);
hold on
plot(1:366, yp_w_no3_3,'--','color','r','linew',3)
xlabel('time(days)','fontsize',13)
ylabel('NO3 (ug/L)','fontsize',13)
title('Á¤¼±°üÃø(Ç¥Ãþ)-polynomial no3.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([-inf inf])
xlim([1 366])

figure;
scatter(1:366,reg_clim_no3_b_3 + regm_no3_b_3);
hold on
plot(1:366, yp_w_no3_b_3,'--','color','r','linew',3)
xlabel('time(days)','fontsize',13)
ylabel('NO3 (ug/L)','fontsize',13)
title('Á¤¼±°üÃø(ÀúÃþ)-polynomial no3.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([-inf inf])
xlim([1 366])

%% ALL
figure;
hold on;
plot(1:366, yp_w_no3_04,'--','color','k','linew',2)
plot(1:10:366, yp_w_no3_04_b(1:10:end),'<-','color','k','linew',2)
plot(1:366, yp_w_no3_1,'--','color','r','linew',2)
plot(1:10:366, yp_w_no3_b_1(1:10:end),'<-','color','r','linew',2)
plot(1:366, yp_w_no3_2,'--','color','g','linew',2)
plot(1:10:366, yp_w_no3_b_2(1:10:end),'<-','color','g','linew',2)
plot(1:366,yp_w_no3_3,'--','color','m','linew',2)
plot(1:10:366,yp_w_no3_b_3(1:10:end),'<-','color','m','linew',2)

plot(1:366, yp_w_no3_04_b,'--','color','k','linew',2)
plot(1:366, yp_w_no3_b_1,'--','color','r','linew',2)
plot(1:366,yp_w_no3_b_2,'--','color','g','linew',2)
plot(1:366,yp_w_no3_b_3,'--','color','m','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('NO3 (ug/L)','fontsize',13)
title('Á¤¼±°üÃø(Ç¥Ãþ&ÀúÃþ)-polynomial.','fontsize',13)
% title('Á¤¼±°üÃø(Ç¥Ãþ)-polynomial.','fontsize',13)
% title('Á¤¼±°üÃø(ÀúÃþ)-polynomial.','fontsize',13)
grid on;
% legend('surf_a_d_v','bot_a_d_v','surf_1_s_t','bot_1_s_t','surf_2_n_d','bot_2_n_d','surf_3_r_d','bot_3_r_d');
% legend('adv','1st','2nd','3rd');
set(gca,'fontsize',13)
xlim([1 366])
ylim([0 inf])

figure;
hold on;
plot(1:366, yp_w_po4,'--','color','k','linew',2)
plot(1:10:366, yp_w_po4_b(1:10:end),'<-','color','k','linew',2)
plot(1:10:366, yp_w_po4_b_1(1:10:end),'<-','color','r','linew',2)
plot(1:10:366, yp_w_po4_b_2(1:10:end),'<-','color','g','linew',2)

% plot(1:366, yp_w_po4_b,'--','color','k','linew',2)
% plot(1:366, yp_w_no3_b_1,'--','color','r','linew',2)
% plot(1:366,yp_w_no3_b_2,'--','color','g','linew',2)

xlabel('time(days)','fontsize',13)
ylabel('PO4 (ug/L)','fontsize',13)
title('Á¤¼±°üÃø(Ç¥Ãþ&ÀúÃþ)-polynomial.','fontsize',13)
% title('Á¤¼±°üÃø(Ç¥Ãþ)-polynomial.','fontsize',13)
% title('Á¤¼±°üÃø(ÀúÃþ)-polynomial.','fontsize',13)
grid on;
% legend('surf_a_d_v','bot_a_d_v','surf_1_s_t','bot_1_s_t','surf_2_n_d','bot_2_n_d','surf_3_r_d','bot_3_r_d');
% legend('adv','1st','2nd','3rd');
set(gca,'fontsize',13)
xlim([1 366])
ylim([0 inf])


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

%plot-po4
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
gray = [128/255 128/255 128/255];
    clearvars temp_s temp_b
   temp_s= po4_sur;
   temp_b= po4_bot;
figure; hold on;
 plot(temp_s,'*','color','b','linew',2);
 plot(temp_b,'*','color','r','linew',2);   
   
t=1:length(po4_sur);
%      temp_t_s(isnan(temp_t_s)) = interp1(t(~isnan(temp_t_s)), temp_t_s(~isnan(temp_t_s)), t(isnan(temp_t_s)) ); 
temp_s(isnan(temp_s)) = interp1( t(~isnan(temp_s)), temp_s(~isnan(temp_s)), t(isnan(temp_s)) ); 
temp_b(isnan(temp_b)) = interp1( t(~isnan(temp_b)), temp_b(~isnan(temp_b)), t(isnan(temp_b)) ); 

interp_nh4_b = temp_b;
interp_nh4_s = temp_s;

 plot(temp_s,'color','b');
 plot(temp_b,'color','r');  
 ylim([0 inf])
xlabel('time(year)','fontsize',13)
ylabel('po4 (ug/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:3:end)]);
set(gca,'xlim',[tx_tick(8) tx_tick(end)]);
set(gca,'xticklabel',1980:3:2019);
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',tx_tick(1:1:end));
legend('surface','bottom')


%plot-NO3
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
gray = [128/255 128/255 128/255];
    clearvars temp_s temp_b
   temp_s= no3_sur;
   temp_b= no3_bot;
figure; hold on;
 plot(temp_s,'*','color','b','linew',2);
 plot(temp_b,'*','color','r','linew',2);   
   
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
ylabel('no3 (ug/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:3:end)]);
set(gca,'xlim',[tx_tick(8) tx_tick(end)]);
set(gca,'xticklabel',1980:3:2019);
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',tx_tick(1:1:end));
legend('surface','bottom')


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

%plot-salt
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
gray = [128/255 128/255 128/255];
    clearvars salt_s salt_b
   salt_s= salt_sur;
   salt_b= salt_bot;
salt_s(find(salt_s > 90)) = NaN;
salt_b(find(salt_b > 90)) = NaN;

figure; hold on;
 plot(salt_s,'*','color','b','linew',2);
 plot(salt_b,'*','color','r','linew',2);
%  plot(salt_sur,'*','color','c','linew',2);
%  plot(salt_bot,'*','color','m','linew',2);   
   
t=1:length(no3_sur);
%      salt_t_s(isnan(salt_t_s)) = interp1(t(~isnan(salt_t_s)), salt_t_s(~isnan(salt_t_s)), t(isnan(salt_t_s)) ); 
salt_s(isnan(salt_s)) = interp1( t(~isnan(salt_s)), salt_s(~isnan(salt_s)), t(isnan(salt_s)) ); 
salt_b(isnan(salt_b)) = interp1( t(~isnan(salt_b)), salt_b(~isnan(salt_b)), t(isnan(salt_b)) ); 

interp_nh4_b = salt_b;
interp_nh4_s = salt_s;

 plot(salt_s,'color','b');
 plot(salt_b,'color','r');  
 ylim([30 inf])
xlabel('time(year)','fontsize',13)
ylabel('salt (^oC)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:3:end)]);
set(gca,'xlim',[tx_tick(1) tx_tick(end)]);
set(gca,'xticklabel',1980:3:2019);
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
text(2,3,['y = ' num2str(y_b_do_1(1),'%.3f') 'x + ' num2str(y_b_do_1(2),'%.3f')],'fontsize',13);


%plot-DO & recon.DO
temp_s_nonan=find(isnan(temp_sur)==0);
yCalc = temp_sur(temp_s_nonan) .* y_s_do_1(1) + y_s_do_1(2);

temp_b_nonan=find(isnan(temp_bot)==0);
yCalc_b = temp_bot(temp_b_nonan) .* y_b_do_1(1) + y_b_do_1(2);

clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
gray = [128/255 128/255 128/255];
    clearvars temp_s temp_b
   temp_s= do_sur;
   temp_b= do_bot;
figure; hold on;
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

set(gca,'xtick',[tx_tick(1:3:end)]);
set(gca,'xlim',[tx_tick(1) tx_tick(end)]);
set(gca,'xticklabel',1980:3:2019);
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
 plot(temp_s_nonan, do_re_2,'*','color','b','linew',2);
    
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

set(gca,'xtick',[tx_tick(1:3:end)]);
set(gca,'xlim',[tx_tick(1) tx_tick(end)]);
set(gca,'xticklabel',1980:3:2019);
set(gcf, 'Units', 'Inches', 'Position', [1, 1, 9, 7], ...
    'PaperUnits', 'Inches', 'PaperSize', [9, 7])

y_s_no3 = yp_w_no3_1;
y_b_no3 = yp_w_no3_b_1;
% save('no3_koem_input_400_16.mat','y_s_no3','y_b_no3','yp_w_no3_*','yp_w_no3_b_*');
save('po4_koem_input_400_16.mat','yp_w_*');

figure;plot(y_b_no3,'b'); hold on; plot(y_s_no3,'r');
