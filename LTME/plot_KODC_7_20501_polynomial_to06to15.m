close all; clear; clc;   % -v3

% % station name tag
% % (for preventing the miss picking using tag from 1 col's name on excel) 

% port
name_tag_1{1} = [num2str(01,'%02d')] 
 
% combining the tag and outter point excluding
% name_tag = name_tag_1'; 
% name_tag{end+1:end+length(name_tag_2)} = name_tag_2; 
% name_tag{end+1:end+length(name_tag_3)} = name_tag_3; 
% size_tag = length(name_tag);


%% pick the row on the excel which has same name with tag
[raw_p1 txt_p1]=xlsread('정선해양조사_205line_from1980to2019.xls','정선해양관측정보','');
[raw_p2 txt_p2]=xlsread('정선해양조사_205line_2020_only.xls','정선해양관측정보','');
% raw_p=cat(1,raw_p1,raw_p2);
txt_p = cat(1,txt_p1,txt_p2(3:end,:));
txt_matc_p = txt_p(3:end,3); % name list
txt_date_p = txt_p(3:end,7); % date list
txt_cut=txt_p(3:end,:);

depth_nan_p = strcmp(txt_cut(:,8),''); % detect NaN
for i = 1:length(txt_cut)
   if depth_nan_p(i) == 1 % when it's nan
       depth_p(i) = NaN;
   else
       depth_p(i) = str2num(char(txt_cut(i,8)));
   end
end

temp_nan_p = strcmp(txt_cut(:,9),''); % detect NaN
for i = 1:length(txt_cut)
   if temp_nan_p(i) == 1 % when it's nan
       temp_p(i) = NaN;
   else
       temp_p(i) = str2num(char(txt_cut(i,9)));
   end
end

salt_nan_p = strcmp(txt_cut(:,11),''); % detect NaN
for i = 1:length(txt_cut)
   if salt_nan_p(i) == 1 % when it's nan
       salt_p(i) = NaN;
   else
       salt_p(i) = str2num(char(txt_cut(i,11)));
   end
end

do_nan_p = strcmp(txt_cut(:,13),''); % detect NaN
for i = 1:length(txt_cut)
   if do_nan_p(i) == 1 % when it's nan
       do_p(i) = NaN;
   else
       do_p(i) = str2num(char(txt_cut(i,13)));
   end
end

no3_nan_p = strcmp(txt_cut(:,18),''); % detect NaN
for i = 1:length(txt_cut)
   if no3_nan_p(i) == 1 % when it's nan
       no3_p(i) = NaN;
   else
       no3_p(i) = str2num(char(txt_cut(i,18)));
   end
end

po4_nan_p = strcmp(txt_cut(:,16),''); % detect NaN
for i = 1:length(txt_cut)
   if po4_nan_p(i) == 1 % when it's nan
       po4_p(i) = NaN;
   else
       po4_p(i) = str2num(char(txt_cut(i,16)));
   end
end

no2_nan_p = strcmp(txt_cut(:,17),''); % detect NaN
for i = 1:length(txt_cut)
   if no2_nan_p(i) == 1 % when it's nan
       no2_p(i) = NaN;
   else
       no2_p(i) = str2num(char(txt_cut(i,17)));
   end
end

si_nan_p = strcmp(txt_cut(:,19),''); % detect NaN
for i = 1:length(txt_cut)
   if si_nan_p(i) == 1 % when it's nan
       si_p(i) = NaN;
   else
       si_p(i) = str2num(char(txt_cut(i,19)));
   end
end

ph_nan_p = strcmp(txt_cut(:,20),''); % detect NaN
for i = 1:length(txt_cut)
   if ph_nan_p(i) == 1 % when it's nan
       ph_p(i) = NaN;
   else
       ph_p(i) = str2num(char(txt_cut(i,20)));
   end
end

secchi_nan_p = strcmp(txt_cut(:,21),''); % detect NaN
for i = 1:length(txt_cut)
   if secchi_nan_p(i) == 1 % when it's nan
       secchi_p(i) = NaN;
   else
       secchi_p(i) = str2num(char(txt_cut(i,21)));
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
    po4_kodc{i} = po4_p(indx_st{i});
    no3_kodc{i} = no3_p(indx_st{i});
    no2_kodc{i} = no2_p(indx_st{i});
    si_kodc{i} = si_p(indx_st{i});
    ph_kodc{i} = ph_p(indx_st{i});
    secchi_kodc{i} = secchi_p(indx_st{i});
end


%% make date to be 'yymm' form
for i = 1:length(name_tag_1)
clearvars temp temp_ymd temp_ymd_c temp_yymmdd temp_yymmdd_c
temp = char(date_st{i});
temp_ymd=temp(:,1:end-9);
temp_yymmdd=temp(:,1:end-6);

for j = 1:length(temp_ymd)
    temp_ymd_c{j} = temp_ymd(j,:);
    temp_yymmdd_c{j} = temp_yymmdd(j,:);
                    
end
date_ymd{i,1} = temp_ymd_c;
date_yymmdd{i,1} = temp_yymmdd_c;
end

% make 1980~present
k=0
for i = 1980:2020
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

k=0; m=0;
for i = 1:41
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
indx_date_s = cell(1,length(ref_yymm));
indx_date_b = cell(1,length(ref_yymm));
for j = 1:length(name_tag_1) % st. axis
    for i = 1:length(ref_yymm) % date axis
       if  sum(strcmp(ref_yymm{i}, date_ymd{j})) ~= 0
           clearvars indx_date
           indx_date = find([strcmp(ref_yymm{i}, date_ymd{j})] == 1); 
           indx_date_s{j,i} = indx_date(find(depth_kodc{j}(1,indx_date) == 0));
           indx_date_b{j,i} = indx_date(find(depth_kodc{j}(1,indx_date) == max(depth_kodc{j}(1,indx_date))));
       end
    end 
end

% % yymmdd
% indx_date_s = cell(1,length(ref_yymmdd));
% indx_date_b = cell(1,length(ref_yymmdd));
% for j = 1:length(name_tag_1) % st. axis
%     for i = 1:length(ref_yymmdd) % date axis
%        if  sum(strcmp(ref_yymmdd{i}, date_yymmdd{j})) ~= 0
%            clearvars indx_date
%            indx_date = find([strcmp(ref_yymmdd{i}, date_yymmdd{j})] == 1); 
%            indx_date_s{j,i} = indx_date(find(depth_kodc{j}(1,indx_date) == 0));
%            indx_date_b{j,i} = indx_date(find(depth_kodc{j}(1,indx_date) == max(depth_kodc{j}(1,indx_date))));
% %        elseif sum(strcmp(ref_yymmdd{i}, date_yymmdd{j})) == 0
% %            indx_date_s{j,i} = NaN;
% %            indx_date_b{j,i} = NaN;
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

%no2
for i = 1:size(indx_date_s,1)
    clearvars temp_kodc_k
    temp_kodc_k = no2_kodc{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_kodc_k(indx_date_s{i,j}))) ~= 0     
        no2_sur(i,j) = nanmean(temp_kodc_k(indx_date_s{i,j}));
    elseif sum(size(temp_kodc_k(indx_date_s{i,j}))) == 0 
        no2_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_kodc_k
    temp_kodc_k = no2_kodc{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_kodc_k(indx_date_b{i,j}))) ~= 0     
        no2_bot(i,j) = nanmean(temp_kodc_k(indx_date_b{i,j}));
    elseif sum(size(temp_kodc_k(indx_date_b{i,j}))) == 0 
        no2_bot(i,j) = NaN;
    end
    end
end

%si
for i = 1:size(indx_date_s,1)
    clearvars temp_kodc_k
    temp_kodc_k = si_kodc{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_kodc_k(indx_date_s{i,j}))) ~= 0     
        si_sur(i,j) = nanmean(temp_kodc_k(indx_date_s{i,j}));
    elseif sum(size(temp_kodc_k(indx_date_s{i,j}))) == 0 
        si_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_kodc_k
    temp_kodc_k = si_kodc{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_kodc_k(indx_date_b{i,j}))) ~= 0     
        si_bot(i,j) = nanmean(temp_kodc_k(indx_date_b{i,j}));
    elseif sum(size(temp_kodc_k(indx_date_b{i,j}))) == 0 
        si_bot(i,j) = NaN;
    end
    end
end

%secchi
for i = 1:size(indx_date_s,1)
    clearvars temp_kodc_k
    temp_kodc_k = secchi_kodc{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_kodc_k(indx_date_s{i,j}))) ~= 0     
        secchi_sur(i,j) = nanmean(temp_kodc_k(indx_date_s{i,j}));
    elseif sum(size(temp_kodc_k(indx_date_s{i,j}))) == 0 
        secchi_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_kodc_k
    temp_kodc_k = secchi_kodc{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_kodc_k(indx_date_b{i,j}))) ~= 0     
        secchi_bot(i,j) = nanmean(temp_kodc_k(indx_date_b{i,j}));
    elseif sum(size(temp_kodc_k(indx_date_b{i,j}))) == 0 
        secchi_bot(i,j) = NaN;
    end
    end
end

%ph
for i = 1:size(indx_date_s,1)
    clearvars temp_kodc_k
    temp_kodc_k = ph_kodc{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_kodc_k(indx_date_s{i,j}))) ~= 0     
        ph_sur(i,j) = nanmean(temp_kodc_k(indx_date_s{i,j}));
    elseif sum(size(temp_kodc_k(indx_date_s{i,j}))) == 0 
        ph_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_kodc_k
    temp_kodc_k = ph_kodc{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_kodc_k(indx_date_b{i,j}))) ~= 0     
        ph_bot(i,j) = nanmean(temp_kodc_k(indx_date_b{i,j}));
    elseif sum(size(temp_kodc_k(indx_date_b{i,j}))) == 0 
        ph_bot(i,j) = NaN;
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
no3_sur =  no3_sur(1,:)';  %umol/L -> ug/L
no3_bot = no3_bot(1,:)';
temp_sur = temp_sur(1,:)';
temp_bot = temp_bot(1,:)';
salt_sur = salt_sur(1,:)';
salt_bot = salt_bot(1,:)';
po4_sur = po4_sur(1,:)';
po4_bot = po4_bot(1,:)';
no2_sur = no2_sur(1,:)';
no2_bot = no2_bot(1,:)';
ph_sur = ph_sur(1,:)';
ph_bot = ph_bot(1,:)';
secchi_sur = secchi_sur(1,:)';
secchi_bot = secchi_bot(1,:)';
si_sur = si_sur(1,:)';
si_bot = si_bot(1,:)';

%% cut 3 regime
cut_0 = find(strcmp(ref_yymm, {'1996-12'}) ~= 0) % 2006-12
cut_1 = find(strcmp(ref_yymm, {'2006-12'}) ~= 0) % 2006-12
cut_2 = find(strcmp(ref_yymm, {'2015-12'}) ~= 0) % 2015-12

do_sur_1 = do_sur(cut_0+1:cut_1);
do_sur_2 = do_sur(cut_1+1:cut_2);
do_sur_3 = do_sur(cut_2+1:end);

do_bot_1 = do_bot(cut_0+1:cut_1);
do_bot_2 = do_bot(cut_1+1:cut_2);
do_bot_3 = do_bot(cut_2+1:end);

no3_sur_1 = no3_sur(cut_0+1:cut_1);
no3_sur_2 = no3_sur(cut_1+1:cut_2);
no3_sur_3 = no3_sur(cut_2+1:end);

no3_bot_1 = no3_bot(cut_0+1:cut_1);
no3_bot_2 = no3_bot(cut_1+1:cut_2);
no3_bot_3 = no3_bot(cut_2+1:end);

temp_sur_1 = temp_sur(cut_0+1:cut_1);
temp_sur_2 = temp_sur(cut_1+1:cut_2);
temp_sur_3 = temp_sur(cut_2+1:end);

temp_bot_1 = temp_bot(cut_0+1:cut_1);
temp_bot_2 = temp_bot(cut_1+1:cut_2);
temp_bot_3 = temp_bot(cut_2+1:end);

salt_sur_1 = salt_sur(cut_0+1:cut_1);
salt_sur_2 = salt_sur(cut_1+1:cut_2);
salt_sur_3 = salt_sur(cut_2+1:end);

salt_bot_1 = salt_bot(cut_0+1:cut_1);
salt_bot_2 = salt_bot(cut_1+1:cut_2);
salt_bot_3 = salt_bot(cut_2+1:end);

po4_sur_1 = po4_sur(cut_0+1:cut_1);
po4_sur_2 = po4_sur(cut_1+1:cut_2);
po4_sur_3 = po4_sur(cut_2+1:end);

po4_bot_1 = po4_bot(cut_0+1:cut_1);
po4_bot_2 = po4_bot(cut_1+1:cut_2);
po4_bot_3 = po4_bot(cut_2+1:end);

no2_sur_1 = no2_sur(cut_0+1:cut_1);
no2_sur_2 = no2_sur(cut_1+1:cut_2);
no2_sur_3 = no2_sur(cut_2+1:end);

no2_bot_1 = no2_bot(cut_0+1:cut_1);
no2_bot_2 = no2_bot(cut_1+1:cut_2);
no2_bot_3 = no2_bot(cut_2+1:end);

ph_sur_1 = ph_sur(cut_0+1:cut_1);
ph_sur_2 = ph_sur(cut_1+1:cut_2);
ph_sur_3 = ph_sur(cut_2+1:end);

ph_bot_1 = ph_bot(cut_0+1:cut_1);
ph_bot_2 = ph_bot(cut_1+1:cut_2);
ph_bot_3 = ph_bot(cut_2+1:end);

secchi_sur_1 = secchi_sur(cut_0+1:cut_1);
secchi_sur_2 = secchi_sur(cut_1+1:cut_2);
secchi_sur_3 = secchi_sur(cut_2+1:end);

secchi_bot_1 = secchi_bot(cut_0+1:cut_1);
secchi_bot_2 = secchi_bot(cut_1+1:cut_2);
secchi_bot_3 = secchi_bot(cut_2+1:end);

si_sur_1 = si_sur(cut_0+1:cut_1);
si_sur_2 = si_sur(cut_1+1:cut_2);
si_sur_3 = si_sur(cut_2+1:end);

si_bot_1 = si_bot(cut_0+1:cut_1);
si_bot_2 = si_bot(cut_1+1:cut_2);
si_bot_3 = si_bot(cut_2+1:end);

%% extract over 3sig
sig = 3; %% sigma
%% extract sigma
sig = 3; 

clearvars tempd
tempd=no3_sur_1;
tempd(tempd > nanmean(no3_sur_1) + sig*nanstd(no3_sur_1)) =NaN;
tempd(tempd < nanmean(no3_sur_1) - sig*nanstd(no3_sur_1)) =NaN;
no3_sur_1_s = tempd;

clearvars tempd
tempd=no3_sur_2;
tempd(tempd > nanmean(no3_sur_2) + sig*nanstd(no3_sur_2)) =NaN;
tempd(tempd < nanmean(no3_sur_2) - sig*nanstd(no3_sur_2)) =NaN;
no3_sur_2_s = tempd;

clearvars tempd
tempd=no3_sur_3;
tempd(tempd > nanmean(no3_sur_3) + sig*nanstd(no3_sur_3)) =NaN;
tempd(tempd < nanmean(no3_sur_3) - sig*nanstd(no3_sur_3)) =NaN;
no3_sur_3_s = tempd;

clearvars tempd
tempd=no3_bot_1;
tempd(tempd > nanmean(no3_bot_1) + sig*nanstd(no3_bot_1)) =NaN;
tempd(tempd < nanmean(no3_bot_1) - sig*nanstd(no3_bot_1)) =NaN;
no3_bot_1_s = tempd;

clearvars tempd
tempd=no3_bot_2;
tempd(tempd > nanmean(no3_bot_2) + sig*nanstd(no3_bot_2)) =NaN;
tempd(tempd < nanmean(no3_bot_2) - sig*nanstd(no3_bot_2)) =NaN;
no3_bot_2_s = tempd;

clearvars tempd
tempd=no3_bot_3;
tempd(tempd > nanmean(no3_bot_3) + sig*nanstd(no3_bot_3)) =NaN;
tempd(tempd < nanmean(no3_bot_3) - sig*nanstd(no3_bot_3)) =NaN;
no3_bot_3_s = tempd;

clearvars tempd
tempd=po4_sur_1;
tempd(tempd > nanmean(po4_sur_1) + sig*nanstd(po4_sur_1)) =NaN;
tempd(tempd < nanmean(po4_sur_1) - sig*nanstd(po4_sur_1)) =NaN;
po4_sur_1_s = tempd;

clearvars tempd
tempd=po4_sur_2;
tempd(tempd > nanmean(po4_sur_2) + sig*nanstd(po4_sur_2)) =NaN;
tempd(tempd < nanmean(po4_sur_2) - sig*nanstd(po4_sur_2)) =NaN;
po4_sur_2_s = tempd;

clearvars tempd
tempd=po4_sur_3;
tempd(tempd > nanmean(po4_sur_3) + sig*nanstd(po4_sur_3)) =NaN;
tempd(tempd < nanmean(po4_sur_3) - sig*nanstd(po4_sur_3)) =NaN;
po4_sur_3_s = tempd;

clearvars tempd
tempd=po4_bot_1;
tempd(tempd > nanmean(po4_bot_1) + sig*nanstd(po4_bot_1)) =NaN;
tempd(tempd < nanmean(po4_bot_1) - sig*nanstd(po4_bot_1)) =NaN;
po4_bot_1_s = tempd;

clearvars tempd
tempd=po4_bot_2;
tempd(tempd > nanmean(po4_bot_2) + sig*nanstd(po4_bot_2)) =NaN;
tempd(tempd < nanmean(po4_bot_2) - sig*nanstd(po4_bot_2)) =NaN;
po4_bot_2_s = tempd;

clearvars tempd
tempd=po4_bot_3;
tempd(tempd > nanmean(po4_bot_3) + sig*nanstd(po4_bot_3)) =NaN;
tempd(tempd < nanmean(po4_bot_3) - sig*nanstd(po4_bot_3)) =NaN;
po4_bot_3_s = tempd;

clearvars tempd
tempd=do_sur_1;
tempd(tempd > nanmean(do_sur_1) + sig*nanstd(do_sur_1)) =NaN;
tempd(tempd < nanmean(do_sur_1) - sig*nanstd(do_sur_1)) =NaN;
do_sur_1_s = tempd;

clearvars tempd
tempd=do_sur_2;
tempd(tempd > nanmean(do_sur_2) + sig*nanstd(do_sur_2)) =NaN;
tempd(tempd < nanmean(do_sur_2) - sig*nanstd(do_sur_2)) =NaN;
do_sur_2_s = tempd;

clearvars tempd
tempd=do_sur_3;
tempd(tempd > nanmean(do_sur_3) + sig*nanstd(do_sur_3)) =NaN;
tempd(tempd < nanmean(do_sur_3) - sig*nanstd(do_sur_3)) =NaN;
do_sur_3_s = tempd;

clearvars tempd
tempd=do_bot_1;
tempd(tempd > nanmean(do_bot_1) + sig*nanstd(do_bot_1)) =NaN;
tempd(tempd < nanmean(do_bot_1) - sig*nanstd(do_bot_1)) =NaN;
do_bot_1_s = tempd;

clearvars tempd
tempd=do_bot_2;
tempd(tempd > nanmean(do_bot_2) + sig*nanstd(do_bot_2)) =NaN;
tempd(tempd < nanmean(do_bot_2) - sig*nanstd(do_bot_2)) =NaN;
do_bot_2_s = tempd;

clearvars tempd
tempd=do_bot_3;
tempd(tempd > nanmean(do_bot_3) + sig*nanstd(do_bot_3)) =NaN;
tempd(tempd < nanmean(do_bot_3) - sig*nanstd(do_bot_3)) =NaN;
do_bot_3_s = tempd;

clearvars tempd
tempd=temp_sur_1;
tempd(tempd > nanmean(temp_sur_1) + sig*nanstd(temp_sur_1)) =NaN;
tempd(tempd < nanmean(temp_sur_1) - sig*nanstd(temp_sur_1)) =NaN;
temp_sur_1_s = tempd;

clearvars tempd
tempd=temp_sur_2;
tempd(tempd > nanmean(temp_sur_2) + sig*nanstd(temp_sur_2)) =NaN;
tempd(tempd < nanmean(temp_sur_2) - sig*nanstd(temp_sur_2)) =NaN;
temp_sur_2_s = tempd;

clearvars tempd
tempd=temp_sur_3;
tempd(tempd > nanmean(temp_sur_3) + sig*nanstd(temp_sur_3)) =NaN;
tempd(tempd < nanmean(temp_sur_3) - sig*nanstd(temp_sur_3)) =NaN;
temp_sur_3_s = tempd;

clearvars tempd
tempd=temp_bot_1;
tempd(tempd > nanmean(temp_bot_1) + sig*nanstd(temp_bot_1)) =NaN;
tempd(tempd < nanmean(temp_bot_1) - sig*nanstd(temp_bot_1)) =NaN;
temp_bot_1_s = tempd;

clearvars tempd
tempd=temp_bot_2;
tempd(tempd > nanmean(temp_bot_2) + sig*nanstd(temp_bot_2)) =NaN;
tempd(tempd < nanmean(temp_bot_2) - sig*nanstd(temp_bot_2)) =NaN;
temp_bot_2_s = tempd;

clearvars tempd
tempd=temp_bot_3;
tempd(tempd > nanmean(temp_bot_3) + sig*nanstd(temp_bot_3)) =NaN;
tempd(tempd < nanmean(temp_bot_3) - sig*nanstd(temp_bot_3)) =NaN;
temp_bot_3_s = tempd;

clearvars tempd
tempd=salt_sur_1;
tempd(tempd > nanmean(salt_sur_1) + sig*nanstd(salt_sur_1)) =NaN;
tempd(tempd < nanmean(salt_sur_1) - sig*nanstd(salt_sur_1)) =NaN;
salt_sur_1_s = tempd;

clearvars tempd
tempd=salt_sur_2;
tempd(tempd > nanmean(salt_sur_2) + sig*nanstd(salt_sur_2)) =NaN;
tempd(tempd < nanmean(salt_sur_2) - sig*nanstd(salt_sur_2)) =NaN;
salt_sur_2_s = tempd;

clearvars tempd
tempd=salt_sur_3;
tempd(tempd > nanmean(salt_sur_3) + sig*nanstd(salt_sur_3)) =NaN;
tempd(tempd < nanmean(salt_sur_3) - sig*nanstd(salt_sur_3)) =NaN;
salt_sur_3_s = tempd;

clearvars tempd
tempd=salt_bot_1;
tempd(tempd > nanmean(salt_bot_1) + sig*nanstd(salt_bot_1)) =NaN;
tempd(tempd < nanmean(salt_bot_1) - sig*nanstd(salt_bot_1)) =NaN;
salt_bot_1_s = tempd;

clearvars tempd
tempd=salt_bot_2;
tempd(tempd > nanmean(salt_bot_2) + sig*nanstd(salt_bot_2)) =NaN;
tempd(tempd < nanmean(salt_bot_2) - sig*nanstd(salt_bot_2)) =NaN;
salt_bot_2_s = tempd;

clearvars tempd
tempd=salt_bot_3;
tempd(tempd > nanmean(salt_bot_3) + sig*nanstd(salt_bot_3)) =NaN;
tempd(tempd < nanmean(salt_bot_3) - sig*nanstd(salt_bot_3)) =NaN;
salt_bot_3_s = tempd;

clearvars tempd
tempd=secchi_sur_1;
tempd(tempd > nanmean(secchi_sur_1) + sig*nanstd(secchi_sur_1)) =NaN;
tempd(tempd < nanmean(secchi_sur_1) - sig*nanstd(secchi_sur_1)) =NaN;
secchi_sur_1_s = tempd;

clearvars tempd
tempd=secchi_sur_2;
tempd(tempd > nanmean(secchi_sur_2) + sig*nanstd(secchi_sur_2)) =NaN;
tempd(tempd < nanmean(secchi_sur_2) - sig*nanstd(secchi_sur_2)) =NaN;
secchi_sur_2_s = tempd;

clearvars tempd
tempd=secchi_sur_3;
tempd(tempd > nanmean(secchi_sur_3) + sig*nanstd(secchi_sur_3)) =NaN;
tempd(tempd < nanmean(secchi_sur_3) - sig*nanstd(secchi_sur_3)) =NaN;
secchi_sur_3_s = tempd;

clearvars tempd
tempd=secchi_bot_1;
tempd(tempd > nanmean(secchi_bot_1) + sig*nanstd(secchi_bot_1)) =NaN;
tempd(tempd < nanmean(secchi_bot_1) - sig*nanstd(secchi_bot_1)) =NaN;
secchi_bot_1_s = tempd;

clearvars tempd
tempd=secchi_bot_2;
tempd(tempd > nanmean(secchi_bot_2) + sig*nanstd(secchi_bot_2)) =NaN;
tempd(tempd < nanmean(secchi_bot_2) - sig*nanstd(secchi_bot_2)) =NaN;
secchi_bot_2_s = tempd;

clearvars tempd
tempd=secchi_bot_3;
tempd(tempd > nanmean(secchi_bot_3) + sig*nanstd(secchi_bot_3)) =NaN;
tempd(tempd < nanmean(secchi_bot_3) - sig*nanstd(secchi_bot_3)) =NaN;
secchi_bot_3_s = tempd;

clearvars tempd
tempd=ph_sur_1;
tempd(tempd > nanmean(ph_sur_1) + sig*nanstd(ph_sur_1)) =NaN;
tempd(tempd < nanmean(ph_sur_1) - sig*nanstd(ph_sur_1)) =NaN;
ph_sur_1_s = tempd;

clearvars tempd
tempd=ph_sur_2;
tempd(tempd > nanmean(ph_sur_2) + sig*nanstd(ph_sur_2)) =NaN;
tempd(tempd < nanmean(ph_sur_2) - sig*nanstd(ph_sur_2)) =NaN;
ph_sur_2_s = tempd;

clearvars tempd
tempd=ph_sur_3;
tempd(tempd > nanmean(ph_sur_3) + sig*nanstd(ph_sur_3)) =NaN;
tempd(tempd < nanmean(ph_sur_3) - sig*nanstd(ph_sur_3)) =NaN;
ph_sur_3_s = tempd;

clearvars tempd
tempd=ph_bot_1;
tempd(tempd > nanmean(ph_bot_1) + sig*nanstd(ph_bot_1)) =NaN;
tempd(tempd < nanmean(ph_bot_1) - sig*nanstd(ph_bot_1)) =NaN;
ph_bot_1_s = tempd;

clearvars tempd
tempd=ph_bot_2;
tempd(tempd > nanmean(ph_bot_2) + sig*nanstd(ph_bot_2)) =NaN;
tempd(tempd < nanmean(ph_bot_2) - sig*nanstd(ph_bot_2)) =NaN;
ph_bot_2_s = tempd;

clearvars tempd
tempd=ph_bot_3;
tempd(tempd > nanmean(ph_bot_3) + sig*nanstd(ph_bot_3)) =NaN;
tempd(tempd < nanmean(ph_bot_3) - sig*nanstd(ph_bot_3)) =NaN;
ph_bot_3_s = tempd;

clearvars tempd
tempd=no2_sur_1;
tempd(tempd > nanmean(no2_sur_1) + sig*nanstd(no2_sur_1)) =NaN;
tempd(tempd < nanmean(no2_sur_1) - sig*nanstd(no2_sur_1)) =NaN;
no2_sur_1_s = tempd;

clearvars tempd
tempd=no2_sur_2;
tempd(tempd > nanmean(no2_sur_2) + sig*nanstd(no2_sur_2)) =NaN;
tempd(tempd < nanmean(no2_sur_2) - sig*nanstd(no2_sur_2)) =NaN;
no2_sur_2_s = tempd;

clearvars tempd
tempd=no2_sur_3;
tempd(tempd > nanmean(no2_sur_3) + sig*nanstd(no2_sur_3)) =NaN;
tempd(tempd < nanmean(no2_sur_3) - sig*nanstd(no2_sur_3)) =NaN;
no2_sur_3_s = tempd;

clearvars tempd
tempd=no2_bot_1;
tempd(tempd > nanmean(no2_bot_1) + sig*nanstd(no2_bot_1)) =NaN;
tempd(tempd < nanmean(no2_bot_1) - sig*nanstd(no2_bot_1)) =NaN;
no2_bot_1_s = tempd;

clearvars tempd
tempd=no2_bot_2;
tempd(tempd > nanmean(no2_bot_2) + sig*nanstd(no2_bot_2)) =NaN;
tempd(tempd < nanmean(no2_bot_2) - sig*nanstd(no2_bot_2)) =NaN;
no2_bot_2_s = tempd;

clearvars tempd
tempd=no2_bot_3;
tempd(tempd > nanmean(no2_bot_3) + sig*nanstd(no2_bot_3)) =NaN;
tempd(tempd < nanmean(no2_bot_3) - sig*nanstd(no2_bot_3)) =NaN;
no2_bot_3_s = tempd;

clearvars tempd
tempd=si_sur_1;
tempd(tempd > nanmean(si_sur_1) + sig*nanstd(si_sur_1)) =NaN;
tempd(tempd < nanmean(si_sur_1) - sig*nanstd(si_sur_1)) =NaN;
si_sur_1_s = tempd;

clearvars tempd
tempd=si_sur_2;
tempd(tempd > nanmean(si_sur_2) + sig*nanstd(si_sur_2)) =NaN;
tempd(tempd < nanmean(si_sur_2) - sig*nanstd(si_sur_2)) =NaN;
si_sur_2_s = tempd;

clearvars tempd
tempd=si_sur_3;
tempd(tempd > nanmean(si_sur_3) + sig*nanstd(si_sur_3)) =NaN;
tempd(tempd < nanmean(si_sur_3) - sig*nanstd(si_sur_3)) =NaN;
si_sur_3_s = tempd;

clearvars tempd
tempd=si_bot_1;
tempd(tempd > nanmean(si_bot_1) + sig*nanstd(si_bot_1)) =NaN;
tempd(tempd < nanmean(si_bot_1) - sig*nanstd(si_bot_1)) =NaN;
si_bot_1_s = tempd;

clearvars tempd
tempd=si_bot_2;
tempd(tempd > nanmean(si_bot_2) + sig*nanstd(si_bot_2)) =NaN;
tempd(tempd < nanmean(si_bot_2) - sig*nanstd(si_bot_2)) =NaN;
si_bot_2_s = tempd;

clearvars tempd
tempd=si_bot_3;
tempd(tempd > nanmean(si_bot_3) + sig*nanstd(si_bot_3)) =NaN;
tempd(tempd < nanmean(si_bot_3) - sig*nanstd(si_bot_3)) =NaN;
si_bot_3_s = tempd;

% define
mon_clim_sur_ph_1 = NaN(1,12);
mon_clim_sur_ph_2 = NaN(1,12);
mon_clim_sur_ph_3 = NaN(1,12);

mon_clim_bot_ph_1 = NaN(1,12);
mon_clim_bot_ph_2 = NaN(1,12);
mon_clim_bot_ph_3 = NaN(1,12);

mon_clim_sur_no2_1 = NaN(1,12);
mon_clim_sur_no2_2 = NaN(1,12);
mon_clim_sur_no2_3 = NaN(1,12);

mon_clim_bot_no2_1 = NaN(1,12);
mon_clim_bot_no2_2 = NaN(1,12);
mon_clim_bot_no2_3 = NaN(1,12);

mon_clim_sur_secchi_1 = NaN(1,12);
mon_clim_sur_secchi_2 = NaN(1,12);
mon_clim_sur_secchi_3 = NaN(1,12); 

mon_clim_bot_secchi_1 = NaN(1,12);
mon_clim_bot_secchi_2 = NaN(1,12);
mon_clim_bot_secchi_3 = NaN(1,12); 

mon_clim_sur_si_1 = NaN(1,12); 
mon_clim_sur_si_2 = NaN(1,12); 
mon_clim_sur_si_3 = NaN(1,12); 

mon_clim_bot_si_1 = NaN(1,12); 
mon_clim_bot_si_2 = NaN(1,12); 
mon_clim_bot_si_3 = NaN(1,12);  

mon_clim_sur_no3_1 = NaN(1,12); 
mon_clim_sur_no3_2 = NaN(1,12); 
mon_clim_sur_no3_3 = NaN(1,12);  

mon_clim_bot_no3_1 = NaN(1,12); 
mon_clim_bot_no3_2 = NaN(1,12); 
mon_clim_bot_no3_3 = NaN(1,12);  

mon_clim_sur_po4_1 = NaN(1,12); 
mon_clim_sur_po4_2 = NaN(1,12); 
mon_clim_sur_po4_3 = NaN(1,12);  

mon_clim_bot_po4_1 = NaN(1,12); 
mon_clim_bot_po4_2 = NaN(1,12);  
mon_clim_bot_po4_3 = NaN(1,12); 

mon_clim_sur_do_1 = NaN(1,12); 
mon_clim_sur_do_2 = NaN(1,12); 
mon_clim_sur_do_3 = NaN(1,12);  

mon_clim_bot_do_1 = NaN(1,12); 
mon_clim_bot_do_2 = NaN(1,12); 
mon_clim_bot_do_3 = NaN(1,12);  

mon_clim_sur_temp_1 = NaN(1,12); 
mon_clim_sur_temp_2 = NaN(1,12); 
mon_clim_sur_temp_3 = NaN(1,12);  

mon_clim_bot_temp_1 = NaN(1,12); 
mon_clim_bot_temp_2 = NaN(1,12); 
mon_clim_bot_temp_3 = NaN(1,12); 

mon_clim_sur_salt_1 = NaN(1,12); 
mon_clim_sur_salt_2 = NaN(1,12); 
mon_clim_sur_salt_3 = NaN(1,12);  

mon_clim_bot_salt_1 = NaN(1,12); 
mon_clim_bot_salt_2 = NaN(1,12); 
mon_clim_bot_salt_3 = NaN(1,12);  

% make monthly climate
for i = 2:2:12 %month   
mon_clim_sur_ph_1(i) = nanmean(ph_sur_1_s(i:12:end));
mon_clim_sur_ph_2(i) = nanmean(ph_sur_2_s(i:12:end));
mon_clim_sur_ph_3(i) = nanmean(ph_sur_3_s(i:12:end)); 

mon_clim_bot_ph_1(i) = nanmean(ph_bot_1_s(i:12:end));
mon_clim_bot_ph_2(i) = nanmean(ph_bot_2_s(i:12:end));
mon_clim_bot_ph_3(i) = nanmean(ph_bot_3_s(i:12:end)); 

mon_clim_sur_no2_1(i) = nanmean(no2_sur_1_s(i:12:end));
mon_clim_sur_no2_2(i) = nanmean(no2_sur_2_s(i:12:end));
mon_clim_sur_no2_3(i) = nanmean(no2_sur_3_s(i:12:end)); 

mon_clim_bot_no2_1(i) = nanmean(no2_bot_1_s(i:12:end));
mon_clim_bot_no2_2(i) = nanmean(no2_bot_2_s(i:12:end));
mon_clim_bot_no2_3(i) = nanmean(no2_bot_3_s(i:12:end)); 

mon_clim_sur_secchi_1(i) = nanmean(secchi_sur_1_s(i:12:end));
mon_clim_sur_secchi_2(i) = nanmean(secchi_sur_2_s(i:12:end));
mon_clim_sur_secchi_3(i) = nanmean(secchi_sur_3_s(i:12:end)); 

mon_clim_bot_secchi_1(i) = nanmean(secchi_bot_1_s(i:12:end));
mon_clim_bot_secchi_2(i) = nanmean(secchi_bot_2_s(i:12:end));
mon_clim_bot_secchi_3(i) = nanmean(secchi_bot_3_s(i:12:end)); 

mon_clim_sur_si_1(i) = nanmean(si_sur_1_s(i:12:end));
mon_clim_sur_si_2(i) = nanmean(si_sur_2_s(i:12:end));
mon_clim_sur_si_3(i) = nanmean(si_sur_3_s(i:12:end)); 

mon_clim_bot_si_1(i) = nanmean(si_bot_1_s(i:12:end));
mon_clim_bot_si_2(i) = nanmean(si_bot_2_s(i:12:end));
mon_clim_bot_si_3(i) = nanmean(si_bot_3_s(i:12:end)); 

mon_clim_sur_no3_1(i) = nanmean(no3_sur_1_s(i:12:end));
mon_clim_sur_no3_2(i) = nanmean(no3_sur_2_s(i:12:end));
mon_clim_sur_no3_3(i) = nanmean(no3_sur_3_s(i:12:end)); 

mon_clim_bot_no3_1(i) = nanmean(no3_bot_1_s(i:12:end));
mon_clim_bot_no3_2(i) = nanmean(no3_bot_2_s(i:12:end));
mon_clim_bot_no3_3(i) = nanmean(no3_bot_3_s(i:12:end)); 

mon_clim_sur_po4_1(i) = nanmean(po4_sur_1_s(i:12:end));
mon_clim_sur_po4_2(i) = nanmean(po4_sur_2_s(i:12:end));
mon_clim_sur_po4_3(i) = nanmean(po4_sur_3_s(i:12:end)); 

mon_clim_bot_po4_1(i) = nanmean(po4_bot_1_s(i:12:end));
mon_clim_bot_po4_2(i) = nanmean(po4_bot_2_s(i:12:end));
mon_clim_bot_po4_3(i) = nanmean(po4_bot_3_s(i:12:end)); 

mon_clim_sur_do_1(i) = nanmean(do_sur_1_s(i:12:end));
mon_clim_sur_do_2(i) = nanmean(do_sur_2_s(i:12:end));
mon_clim_sur_do_3(i) = nanmean(do_sur_3_s(i:12:end)); 

mon_clim_bot_do_1(i) = nanmean(do_bot_1_s(i:12:end));
mon_clim_bot_do_2(i) = nanmean(do_bot_2_s(i:12:end));
mon_clim_bot_do_3(i) = nanmean(do_bot_3_s(i:12:end)); 

mon_clim_sur_temp_1(i) = nanmean(temp_sur_1_s(i:12:end));
mon_clim_sur_temp_2(i) = nanmean(temp_sur_2_s(i:12:end));
mon_clim_sur_temp_3(i) = nanmean(temp_sur_3_s(i:12:end)); 

mon_clim_bot_temp_1(i) = nanmean(temp_bot_1_s(i:12:end));
mon_clim_bot_temp_2(i) = nanmean(temp_bot_2_s(i:12:end));
mon_clim_bot_temp_3(i) = nanmean(temp_bot_3_s(i:12:end)); 

mon_clim_sur_salt_1(i) = nanmean(salt_sur_1_s(i:12:end));
mon_clim_sur_salt_2(i) = nanmean(salt_sur_2_s(i:12:end));
mon_clim_sur_salt_3(i) = nanmean(salt_sur_3_s(i:12:end)); 

mon_clim_bot_salt_1(i) = nanmean(salt_bot_1_s(i:12:end));
mon_clim_bot_salt_2(i) = nanmean(salt_bot_2_s(i:12:end));
mon_clim_bot_salt_3(i) = nanmean(salt_bot_3_s(i:12:end)); 
end

%% make it new matrix for edge interp
mon_clim_sur_secchi_1 = [mon_clim_sur_secchi_1, mon_clim_sur_secchi_1, mon_clim_sur_secchi_1]
mon_clim_sur_secchi_2 = [mon_clim_sur_secchi_2, mon_clim_sur_secchi_2, mon_clim_sur_secchi_2]
mon_clim_sur_secchi_3 = [mon_clim_sur_secchi_3, mon_clim_sur_secchi_3, mon_clim_sur_secchi_3]

mon_clim_bot_secchi_1 = [mon_clim_bot_secchi_1, mon_clim_bot_secchi_1, mon_clim_bot_secchi_1]
mon_clim_bot_secchi_2 = [mon_clim_bot_secchi_2, mon_clim_bot_secchi_2, mon_clim_bot_secchi_2]
mon_clim_bot_secchi_3 = [mon_clim_bot_secchi_3, mon_clim_bot_secchi_3, mon_clim_bot_secchi_3]

mon_clim_sur_no3_1 = [mon_clim_sur_no3_1, mon_clim_sur_no3_1, mon_clim_sur_no3_1]
mon_clim_sur_no3_2 = [mon_clim_sur_no3_2, mon_clim_sur_no3_2, mon_clim_sur_no3_2]
mon_clim_sur_no3_3 = [mon_clim_sur_no3_3, mon_clim_sur_no3_3, mon_clim_sur_no3_3]

mon_clim_bot_no3_1 = [mon_clim_bot_no3_1, mon_clim_bot_no3_1, mon_clim_bot_no3_1]
mon_clim_bot_no3_2 = [mon_clim_bot_no3_2, mon_clim_bot_no3_2, mon_clim_bot_no3_2]
mon_clim_bot_no3_3 = [mon_clim_bot_no3_3, mon_clim_bot_no3_3, mon_clim_bot_no3_3]

mon_clim_sur_no2_1 = [mon_clim_sur_no2_1, mon_clim_sur_no2_1, mon_clim_sur_no2_1]
mon_clim_sur_no2_2 = [mon_clim_sur_no2_2, mon_clim_sur_no2_2, mon_clim_sur_no2_2]
mon_clim_sur_no2_3 = [mon_clim_sur_no2_3, mon_clim_sur_no2_3, mon_clim_sur_no2_3]

mon_clim_bot_no2_1 = [mon_clim_bot_no2_1, mon_clim_bot_no2_1, mon_clim_bot_no2_1]
mon_clim_bot_no2_2 = [mon_clim_bot_no2_2, mon_clim_bot_no2_2, mon_clim_bot_no2_2]
mon_clim_bot_no2_3 = [mon_clim_bot_no2_3, mon_clim_bot_no2_3, mon_clim_bot_no2_3]

mon_clim_sur_po4_1 = [mon_clim_sur_po4_1, mon_clim_sur_po4_1, mon_clim_sur_po4_1]
mon_clim_sur_po4_2 = [mon_clim_sur_po4_2, mon_clim_sur_po4_2, mon_clim_sur_po4_2]
mon_clim_sur_po4_3 = [mon_clim_sur_po4_3, mon_clim_sur_po4_3, mon_clim_sur_po4_3]

mon_clim_bot_po4_1 = [mon_clim_bot_po4_1, mon_clim_bot_po4_1, mon_clim_bot_po4_1]
mon_clim_bot_po4_2 = [mon_clim_bot_po4_2, mon_clim_bot_po4_2, mon_clim_bot_po4_2]
mon_clim_bot_po4_3 = [mon_clim_bot_po4_3, mon_clim_bot_po4_3, mon_clim_bot_po4_3]

mon_clim_sur_do_1 = [mon_clim_sur_do_1, mon_clim_sur_do_1, mon_clim_sur_do_1]
mon_clim_sur_do_2 = [mon_clim_sur_do_2, mon_clim_sur_do_2, mon_clim_sur_do_2]
mon_clim_sur_do_3 = [mon_clim_sur_do_3, mon_clim_sur_do_3, mon_clim_sur_do_3]

mon_clim_bot_do_1 = [mon_clim_bot_do_1, mon_clim_bot_do_1, mon_clim_bot_do_1]
mon_clim_bot_do_2 = [mon_clim_bot_do_2, mon_clim_bot_do_2, mon_clim_bot_do_2]
mon_clim_bot_do_3 = [mon_clim_bot_do_3, mon_clim_bot_do_3, mon_clim_bot_do_3]

mon_clim_sur_si_1 = [mon_clim_sur_si_1, mon_clim_sur_si_1, mon_clim_sur_si_1]
mon_clim_sur_si_2 = [mon_clim_sur_si_2, mon_clim_sur_si_2, mon_clim_sur_si_2]
mon_clim_sur_si_3 = [mon_clim_sur_si_3, mon_clim_sur_si_3, mon_clim_sur_si_3]

mon_clim_bot_si_1 = [mon_clim_bot_si_1, mon_clim_bot_si_1, mon_clim_bot_si_1]
mon_clim_bot_si_2 = [mon_clim_bot_si_2, mon_clim_bot_si_2, mon_clim_bot_si_2]
mon_clim_bot_si_3 = [mon_clim_bot_si_3, mon_clim_bot_si_3, mon_clim_bot_si_3]

mon_clim_sur_temp_1 = [mon_clim_sur_temp_1, mon_clim_sur_temp_1, mon_clim_sur_temp_1]
mon_clim_sur_temp_2 = [mon_clim_sur_temp_2, mon_clim_sur_temp_2, mon_clim_sur_temp_2]
mon_clim_sur_temp_3 = [mon_clim_sur_temp_3, mon_clim_sur_temp_3, mon_clim_sur_temp_3]

mon_clim_bot_temp_1 = [mon_clim_bot_temp_1, mon_clim_bot_temp_1, mon_clim_bot_temp_1]
mon_clim_bot_temp_2 = [mon_clim_bot_temp_2, mon_clim_bot_temp_2, mon_clim_bot_temp_2]
mon_clim_bot_temp_3 = [mon_clim_bot_temp_3, mon_clim_bot_temp_3, mon_clim_bot_temp_3]

mon_clim_sur_salt_1 = [mon_clim_sur_salt_1, mon_clim_sur_salt_1, mon_clim_sur_salt_1]
mon_clim_sur_salt_2 = [mon_clim_sur_salt_2, mon_clim_sur_salt_2, mon_clim_sur_salt_2]
mon_clim_sur_salt_3 = [mon_clim_sur_salt_3, mon_clim_sur_salt_3, mon_clim_sur_salt_3]

mon_clim_bot_salt_1 = [mon_clim_bot_salt_1, mon_clim_bot_salt_1, mon_clim_bot_salt_1]
mon_clim_bot_salt_2 = [mon_clim_bot_salt_2, mon_clim_bot_salt_2, mon_clim_bot_salt_2]
mon_clim_bot_salt_3 = [mon_clim_bot_salt_3, mon_clim_bot_salt_3, mon_clim_bot_salt_3]

mon_clim_sur_ph_1 = [mon_clim_sur_ph_1, mon_clim_sur_ph_1, mon_clim_sur_ph_1]
mon_clim_sur_ph_2 = [mon_clim_sur_ph_2, mon_clim_sur_ph_2, mon_clim_sur_ph_2]
mon_clim_sur_ph_3 = [mon_clim_sur_ph_3, mon_clim_sur_ph_3, mon_clim_sur_ph_3]

mon_clim_bot_ph_1 = [mon_clim_bot_ph_1, mon_clim_bot_ph_1, mon_clim_bot_ph_1]
mon_clim_bot_ph_2 = [mon_clim_bot_ph_2, mon_clim_bot_ph_2, mon_clim_bot_ph_2]
mon_clim_bot_ph_3 = [mon_clim_bot_ph_3, mon_clim_bot_ph_3, mon_clim_bot_ph_3]
%% interp
t=1:length(mon_clim_sur_ph_1);
if length(t(~isnan(mon_clim_sur_ph_1))) > 0
t=1:length(mon_clim_sur_ph_1);
mon_clim_sur_ph_1(isnan(mon_clim_sur_ph_1)) = interp1(t(~isnan(mon_clim_sur_ph_1)),mon_clim_sur_ph_1(~isnan(mon_clim_sur_ph_1)),t(isnan(mon_clim_sur_ph_1)));
t=1:length(mon_clim_sur_ph_2);
mon_clim_sur_ph_2(isnan(mon_clim_sur_ph_2)) = interp1(t(~isnan(mon_clim_sur_ph_2)),mon_clim_sur_ph_2(~isnan(mon_clim_sur_ph_2)),t(isnan(mon_clim_sur_ph_2)));
t=1:length(mon_clim_sur_ph_3);
mon_clim_sur_ph_3(isnan(mon_clim_sur_ph_3)) = interp1(t(~isnan(mon_clim_sur_ph_3)),mon_clim_sur_ph_3(~isnan(mon_clim_sur_ph_3)),t(isnan(mon_clim_sur_ph_3)));


t=1:length(mon_clim_bot_ph_1);
mon_clim_bot_ph_1(isnan(mon_clim_bot_ph_1)) = interp1(t(~isnan(mon_clim_bot_ph_1)),mon_clim_bot_ph_1(~isnan(mon_clim_bot_ph_1)),t(isnan(mon_clim_bot_ph_1)));
t=1:length(mon_clim_bot_ph_2);
mon_clim_bot_ph_2(isnan(mon_clim_bot_ph_2)) = interp1(t(~isnan(mon_clim_bot_ph_2)),mon_clim_bot_ph_2(~isnan(mon_clim_bot_ph_2)),t(isnan(mon_clim_bot_ph_2)));
t=1:length(mon_clim_bot_ph_3);
mon_clim_bot_ph_3(isnan(mon_clim_bot_ph_3)) = interp1(t(~isnan(mon_clim_bot_ph_3)),mon_clim_bot_ph_3(~isnan(mon_clim_bot_ph_3)),t(isnan(mon_clim_bot_ph_3)));
end

t=1:length(mon_clim_sur_secchi_1);
mon_clim_sur_secchi_1(isnan(mon_clim_sur_secchi_1)) = interp1(t(~isnan(mon_clim_sur_secchi_1)),mon_clim_sur_secchi_1(~isnan(mon_clim_sur_secchi_1)),t(isnan(mon_clim_sur_secchi_1)));
t=1:length(mon_clim_sur_secchi_2);
mon_clim_sur_secchi_2(isnan(mon_clim_sur_secchi_2)) = interp1(t(~isnan(mon_clim_sur_secchi_2)),mon_clim_sur_secchi_2(~isnan(mon_clim_sur_secchi_2)),t(isnan(mon_clim_sur_secchi_2)));
t=1:length(mon_clim_sur_secchi_3);
mon_clim_sur_secchi_3(isnan(mon_clim_sur_secchi_3)) = interp1(t(~isnan(mon_clim_sur_secchi_3)),mon_clim_sur_secchi_3(~isnan(mon_clim_sur_secchi_3)),t(isnan(mon_clim_sur_secchi_3)));

t=1:length(mon_clim_bot_secchi_1);
mon_clim_bot_secchi_1(isnan(mon_clim_bot_secchi_1)) = interp1(t(~isnan(mon_clim_bot_secchi_1)),mon_clim_bot_secchi_1(~isnan(mon_clim_bot_secchi_1)),t(isnan(mon_clim_bot_secchi_1)));
t=1:length(mon_clim_bot_secchi_2);
mon_clim_bot_secchi_2(isnan(mon_clim_bot_secchi_2)) = interp1(t(~isnan(mon_clim_bot_secchi_2)),mon_clim_bot_secchi_2(~isnan(mon_clim_bot_secchi_2)),t(isnan(mon_clim_bot_secchi_2)));
t=1:length(mon_clim_bot_secchi_3);
mon_clim_bot_secchi_3(isnan(mon_clim_bot_secchi_3)) = interp1(t(~isnan(mon_clim_bot_secchi_3)),mon_clim_bot_secchi_3(~isnan(mon_clim_bot_secchi_3)),t(isnan(mon_clim_bot_secchi_3)));


t=1:length(mon_clim_sur_no2_1);
mon_clim_sur_no2_1(isnan(mon_clim_sur_no2_1)) = interp1(t(~isnan(mon_clim_sur_no2_1)),mon_clim_sur_no2_1(~isnan(mon_clim_sur_no2_1)),t(isnan(mon_clim_sur_no2_1)));
t=1:length(mon_clim_sur_no2_2);
mon_clim_sur_no2_2(isnan(mon_clim_sur_no2_2)) = interp1(t(~isnan(mon_clim_sur_no2_2)),mon_clim_sur_no2_2(~isnan(mon_clim_sur_no2_2)),t(isnan(mon_clim_sur_no2_2)));
t=1:length(mon_clim_sur_no2_3);
mon_clim_sur_no2_3(isnan(mon_clim_sur_no2_3)) = interp1(t(~isnan(mon_clim_sur_no2_3)),mon_clim_sur_no2_3(~isnan(mon_clim_sur_no2_3)),t(isnan(mon_clim_sur_no2_3)));

t=1:length(mon_clim_bot_no2_1);
mon_clim_bot_no2_1(isnan(mon_clim_bot_no2_1)) = interp1(t(~isnan(mon_clim_bot_no2_1)),mon_clim_bot_no2_1(~isnan(mon_clim_bot_no2_1)),t(isnan(mon_clim_bot_no2_1)));
t=1:length(mon_clim_bot_no2_2);
mon_clim_bot_no2_2(isnan(mon_clim_bot_no2_2)) = interp1(t(~isnan(mon_clim_bot_no2_2)),mon_clim_bot_no2_2(~isnan(mon_clim_bot_no2_2)),t(isnan(mon_clim_bot_no2_2)));
t=1:length(mon_clim_bot_no2_3);
mon_clim_bot_no2_3(isnan(mon_clim_bot_no2_3)) = interp1(t(~isnan(mon_clim_bot_no2_3)),mon_clim_bot_no2_3(~isnan(mon_clim_bot_no2_3)),t(isnan(mon_clim_bot_no2_3)));

t=1:length(mon_clim_sur_no3_1);
mon_clim_sur_no3_1(isnan(mon_clim_sur_no3_1)) = interp1(t(~isnan(mon_clim_sur_no3_1)),mon_clim_sur_no3_1(~isnan(mon_clim_sur_no3_1)),t(isnan(mon_clim_sur_no3_1)));
t=1:length(mon_clim_sur_no3_2);
mon_clim_sur_no3_2(isnan(mon_clim_sur_no3_2)) = interp1(t(~isnan(mon_clim_sur_no3_2)),mon_clim_sur_no3_2(~isnan(mon_clim_sur_no3_2)),t(isnan(mon_clim_sur_no3_2)));
t=1:length(mon_clim_sur_no3_3);
mon_clim_sur_no3_3(isnan(mon_clim_sur_no3_3)) = interp1(t(~isnan(mon_clim_sur_no3_3)),mon_clim_sur_no3_3(~isnan(mon_clim_sur_no3_3)),t(isnan(mon_clim_sur_no3_3)));

t=1:length(mon_clim_bot_no3_1);
mon_clim_bot_no3_1(isnan(mon_clim_bot_no3_1)) = interp1(t(~isnan(mon_clim_bot_no3_1)),mon_clim_bot_no3_1(~isnan(mon_clim_bot_no3_1)),t(isnan(mon_clim_bot_no3_1)));
t=1:length(mon_clim_bot_no3_2);
mon_clim_bot_no3_2(isnan(mon_clim_bot_no3_2)) = interp1(t(~isnan(mon_clim_bot_no3_2)),mon_clim_bot_no3_2(~isnan(mon_clim_bot_no3_2)),t(isnan(mon_clim_bot_no3_2)));
t=1:length(mon_clim_bot_no3_3);
mon_clim_bot_no3_3(isnan(mon_clim_bot_no3_3)) = interp1(t(~isnan(mon_clim_bot_no3_3)),mon_clim_bot_no3_3(~isnan(mon_clim_bot_no3_3)),t(isnan(mon_clim_bot_no3_3)));

t=1:length(mon_clim_sur_po4_1);
mon_clim_sur_po4_1(isnan(mon_clim_sur_po4_1)) = interp1(t(~isnan(mon_clim_sur_po4_1)),mon_clim_sur_po4_1(~isnan(mon_clim_sur_po4_1)),t(isnan(mon_clim_sur_po4_1)));
t=1:length(mon_clim_sur_po4_2);
mon_clim_sur_po4_2(isnan(mon_clim_sur_po4_2)) = interp1(t(~isnan(mon_clim_sur_po4_2)),mon_clim_sur_po4_2(~isnan(mon_clim_sur_po4_2)),t(isnan(mon_clim_sur_po4_2)));
t=1:length(mon_clim_sur_po4_3);
mon_clim_sur_po4_3(isnan(mon_clim_sur_po4_3)) = interp1(t(~isnan(mon_clim_sur_po4_3)),mon_clim_sur_po4_3(~isnan(mon_clim_sur_po4_3)),t(isnan(mon_clim_sur_po4_3)));

t=1:length(mon_clim_bot_po4_1);
mon_clim_bot_po4_1(isnan(mon_clim_bot_po4_1)) = interp1(t(~isnan(mon_clim_bot_po4_1)),mon_clim_bot_po4_1(~isnan(mon_clim_bot_po4_1)),t(isnan(mon_clim_bot_po4_1)));
t=1:length(mon_clim_bot_po4_2);
mon_clim_bot_po4_2(isnan(mon_clim_bot_po4_2)) = interp1(t(~isnan(mon_clim_bot_po4_2)),mon_clim_bot_po4_2(~isnan(mon_clim_bot_po4_2)),t(isnan(mon_clim_bot_po4_2)));
t=1:length(mon_clim_bot_po4_3);
mon_clim_bot_po4_3(isnan(mon_clim_bot_po4_3)) = interp1(t(~isnan(mon_clim_bot_po4_3)),mon_clim_bot_po4_3(~isnan(mon_clim_bot_po4_3)),t(isnan(mon_clim_bot_po4_3)));

t=1:length(mon_clim_sur_si_1);
mon_clim_sur_si_1(isnan(mon_clim_sur_si_1)) = interp1(t(~isnan(mon_clim_sur_si_1)),mon_clim_sur_si_1(~isnan(mon_clim_sur_si_1)),t(isnan(mon_clim_sur_si_1)));
t=1:length(mon_clim_sur_si_2);
mon_clim_sur_si_2(isnan(mon_clim_sur_si_2)) = interp1(t(~isnan(mon_clim_sur_si_2)),mon_clim_sur_si_2(~isnan(mon_clim_sur_si_2)),t(isnan(mon_clim_sur_si_2)));
t=1:length(mon_clim_sur_si_3);
mon_clim_sur_si_3(isnan(mon_clim_sur_si_3)) = interp1(t(~isnan(mon_clim_sur_si_3)),mon_clim_sur_si_3(~isnan(mon_clim_sur_si_3)),t(isnan(mon_clim_sur_si_3)));

t=1:length(mon_clim_bot_si_1);
mon_clim_bot_si_1(isnan(mon_clim_bot_si_1)) = interp1(t(~isnan(mon_clim_bot_si_1)),mon_clim_bot_si_1(~isnan(mon_clim_bot_si_1)),t(isnan(mon_clim_bot_si_1)));
t=1:length(mon_clim_bot_si_2);
mon_clim_bot_si_2(isnan(mon_clim_bot_si_2)) = interp1(t(~isnan(mon_clim_bot_si_2)),mon_clim_bot_si_2(~isnan(mon_clim_bot_si_2)),t(isnan(mon_clim_bot_si_2)));
t=1:length(mon_clim_bot_si_3);
mon_clim_bot_si_3(isnan(mon_clim_bot_si_3)) = interp1(t(~isnan(mon_clim_bot_si_3)),mon_clim_bot_si_3(~isnan(mon_clim_bot_si_3)),t(isnan(mon_clim_bot_si_3)));

t=1:length(mon_clim_sur_do_1);
mon_clim_sur_do_1(isnan(mon_clim_sur_do_1)) = interp1(t(~isnan(mon_clim_sur_do_1)),mon_clim_sur_do_1(~isnan(mon_clim_sur_do_1)),t(isnan(mon_clim_sur_do_1)));
t=1:length(mon_clim_sur_do_2);
mon_clim_sur_do_2(isnan(mon_clim_sur_do_2)) = interp1(t(~isnan(mon_clim_sur_do_2)),mon_clim_sur_do_2(~isnan(mon_clim_sur_do_2)),t(isnan(mon_clim_sur_do_2)));
t=1:length(mon_clim_sur_do_3);
mon_clim_sur_do_3(isnan(mon_clim_sur_do_3)) = interp1(t(~isnan(mon_clim_sur_do_3)),mon_clim_sur_do_3(~isnan(mon_clim_sur_do_3)),t(isnan(mon_clim_sur_do_3)));

t=1:length(mon_clim_bot_do_1);
mon_clim_bot_do_1(isnan(mon_clim_bot_do_1)) = interp1(t(~isnan(mon_clim_bot_do_1)),mon_clim_bot_do_1(~isnan(mon_clim_bot_do_1)),t(isnan(mon_clim_bot_do_1)));
t=1:length(mon_clim_bot_do_2);
mon_clim_bot_do_2(isnan(mon_clim_bot_do_2)) = interp1(t(~isnan(mon_clim_bot_do_2)),mon_clim_bot_do_2(~isnan(mon_clim_bot_do_2)),t(isnan(mon_clim_bot_do_2)));
t=1:length(mon_clim_bot_do_3);
mon_clim_bot_do_3(isnan(mon_clim_bot_do_3)) = interp1(t(~isnan(mon_clim_bot_do_3)),mon_clim_bot_do_3(~isnan(mon_clim_bot_do_3)),t(isnan(mon_clim_bot_do_3)));

t=1:length(mon_clim_sur_temp_1);
mon_clim_sur_temp_1(isnan(mon_clim_sur_temp_1)) = interp1(t(~isnan(mon_clim_sur_temp_1)),mon_clim_sur_temp_1(~isnan(mon_clim_sur_temp_1)),t(isnan(mon_clim_sur_temp_1)));
t=1:length(mon_clim_sur_temp_2);
mon_clim_sur_temp_2(isnan(mon_clim_sur_temp_2)) = interp1(t(~isnan(mon_clim_sur_temp_2)),mon_clim_sur_temp_2(~isnan(mon_clim_sur_temp_2)),t(isnan(mon_clim_sur_temp_2)));
t=1:length(mon_clim_sur_temp_3);
mon_clim_sur_temp_3(isnan(mon_clim_sur_temp_3)) = interp1(t(~isnan(mon_clim_sur_temp_3)),mon_clim_sur_temp_3(~isnan(mon_clim_sur_temp_3)),t(isnan(mon_clim_sur_temp_3)));

t=1:length(mon_clim_bot_temp_1);
mon_clim_bot_temp_1(isnan(mon_clim_bot_temp_1)) = interp1(t(~isnan(mon_clim_bot_temp_1)),mon_clim_bot_temp_1(~isnan(mon_clim_bot_temp_1)),t(isnan(mon_clim_bot_temp_1)));
t=1:length(mon_clim_bot_temp_2);
mon_clim_bot_temp_2(isnan(mon_clim_bot_temp_2)) = interp1(t(~isnan(mon_clim_bot_temp_2)),mon_clim_bot_temp_2(~isnan(mon_clim_bot_temp_2)),t(isnan(mon_clim_bot_temp_2)));
t=1:length(mon_clim_bot_temp_3);
mon_clim_bot_temp_3(isnan(mon_clim_bot_temp_3)) = interp1(t(~isnan(mon_clim_bot_temp_3)),mon_clim_bot_temp_3(~isnan(mon_clim_bot_temp_3)),t(isnan(mon_clim_bot_temp_3)));

t=1:length(mon_clim_sur_salt_1);
mon_clim_sur_salt_1(isnan(mon_clim_sur_salt_1)) = interp1(t(~isnan(mon_clim_sur_salt_1)),mon_clim_sur_salt_1(~isnan(mon_clim_sur_salt_1)),t(isnan(mon_clim_sur_salt_1)));
t=1:length(mon_clim_sur_salt_2);
mon_clim_sur_salt_2(isnan(mon_clim_sur_salt_2)) = interp1(t(~isnan(mon_clim_sur_salt_2)),mon_clim_sur_salt_2(~isnan(mon_clim_sur_salt_2)),t(isnan(mon_clim_sur_salt_2)));
t=1:length(mon_clim_sur_salt_3);
mon_clim_sur_salt_3(isnan(mon_clim_sur_salt_3)) = interp1(t(~isnan(mon_clim_sur_salt_3)),mon_clim_sur_salt_3(~isnan(mon_clim_sur_salt_3)),t(isnan(mon_clim_sur_salt_3)));

t=1:length(mon_clim_bot_salt_1);
mon_clim_bot_salt_1(isnan(mon_clim_bot_salt_1)) = interp1(t(~isnan(mon_clim_bot_salt_1)),mon_clim_bot_salt_1(~isnan(mon_clim_bot_salt_1)),t(isnan(mon_clim_bot_salt_1)));
t=1:length(mon_clim_bot_salt_2);
mon_clim_bot_salt_2(isnan(mon_clim_bot_salt_2)) = interp1(t(~isnan(mon_clim_bot_salt_2)),mon_clim_bot_salt_2(~isnan(mon_clim_bot_salt_2)),t(isnan(mon_clim_bot_salt_2)));
t=1:length(mon_clim_bot_salt_3);
mon_clim_bot_salt_3(isnan(mon_clim_bot_salt_3)) = interp1(t(~isnan(mon_clim_bot_salt_3)),mon_clim_bot_salt_3(~isnan(mon_clim_bot_salt_3)),t(isnan(mon_clim_bot_salt_3)));


%% cut edge
mon_clim_sur_secchi_1_in= mon_clim_sur_secchi_1(13:24); mon_clim_sur_secchi_2_in= mon_clim_sur_secchi_2(13:24); mon_clim_sur_secchi_3_in= mon_clim_sur_secchi_3(13:24);
mon_clim_bot_secchi_1_in= mon_clim_bot_secchi_1(13:24); mon_clim_bot_secchi_2_in= mon_clim_bot_secchi_2(13:24); mon_clim_bot_secchi_3_in= mon_clim_bot_secchi_3(13:24);

mon_clim_sur_no3_1_in= mon_clim_sur_no3_1(13:24); mon_clim_sur_no3_2_in= mon_clim_sur_no3_2(13:24); mon_clim_sur_no3_3_in= mon_clim_sur_no3_3(13:24);
mon_clim_bot_no3_1_in= mon_clim_bot_no3_1(13:24); mon_clim_bot_no3_2_in= mon_clim_bot_no3_2(13:24); mon_clim_bot_no3_3_in= mon_clim_bot_no3_3(13:24);

mon_clim_sur_si_1_in= mon_clim_sur_si_1(13:24); mon_clim_sur_si_2_in= mon_clim_sur_si_2(13:24); mon_clim_sur_si_3_in= mon_clim_sur_si_3(13:24);
mon_clim_bot_si_1_in= mon_clim_bot_si_1(13:24); mon_clim_bot_si_2_in= mon_clim_bot_si_2(13:24); mon_clim_bot_si_3_in= mon_clim_bot_si_3(13:24);

mon_clim_sur_do_1_in= mon_clim_sur_do_1(13:24); mon_clim_sur_do_2_in= mon_clim_sur_do_2(13:24); mon_clim_sur_do_3_in= mon_clim_sur_do_3(13:24);
mon_clim_bot_do_1_in= mon_clim_bot_do_1(13:24); mon_clim_bot_do_2_in= mon_clim_bot_do_2(13:24); mon_clim_bot_do_3_in= mon_clim_bot_do_3(13:24);

mon_clim_sur_temp_1_in= mon_clim_sur_temp_1(13:24); mon_clim_sur_temp_2_in= mon_clim_sur_temp_2(13:24); mon_clim_sur_temp_3_in= mon_clim_sur_temp_3(13:24);
mon_clim_bot_temp_1_in= mon_clim_bot_temp_1(13:24); mon_clim_bot_temp_2_in= mon_clim_bot_temp_2(13:24); mon_clim_bot_temp_3_in= mon_clim_bot_temp_3(13:24);

mon_clim_sur_salt_1_in= mon_clim_sur_salt_1(13:24); mon_clim_sur_salt_2_in= mon_clim_sur_salt_2(13:24); mon_clim_sur_salt_3_in= mon_clim_sur_salt_3(13:24);
mon_clim_bot_salt_1_in= mon_clim_bot_salt_1(13:24); mon_clim_bot_salt_2_in= mon_clim_bot_salt_2(13:24); mon_clim_bot_salt_3_in= mon_clim_bot_salt_3(13:24);

mon_clim_sur_ph_1_in= mon_clim_sur_ph_1(13:24); mon_clim_sur_ph_2_in= mon_clim_sur_ph_2(13:24); mon_clim_sur_ph_3_in= mon_clim_sur_ph_3(13:24);
mon_clim_bot_ph_1_in= mon_clim_bot_ph_1(13:24); mon_clim_bot_ph_2_in= mon_clim_bot_ph_2(13:24); mon_clim_bot_ph_3_in= mon_clim_bot_ph_3(13:24);

mon_clim_sur_no2_1_in= mon_clim_sur_no2_1(13:24); mon_clim_sur_no2_2_in= mon_clim_sur_no2_2(13:24); mon_clim_sur_no2_3_in= mon_clim_sur_no2_3(13:24);
mon_clim_bot_no2_1_in= mon_clim_bot_no2_1(13:24); mon_clim_bot_no2_2_in= mon_clim_bot_no2_2(13:24); mon_clim_bot_no2_3_in= mon_clim_bot_no2_3(13:24);

mon_clim_sur_po4_1_in= mon_clim_sur_po4_1(13:24); mon_clim_sur_po4_2_in= mon_clim_sur_po4_2(13:24); mon_clim_sur_po4_3_in= mon_clim_sur_po4_3(13:24);
mon_clim_bot_po4_1_in= mon_clim_bot_po4_1(13:24); mon_clim_bot_po4_2_in= mon_clim_bot_po4_2(13:24); mon_clim_bot_po4_3_in= mon_clim_bot_po4_3(13:24);

save('kodc_input_fix_to06to15_3sig(3regime)_20501_v7.mat');


figure; hold on;
plot(mon_clim_sur_salt_1_in,'r','linew',2); plot(mon_clim_sur_salt_2_in,'g','linew',2); 
plot(mon_clim_sur_salt_3_in,'b','linew',2);
title('205-01 NIFS salinity')
xlim([1 12]);
xticks(1:12);
grid on; set(gca,'fontsize',15); xlabel('month');

salt_bot(salt_bot > 99)=NaN;

part_sur_jan=salt_sur(1:12:end);
part_bot_jan=salt_bot(1:12:end);
part_sur_mar=salt_sur(3:12:end);
part_bot_mar=salt_bot(3:12:end);
part_sur_feb=salt_sur(2:12:end);
part_bot_feb=salt_bot(2:12:end);
nandx=find(isnan(part_sur_feb)==1);
nandx_bot=find(isnan(part_bot_feb)==1);

part_sur_jan(nandx)
part_sur_mar(nandx)
part_bot_jan(nandx)
part_bot_mar(nandx)

merg_sur=part_sur_feb;
merg_bot=part_bot_feb;
merg_sur(nandx(1))=part_sur_jan(nandx(1));
merg_bot(nandx(1))=part_bot_jan(nandx(1));
merg_sur(nandx(2:3))=part_sur_mar(nandx(2:3));
merg_bot(nandx(2:3))=part_bot_mar(nandx(2:3));

figure; hold on;
plot(salt_sur(2:12:end),'ko');
plot(salt_bot(2:12:end),'ro');
plot(salt_sur(2:12:end),'k');
plot(salt_bot(2:12:end),'r');
xticks(1:41);
xticklabels(1980:2020);
xtickangle(45)
grid on; xlim([1 inf]); title('205-01 NIFS salinity feb')
legend('surf.','bot.');

figure; hold on;
plot(merg_sur,'ko');
plot(merg_bot,'ro');
plot(merg_sur,'k');
plot(merg_bot,'r');
xticks(1:41);
xticklabels(1980:2020);
xtickangle(45)
grid on; xlim([1 inf]); title('205-01 NIFS salinity feb (jan, mar)')
plot(nandx,merg_sur(nandx),'ok','MarkerFaceColor','k');
plot(nandx,merg_bot(nandx),'or','MarkerFaceColor','r');
legend('surf.','bot.');




