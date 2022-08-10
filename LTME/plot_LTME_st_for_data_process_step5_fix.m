close all; clear; clc;   % -v3
cd D:\장기생태\observation\관측\장기생태프로젝트\남해_DB
% % station name tag
% % (for preventing the miss picking using tag from 1 col's name on excel) 

% port
for i = 1:9
name_tag{i} = ['광양_' num2str(i)] 
end

% combining the tag and outter point excluding 
size_tag = length(name_tag);


%% pick the row on the excel which has same name with tag
[raw txt]=xlsread('표준입력시스템 DB_해양환경_수질환경_광양_2013_fix.xlsx','표본_일반수질항목','');
txt_matc_p = txt(4:end,1); % name list
txt_date_p = txt(4:end,2); % date list
txt_lat_p = txt(4:end,4); % name list
txt_lon_p = txt(4:end,5); % date list
temp_p = raw(1:end,3); 
salt_p = raw(1:end,4); 
no3_p = raw(1:end,10); 
nh4_p = raw(1:end,11); 
dep_p = raw(1:end,1); 



%% pick matched name with tag
% %port
% for i = 1:length(name_tag)
%    if  sum(strcmp(name_tag{i}, txt_matc_p)) ~= 0
%        indx{i} = find([strcmp(name_tag{i}, txt_matc_p)] == 1)     
%    end
% end

%merge
clearvars indx
for i = 1:length(name_tag)
   if  sum(strcmp(name_tag{i}, txt_matc_p)) ~= 0
       indx{i} = find([strcmp(name_tag{i}, txt_matc_p)] == 1)
       date_st{i} = txt_date_p(indx{i})
   end
end

% txt_date_p(indx{1})
% for i = 1:length
% surf_dx(i,:) = find(dep_p(indx{i}) == 0);
% end

for i = 1:length(name_tag)
    depth{i} = dep_p(indx{i});
    temp{i} = temp_p(indx{i});
    salt{i} = salt_p(indx{i});
    nh4{i} = nh4_p(indx{i});
    no3{i} = no3_p(indx{i});
end


% no3_sur_p(indx{9})
dep_p(indx{9})


%% make date to be 'yymm' form
for i = 1:length(name_tag)
clearvars temp_d temp_ymd temp_ymd_c temp_yymmdd temp_yymmdd_c
temp_d = char(date_st{i});
temp_yymmdd=temp_d;

for j = 1:length(temp_yymmdd)
    temp_yymmdd_c{j} = temp_yymmdd(j,:);          
end
date_yymmdd{i,1} = temp_yymmdd_c;
end

%% make 1997 to 2018 'yymm' form
% make 1980~present
k=0
for i = 2013:2013
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

k=0; m=0;
for i = 1:1
    l=0
        ref_yy(i,:)=[num2str(i+2012)];
    for n = 1:12
        m = m+1;
        ref_yymm{m}=[num2str(i+2012) '-' num2str(n,'%02d')];
    for j = 1:eom_d(i,n)
        l=l+1; % index for 1~31(or 30 or 29)
        k=k+1; % index for 1~366
        ref_yymmdd{k}=[num2str(i+2012) '-' num2str(n,'%02d') '-'  num2str(j,'%02d')];
        ref_mmdd{k}=[num2str(n,'%02d') '-'  num2str(j,'%02d')];
        ref_mm{k}=[num2str(n,'%02d')];
    end
    end
end

clearvars indx_date_s indx_date_b
indx_date_s = cell(1,length(ref_yymmdd));
indx_date_b = cell(1,length(ref_yymmdd));
% matched date 'yymm' form on surf & bottom
for j = 1:length(name_tag) % st. axis
    for i = 1:length(ref_yymmdd) % date axis
       if  sum(strcmp(ref_yymmdd{i}, date_st{j})) ~= 0
           clearvars indx_date
           indx_date = find([strcmp(ref_yymmdd{i}, date_yymmdd{j})] == 1); 
           indx_date_s{j,i} = indx_date(find(depth{j}(indx_date) == 0)); 
           indx_date_b{j,i} = indx_date(find(depth{j}(indx_date) == max(depth{j}(indx_date))));
       end
    end 
end

% make climate
%temp
for i = 1:size(indx_date_s,1)
    clearvars temp_k
    temp_k = temp{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_k(indx_date_s{i,j}))) ~= 0     
        temp_sur(i,j) = nanmean(temp_k(indx_date_s{i,j})); % matched date surf. (if more then 1 will be averaged)
    elseif sum(size(temp_k(indx_date_s{i,j}))) == 0  % matched date bot. (if more then 1 will be averaged)
        temp_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_k
    temp_k = temp{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_k(indx_date_b{i,j}))) ~= 0     
        temp_bot(i,j) = nanmean(temp_k(indx_date_b{i,j}));
    elseif sum(size(temp_k(indx_date_b{i,j}))) == 0 
        temp_bot(i,j) = NaN;
    end
    end
end

%salt
for i = 1:size(indx_date_s,1)
    clearvars temp_k
    temp_k = salt{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_k(indx_date_s{i,j}))) ~= 0     
        salt_sur(i,j) = nanmean(temp_k(indx_date_s{i,j}));
    elseif sum(size(temp_k(indx_date_s{i,j}))) == 0 
        salt_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_k
    temp_k = salt{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_k(indx_date_b{i,j}))) ~= 0     
        salt_bot(i,j) = nanmean(temp_k(indx_date_b{i,j}));
    elseif sum(size(temp_k(indx_date_b{i,j}))) == 0 
        salt_bot(i,j) = NaN;
    end
    end
end

%nh4
for i = 1:size(indx_date_s,1)
    clearvars temp_k
    temp_k = nh4{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_k(indx_date_s{i,j}))) ~= 0     
        nh4_sur(i,j) = nanmean(temp_k(indx_date_s{i,j}));
    elseif sum(size(temp_k(indx_date_s{i,j}))) == 0 
        nh4_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_k
    temp_k = nh4{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_k(indx_date_b{i,j}))) ~= 0     
        nh4_bot(i,j) = nanmean(temp_k(indx_date_b{i,j}));
    elseif sum(size(temp_k(indx_date_b{i,j}))) == 0 
        nh4_bot(i,j) = NaN;
    end
    end
end


%no3
for i = 1:size(indx_date_s,1)
    clearvars temp_k
    temp_k = no3{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_k(indx_date_s{i,j}))) ~= 0     
        no3_sur(i,j) = nanmean(temp_k(indx_date_s{i,j}));
    elseif sum(size(temp_k(indx_date_s{i,j}))) == 0 
        no3_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_k
    temp_k = no3{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_k(indx_date_b{i,j}))) ~= 0     
        no3_bot(i,j) = nanmean(temp_k(indx_date_b{i,j}));
    elseif sum(size(temp_k(indx_date_b{i,j}))) == 0 
        no3_bot(i,j) = NaN;
    end
    end
end

t_s_13 = temp_sur; s_s_13 = salt_sur; nh4_s_13 = nh4_sur; no3_s_13 = no3_sur;
t_b_13 = temp_bot; s_b_13 = salt_bot; nh4_b_13 = nh4_bot; no3_b_13 = no3_bot;
clearvars -except '*_s_13' '*_b_13'

%% 2015
% close all; clear; clc;   % -v3
present_yr = 2015;
cd D:\장기생태\observation\관측\장기생태프로젝트\남해_DB
% % station name tag
% % (for preventing the miss picking using tag from 1 col's name on excel) 

% port
for i = 1:9
name_tag{i} = ['광양_' num2str(i)] 
end

% combining the tag and outter point excluding 
size_tag = length(name_tag);

%% pick the row on the excel which has same name with tag
[raw txt]=xlsread('표준입력시스템 DB_해양환경_수질환경_광양_2015.xlsx','표본_일반수질항목','');
txt_matc_p = txt(4:end,1); % name list
txt_date_p = txt(4:end,2); % date list
txt_lat_p = txt(4:end,4); % name list
txt_lon_p = txt(4:end,5); % date list
temp_p = raw(1:end,6); 
salt_p = raw(1:end,7); 
no3_p = raw(1:end,13); 
nh4_p = raw(1:end,14); 
chl_p = raw(1:end,15); 
dep_p = raw(1:end,4); 

%% pick matched name with tag
% %port
% for i = 1:length(name_tag)
%    if  sum(strcmp(name_tag{i}, txt_matc_p)) ~= 0
%        indx{i} = find([strcmp(name_tag{i}, txt_matc_p)] == 1)     
%    end
% end

%merge
clearvars indx
for i = 1:length(name_tag)
   if  sum(strcmp(name_tag{i}, txt_matc_p)) ~= 0
       indx{i} = find([strcmp(name_tag{i}, txt_matc_p)] == 1)
       date_st{i} = txt_date_p(indx{i})
   end
end

% txt_date_p(indx{1})
% for i = 1:length
% surf_dx(i,:) = find(dep_p(indx{i}) == 0);
% end

for i = 1:length(name_tag)
    depth{i} = dep_p(indx{i});
    temp{i} = temp_p(indx{i});
    salt{i} = salt_p(indx{i});
    nh4{i} = nh4_p(indx{i});
    no3{i} = no3_p(indx{i});
    chl{i} = chl_p(indx{i});
end

% no3_sur_p(indx{9})
dep_p(indx{9})

%% make date to be 'yymm' form
for i = 1:length(name_tag)
clearvars temp_d temp_ymd temp_ymd_c temp_yymmdd temp_yymmdd_c
temp_d = char(date_st{i});
temp_yymmdd=temp_d;

for j = 1:size(temp_yymmdd,1)
    temp_yymmdd_c{j} = temp_yymmdd(j,:);          
end
date_yymmdd{i,1} = temp_yymmdd_c;
end

%% make 1997 to 2018 'yymm' form
% make 1980~present
k=0
for i = present_yr:present_yr
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

k=0; m=0;
for i = 1:1
    l=0
        ref_yy(i,:)=[num2str(i+present_yr-1)];
    for n = 1:12
        m = m+1;
        ref_yymm{m}=[num2str(i+present_yr-1) '-' num2str(n,'%02d')];
    for j = 1:eom_d(i,n)
        l=l+1; % index for 1~31(or 30 or 29)
        k=k+1; % index for 1~366
        ref_yymmdd{k}=[num2str(i+present_yr-1) '-' num2str(n,'%02d') '-'  num2str(j,'%02d')];
        ref_mmdd{k}=[num2str(n,'%02d') '-'  num2str(j,'%02d')];
        ref_mm{k}=[num2str(n,'%02d')];
    end
    end
end

clearvars indx_date_s indx_date_b
indx_date_s = cell(1,length(ref_yymmdd));
indx_date_b = cell(1,length(ref_yymmdd));
% matched date 'yymm' form on surf & bottom
for j = 1:length(name_tag) % st. axis
    for i = 1:length(ref_yymmdd) % date axis
       if  sum(strcmp(ref_yymmdd{i}, date_st{j})) ~= 0
           clearvars indx_date
           indx_date = find([strcmp(ref_yymmdd{i}, date_yymmdd{j})] == 1); 
           indx_date_s{j,i} = indx_date(find(depth{j}(indx_date) == 0)); 
           indx_date_b{j,i} = indx_date(find(depth{j}(indx_date) == max(depth{j}(indx_date))));
       end
    end 
end

% make climate
%temp
for i = 1:size(indx_date_s,1)
    clearvars temp_k
    temp_k = temp{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_k(indx_date_s{i,j}))) ~= 0     
        temp_sur(i,j) = nanmean(temp_k(indx_date_s{i,j})); % matched date surf. (if more then 1 will be averaged)
    elseif sum(size(temp_k(indx_date_s{i,j}))) == 0  % matched date bot. (if more then 1 will be averaged)
        temp_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_k
    temp_k = temp{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_k(indx_date_b{i,j}))) ~= 0     
        temp_bot(i,j) = nanmean(temp_k(indx_date_b{i,j}));
    elseif sum(size(temp_k(indx_date_b{i,j}))) == 0 
        temp_bot(i,j) = NaN;
    end
    end
end

%salt
for i = 1:size(indx_date_s,1)
    clearvars temp_k
    temp_k = salt{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_k(indx_date_s{i,j}))) ~= 0     
        salt_sur(i,j) = nanmean(temp_k(indx_date_s{i,j}));
    elseif sum(size(temp_k(indx_date_s{i,j}))) == 0 
        salt_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_k
    temp_k = salt{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_k(indx_date_b{i,j}))) ~= 0     
        salt_bot(i,j) = nanmean(temp_k(indx_date_b{i,j}));
    elseif sum(size(temp_k(indx_date_b{i,j}))) == 0 
        salt_bot(i,j) = NaN;
    end
    end
end

%nh4
for i = 1:size(indx_date_s,1)
    clearvars temp_k
    temp_k = nh4{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_k(indx_date_s{i,j}))) ~= 0     
        nh4_sur(i,j) = nanmean(temp_k(indx_date_s{i,j}));
    elseif sum(size(temp_k(indx_date_s{i,j}))) == 0 
        nh4_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_k
    temp_k = nh4{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_k(indx_date_b{i,j}))) ~= 0     
        nh4_bot(i,j) = nanmean(temp_k(indx_date_b{i,j}));
    elseif sum(size(temp_k(indx_date_b{i,j}))) == 0 
        nh4_bot(i,j) = NaN;
    end
    end
end


%no3
for i = 1:size(indx_date_s,1)
    clearvars temp_k
    temp_k = no3{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_k(indx_date_s{i,j}))) ~= 0     
        no3_sur(i,j) = nanmean(temp_k(indx_date_s{i,j}));
    elseif sum(size(temp_k(indx_date_s{i,j}))) == 0 
        no3_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_k
    temp_k = no3{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_k(indx_date_b{i,j}))) ~= 0     
        no3_bot(i,j) = nanmean(temp_k(indx_date_b{i,j}));
    elseif sum(size(temp_k(indx_date_b{i,j}))) == 0 
        no3_bot(i,j) = NaN;
    end
    end
end

%chl
for i = 1:size(indx_date_s,1)
    clearvars temp_k
    temp_k = chl{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_k(indx_date_s{i,j}))) ~= 0     
        chl_sur(i,j) = nanmean(temp_k(indx_date_s{i,j}));
    elseif sum(size(temp_k(indx_date_s{i,j}))) == 0 
        chl_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_k
    temp_k = chl{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_k(indx_date_b{i,j}))) ~= 0     
        chl_bot(i,j) = nanmean(temp_k(indx_date_b{i,j}));
    elseif sum(size(temp_k(indx_date_b{i,j}))) == 0 
        chl_bot(i,j) = NaN;
    end
    end
end

t_s_15 = temp_sur; s_s_15 = salt_sur; nh4_s_15 = nh4_sur; no3_s_15 = no3_sur; chl_s_15 = chl_sur;
t_b_15 = temp_bot; s_b_15 = salt_bot; nh4_b_15 = nh4_bot; no3_b_15 = no3_bot; chl_b_15 = chl_bot;
clearvars -except '*_s_15' '*_b_15' '*_s_13' '*_b_13'


%% 2016
% close all; clear; clc;   % -v3
present_yr = 2016;
cd D:\장기생태\observation\관측\장기생태프로젝트\남해_DB
% % station name tag
% % (for preventing the miss picking using tag from 1 col's name on excel) 

% port
for i = 1:9
name_tag{i} = ['광양_' num2str(i)] 
end

% combining the tag and outter point excluding 
size_tag = length(name_tag);

%% pick the row on the excel which has same name with tag
[raw txt]=xlsread('표준입력시스템 DB_해양환경_수질환경_광양_2016.xlsx','표본_일반수질항목','');
txt_matc_p = txt(4:end,1); % name list
txt_date_p = txt(4:end,2); % date list
txt_lat_p = txt(4:end,4); % name list
txt_lon_p = txt(4:end,5); % date list
temp_p = raw(1:end,6); 
salt_p = raw(1:end,7); 
no3_p = raw(1:end,13); 
nh4_p = raw(1:end,14); 
chl_p = raw(1:end,15); 
dep_p = raw(1:end,4); 

%% pick matched name with tag
% %port
% for i = 1:length(name_tag)
%    if  sum(strcmp(name_tag{i}, txt_matc_p)) ~= 0
%        indx{i} = find([strcmp(name_tag{i}, txt_matc_p)] == 1)     
%    end
% end

%merge
clearvars indx
for i = 1:length(name_tag)
   if  sum(strcmp(name_tag{i}, txt_matc_p)) ~= 0
       indx{i} = find([strcmp(name_tag{i}, txt_matc_p)] == 1)
       date_st{i} = txt_date_p(indx{i})
   end
end

% txt_date_p(indx{1})
% for i = 1:length
% surf_dx(i,:) = find(dep_p(indx{i}) == 0);
% end

for i = 1:length(name_tag)
    depth{i} = dep_p(indx{i});
    temp{i} = temp_p(indx{i});
    salt{i} = salt_p(indx{i});
    nh4{i} = nh4_p(indx{i});
    no3{i} = no3_p(indx{i});
    chl{i} = chl_p(indx{i});
end

% no3_sur_p(indx{9})
dep_p(indx{9})

%% make date to be 'yymm' form
for i = 1:length(name_tag)
clearvars temp_d temp_ymd temp_ymd_c temp_yymmdd temp_yymmdd_c
temp_d = char(date_st{i});
temp_yymmdd=temp_d;

for j = 1:size(temp_yymmdd,1)
    temp_yymmdd_c{j} = temp_yymmdd(j,:);          
end
date_yymmdd{i,1} = temp_yymmdd_c;
end

%% make 1997 to 2018 'yymm' form
% make 1980~present
k=0
for i = present_yr:present_yr
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

k=0; m=0;
for i = 1:1
    l=0
        ref_yy(i,:)=[num2str(i+present_yr-1)];
    for n = 1:12
        m = m+1;
        ref_yymm{m}=[num2str(i+present_yr-1) '-' num2str(n,'%02d')];
    for j = 1:eom_d(i,n)
        l=l+1; % index for 1~31(or 30 or 29)
        k=k+1; % index for 1~366
        ref_yymmdd{k}=[num2str(i+present_yr-1) '-' num2str(n,'%02d') '-'  num2str(j,'%02d')];
        ref_mmdd{k}=[num2str(n,'%02d') '-'  num2str(j,'%02d')];
        ref_mm{k}=[num2str(n,'%02d')];
    end
    end
end

clearvars indx_date_s indx_date_b
indx_date_s = cell(1,length(ref_yymmdd));
indx_date_b = cell(1,length(ref_yymmdd));
% matched date 'yymm' form on surf & bottom
for j = 1:length(name_tag) % st. axis
    for i = 1:length(ref_yymmdd) % date axis
       if  sum(strcmp(ref_yymmdd{i}, date_st{j})) ~= 0
           clearvars indx_date
           indx_date = find([strcmp(ref_yymmdd{i}, date_yymmdd{j})] == 1); 
           indx_date_s{j,i} = indx_date(find(depth{j}(indx_date) == 0)); 
           indx_date_b{j,i} = indx_date(find(depth{j}(indx_date) == max(depth{j}(indx_date))));
       end
    end 
end

% make climate
%temp
for i = 1:size(indx_date_s,1)
    clearvars temp_k
    temp_k = temp{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_k(indx_date_s{i,j}))) ~= 0     
        temp_sur(i,j) = nanmean(temp_k(indx_date_s{i,j})); % matched date surf. (if more then 1 will be averaged)
    elseif sum(size(temp_k(indx_date_s{i,j}))) == 0  % matched date bot. (if more then 1 will be averaged)
        temp_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_k
    temp_k = temp{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_k(indx_date_b{i,j}))) ~= 0     
        temp_bot(i,j) = nanmean(temp_k(indx_date_b{i,j}));
    elseif sum(size(temp_k(indx_date_b{i,j}))) == 0 
        temp_bot(i,j) = NaN;
    end
    end
end

%salt
for i = 1:size(indx_date_s,1)
    clearvars temp_k
    temp_k = salt{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_k(indx_date_s{i,j}))) ~= 0     
        salt_sur(i,j) = nanmean(temp_k(indx_date_s{i,j}));
    elseif sum(size(temp_k(indx_date_s{i,j}))) == 0 
        salt_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_k
    temp_k = salt{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_k(indx_date_b{i,j}))) ~= 0     
        salt_bot(i,j) = nanmean(temp_k(indx_date_b{i,j}));
    elseif sum(size(temp_k(indx_date_b{i,j}))) == 0 
        salt_bot(i,j) = NaN;
    end
    end
end

%nh4
for i = 1:size(indx_date_s,1)
    clearvars temp_k
    temp_k = nh4{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_k(indx_date_s{i,j}))) ~= 0     
        nh4_sur(i,j) = nanmean(temp_k(indx_date_s{i,j}));
    elseif sum(size(temp_k(indx_date_s{i,j}))) == 0 
        nh4_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_k
    temp_k = nh4{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_k(indx_date_b{i,j}))) ~= 0     
        nh4_bot(i,j) = nanmean(temp_k(indx_date_b{i,j}));
    elseif sum(size(temp_k(indx_date_b{i,j}))) == 0 
        nh4_bot(i,j) = NaN;
    end
    end
end


%no3
for i = 1:size(indx_date_s,1)
    clearvars temp_k
    temp_k = no3{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_k(indx_date_s{i,j}))) ~= 0     
        no3_sur(i,j) = nanmean(temp_k(indx_date_s{i,j}));
    elseif sum(size(temp_k(indx_date_s{i,j}))) == 0 
        no3_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_k
    temp_k = no3{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_k(indx_date_b{i,j}))) ~= 0     
        no3_bot(i,j) = nanmean(temp_k(indx_date_b{i,j}));
    elseif sum(size(temp_k(indx_date_b{i,j}))) == 0 
        no3_bot(i,j) = NaN;
    end
    end
end

%chl
for i = 1:size(indx_date_s,1)
    clearvars temp_k
    temp_k = chl{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_k(indx_date_s{i,j}))) ~= 0     
        chl_sur(i,j) = nanmean(temp_k(indx_date_s{i,j}));
    elseif sum(size(temp_k(indx_date_s{i,j}))) == 0 
        chl_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_k
    temp_k = chl{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_k(indx_date_b{i,j}))) ~= 0     
        chl_bot(i,j) = nanmean(temp_k(indx_date_b{i,j}));
    elseif sum(size(temp_k(indx_date_b{i,j}))) == 0 
        chl_bot(i,j) = NaN;
    end
    end
end

t_s_16 = temp_sur; s_s_16 = salt_sur; nh4_s_16 = nh4_sur; no3_s_16 = no3_sur; chl_s_16 = chl_sur;
t_b_16 = temp_bot; s_b_16 = salt_bot; nh4_b_16 = nh4_bot; no3_b_16 = no3_bot; chl_b_16 = chl_bot;
clearvars -except '*_s_16' '*_b_16' '*_s_15' '*_b_15' '*_s_13' '*_b_13'


%% 2017
% close all; clear; clc;   % -v3
present_yr = 2017;
cd D:\장기생태\observation\관측\장기생태프로젝트\남해_DB
% % station name tag
% % (for preventing the miss picking using tag from 1 col's name on excel) 

% port
for i = 1:9
name_tag{i} = ['광양_' num2str(i)] 
end

% combining the tag and outter point excluding 
size_tag = length(name_tag);

%% pick the row on the excel which has same name with tag
[raw txt]=xlsread('표준입력시스템 DB_해양환경_수질환경_광양_2017.xlsx','표본_일반수질항목','');
txt_matc_p = txt(4:end,1); % name list
txt_date_p = txt(4:end,2); % date list
txt_lat_p = txt(4:end,4); % name list
txt_lon_p = txt(4:end,5); % date list
temp_p = raw(1:end,6); 
salt_p = raw(1:end,7); 
no3_p = raw(1:end,13); 
nh4_p = raw(1:end,14); 
% chl_p = raw(1:end,15);  no-data
dep_p = raw(1:end,4); 

%% pick matched name with tag
% %port
% for i = 1:length(name_tag)
%    if  sum(strcmp(name_tag{i}, txt_matc_p)) ~= 0
%        indx{i} = find([strcmp(name_tag{i}, txt_matc_p)] == 1)     
%    end
% end

%merge
clearvars indx
for i = 1:length(name_tag)
   if  sum(strcmp(name_tag{i}, txt_matc_p)) ~= 0
       indx{i} = find([strcmp(name_tag{i}, txt_matc_p)] == 1)
       date_st{i} = txt_date_p(indx{i})
   end
end

% txt_date_p(indx{1})
% for i = 1:length
% surf_dx(i,:) = find(dep_p(indx{i}) == 0);
% end

for i = 1:length(name_tag)
    depth{i} = dep_p(indx{i});
    temp{i} = temp_p(indx{i});
    salt{i} = salt_p(indx{i});
    nh4{i} = nh4_p(indx{i});
    no3{i} = no3_p(indx{i});
%     chl{i} = chl_p(indx{i});
end

% no3_sur_p(indx{9})
dep_p(indx{9})

%% make date to be 'yymm' form
for i = 1:length(name_tag)
clearvars temp_d temp_ymd temp_ymd_c temp_yymmdd temp_yymmdd_c
temp_d = char(date_st{i});
temp_yymmdd=temp_d;

for j = 1:size(temp_yymmdd,1)
    temp_yymmdd_c{j} = temp_yymmdd(j,:);          
end
date_yymmdd{i,1} = temp_yymmdd_c;
end

%% make 1997 to 2018 'yymm' form
% make 1980~present
k=0
for i = present_yr:present_yr
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

k=0; m=0;
for i = 1:1
    l=0
        ref_yy(i,:)=[num2str(i+present_yr-1)];
    for n = 1:12
        m = m+1;
        ref_yymm{m}=[num2str(i+present_yr-1) '-' num2str(n,'%02d')];
    for j = 1:eom_d(i,n)
        l=l+1; % index for 1~31(or 30 or 29)
        k=k+1; % index for 1~366
        ref_yymmdd{k}=[num2str(i+present_yr-1) '-' num2str(n,'%02d') '-'  num2str(j,'%02d')];
        ref_mmdd{k}=[num2str(n,'%02d') '-'  num2str(j,'%02d')];
        ref_mm{k}=[num2str(n,'%02d')];
    end
    end
end

clearvars indx_date_s indx_date_b
indx_date_s = cell(1,length(ref_yymmdd));
indx_date_b = cell(1,length(ref_yymmdd));
% matched date 'yymm' form on surf & bottom
for j = 1:length(name_tag) % st. axis
    for i = 1:length(ref_yymmdd) % date axis
       if  sum(strcmp(ref_yymmdd{i}, date_st{j})) ~= 0
           clearvars indx_date
           indx_date = find([strcmp(ref_yymmdd{i}, date_yymmdd{j})] == 1); 
           indx_date_s{j,i} = indx_date(find(depth{j}(indx_date) == 0)); 
           indx_date_b{j,i} = indx_date(find(depth{j}(indx_date) == max(depth{j}(indx_date))));
       end
    end 
end

% make climate
%temp
for i = 1:size(indx_date_s,1)
    clearvars temp_k
    temp_k = temp{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_k(indx_date_s{i,j}))) ~= 0     
        temp_sur(i,j) = nanmean(temp_k(indx_date_s{i,j})); % matched date surf. (if more then 1 will be averaged)
    elseif sum(size(temp_k(indx_date_s{i,j}))) == 0  % matched date bot. (if more then 1 will be averaged)
        temp_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_k
    temp_k = temp{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_k(indx_date_b{i,j}))) ~= 0     
        temp_bot(i,j) = nanmean(temp_k(indx_date_b{i,j}));
    elseif sum(size(temp_k(indx_date_b{i,j}))) == 0 
        temp_bot(i,j) = NaN;
    end
    end
end

%salt
for i = 1:size(indx_date_s,1)
    clearvars temp_k
    temp_k = salt{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_k(indx_date_s{i,j}))) ~= 0     
        salt_sur(i,j) = nanmean(temp_k(indx_date_s{i,j}));
    elseif sum(size(temp_k(indx_date_s{i,j}))) == 0 
        salt_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_k
    temp_k = salt{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_k(indx_date_b{i,j}))) ~= 0     
        salt_bot(i,j) = nanmean(temp_k(indx_date_b{i,j}));
    elseif sum(size(temp_k(indx_date_b{i,j}))) == 0 
        salt_bot(i,j) = NaN;
    end
    end
end

%nh4
for i = 1:size(indx_date_s,1)
    clearvars temp_k
    temp_k = nh4{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_k(indx_date_s{i,j}))) ~= 0     
        nh4_sur(i,j) = nanmean(temp_k(indx_date_s{i,j}));
    elseif sum(size(temp_k(indx_date_s{i,j}))) == 0 
        nh4_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_k
    temp_k = nh4{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_k(indx_date_b{i,j}))) ~= 0     
        nh4_bot(i,j) = nanmean(temp_k(indx_date_b{i,j}));
    elseif sum(size(temp_k(indx_date_b{i,j}))) == 0 
        nh4_bot(i,j) = NaN;
    end
    end
end


%no3
for i = 1:size(indx_date_s,1)
    clearvars temp_k
    temp_k = no3{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_k(indx_date_s{i,j}))) ~= 0     
        no3_sur(i,j) = nanmean(temp_k(indx_date_s{i,j}));
    elseif sum(size(temp_k(indx_date_s{i,j}))) == 0 
        no3_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_k
    temp_k = no3{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_k(indx_date_b{i,j}))) ~= 0     
        no3_bot(i,j) = nanmean(temp_k(indx_date_b{i,j}));
    elseif sum(size(temp_k(indx_date_b{i,j}))) == 0 
        no3_bot(i,j) = NaN;
    end
    end
end

% %chl
% for i = 1:size(indx_date_s,1)
%     clearvars temp_k
%     temp_k = chl{i};
%     for j = 1:size(indx_date_s,2)
%     if sum(size(temp_k(indx_date_s{i,j}))) ~= 0     
%         chl_sur(i,j) = nanmean(temp_k(indx_date_s{i,j}));
%     elseif sum(size(temp_k(indx_date_s{i,j}))) == 0 
%         chl_sur(i,j) = NaN;
%     end
%     end
% end
% 
% for i = 1:size(indx_date_b,1)
%     clearvars temp_k
%     temp_k = chl{i};
%     for j = 1:size(indx_date_b,2)
%     if sum(size(temp_k(indx_date_b{i,j}))) ~= 0     
%         chl_bot(i,j) = nanmean(temp_k(indx_date_b{i,j}));
%     elseif sum(size(temp_k(indx_date_b{i,j}))) == 0 
%         chl_bot(i,j) = NaN;
%     end
%     end
% end

t_s_17 = temp_sur; s_s_17 = salt_sur; nh4_s_17 = nh4_sur; no3_s_17 = no3_sur; %chl_s_17 = chl_sur;
t_b_17 = temp_bot; s_b_17 = salt_bot; nh4_b_17 = nh4_bot; no3_b_17 = no3_bot; %chl_b_17 = chl_bot;
clearvars -except '*_s_17' '*_b_17' '*_s_16' '*_b_16' '*_s_15' '*_b_15' '*_s_13' '*_b_13'

% 2013
clearvars eom_d
k=0
for i = 2013:2013
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

for i = 1:12
   eom_13(i) = sum(eom_d(1:i))
end
eom_13 = [1 eom_13(1:end-1)+1]

%2015
clearvars eom_d
k=0
for i = 2015:2015
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

for i = 1:12
   eom_15(i) = sum(eom_d(1:i))
end
eom_15 = [1 eom_15(1:end-1)+1]

%2016
clearvars eom_d
k=0
for i = 2016:2016
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

for i = 1:12
   eom_16(i) = sum(eom_d(1:i))
end
eom_16 = [1 eom_16(1:end-1)+1]

%2017
clearvars eom_d
k=0
for i = 2017:2017
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

for i = 1:12
   eom_17(i) = sum(eom_d(1:i))
end
eom_17 = [1 eom_17(1:end-1)+1]

save('LTME_observation_data.mat')

%% PLOT
close all; clear; clc;
load('LTME_observation_data.mat')
%% nh4
for i = 1:9
fig = figure; hold on;
        plot(nh4_s_13(i,:),'b*','linew',2);
        title(['LTME 2013 daily OBS nh4']);
        xlabel('time(month)','fontsize',13)
        ylabel('nh4 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[eom_13(1:end)]);
        xlim([1 length(nh4_s_13)])
        ylim([0 10])
        set(gca,'xticklabel',1:12,'fontsize',10);
saveas(gcf,[num2str(i) '-st_LTME_2013_daily_OBS_nh4.png']);
close all;
end

for i = 1:9
fig = figure; hold on;
        plot(nh4_s_15(i,:),'b*','linew',2);
        title(['LTME 2015 daily OBS nh4']);
        xlabel('time(month)','fontsize',13)
        ylabel('nh4 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[eom_15(1:end)]);
        xlim([1 length(nh4_s_15)])
        ylim([0 10])
        set(gca,'xticklabel',1:12,'fontsize',10);
saveas(gcf,[num2str(i) '-st_LTME_2015_daily_OBS_nh4.png']);
close all;
end

for i = 1:9
fig = figure; hold on;
        plot(nh4_s_16(i,:),'b*','linew',2);
        title(['LTME 2016 daily OBS nh4']);
        xlabel('time(month)','fontsize',13)
        ylabel('nh4 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[eom_16(1:end)]);
        xlim([1 length(nh4_s_16)])
        ylim([0 10])
        set(gca,'xticklabel',1:12,'fontsize',10);
saveas(gcf,[num2str(i) '-st_LTME_2016_daily_OBS_nh4.png']);
close all;
end

for i = 1:9
fig = figure; hold on;
        plot(nh4_s_17(i,:),'b*','linew',2);
        title(['LTME 2017 daily OBS nh4']);
        xlabel('time(month)','fontsize',13)
        ylabel('nh4 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[eom_17(1:end)]);
        xlim([1 length(nh4_s_17)])
        ylim([0 10])
        set(gca,'xticklabel',1:12,'fontsize',10);
saveas(gcf,[num2str(i) '-st_LTME_2017_daily_OBS_nh4.png']);
close all;
end

%% no3
for i = 1:9
fig = figure; hold on;
        plot(no3_s_13(i,:),'b*','linew',2);
        title(['LTME 2013 daily OBS no3']);
        xlabel('time(month)','fontsize',13)
        ylabel('no3 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[eom_13(1:end)]);
        xlim([1 length(no3_s_13)])
%         ylim([0 100])
        set(gca,'xticklabel',1:12,'fontsize',10);
saveas(gcf,[num2str(i) '-st_LTME_2013_daily_OBS_no3.png']);
% close all;
end

for i = 1:9
fig = figure; hold on;
        plot(no3_s_15(i,:),'b*','linew',2);
        title(['LTME 2015 daily OBS no3']);
        xlabel('time(month)','fontsize',13)
        ylabel('no3 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[eom_15(1:end)]);
        xlim([1 length(no3_s_15)])
%         ylim([0 10])
        set(gca,'xticklabel',1:12,'fontsize',10);
saveas(gcf,[num2str(i) '-st_LTME_2015_daily_OBS_no3.png']);
close all;
end

for i = 1:9
fig = figure; hold on;
        plot(no3_s_16(i,:),'b*','linew',2);
        title(['LTME 2016 daily OBS no3']);
        xlabel('time(month)','fontsize',13)
        ylabel('no3 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[eom_16(1:end)]);
        xlim([1 length(no3_s_16)])
%         ylim([0 10])
        set(gca,'xticklabel',1:12,'fontsize',10);
saveas(gcf,[num2str(i) '-st_LTME_2016_daily_OBS_no3.png']);
close all;
end

for i = 1:9
fig = figure; hold on;
        plot(no3_s_17(i,:),'b*','linew',2);
        title(['LTME 2017 daily OBS no3']);
        xlabel('time(month)','fontsize',13)
        ylabel('no3 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[eom_17(1:end)]);
        xlim([1 length(no3_s_17)])
%         ylim([0 10])
        set(gca,'xticklabel',1:12,'fontsize',10);
saveas(gcf,[num2str(i) '-st_LTME_2017_daily_OBS_no3.png']);
close all;
end

%% chl
for i = 1:9
fig = figure; hold on;
        plot(chl_s_15(i,:),'b*','linew',2);
        title(['LTME 2015 daily OBS chl']);
        xlabel('time(month)','fontsize',13)
        ylabel('chl (mg/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[eom_15(1:end)]);
        xlim([1 length(chl_s_15)])
        ylim([0 10])
        set(gca,'xticklabel',1:12,'fontsize',10);
saveas(gcf,[num2str(i) '-st_LTME_2015_daily_OBS_chl.png']);
close all;
end

for i = 1:9
fig = figure; hold on;
        plot(chl_s_16(i,:),'b*','linew',2);
        title(['LTME 2016 daily OBS chl']);
        xlabel('time(month)','fontsize',13)
        ylabel('chl (mg/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[eom_16(1:end)]);
        xlim([1 length(chl_s_16)])
        ylim([0 10])
        set(gca,'xticklabel',1:12,'fontsize',10);
saveas(gcf,[num2str(i) '-st_LTME_2016_daily_OBS_chl.png']);
close all;
end

%% bot
%% nh4
for i = 1:9
fig = figure; hold on;
        plot(nh4_b_13(i,:),'b*','linew',2);
        title(['LTME 2013 daily OBS nh4']);
        xlabel('time(month)','fontsize',13)
        ylabel('nh4 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[eom_13(1:end)]);
        xlim([1 length(nh4_b_13)])
        ylim([0 10])
        set(gca,'xticklabel',1:12,'fontsize',10);
saveas(gcf,[num2str(i) '-st_LTME_2013_daily_OBS_nh4.png']);
close all;
end

for i = 1:9
fig = figure; hold on;
        plot(nh4_b_15(i,:),'b*','linew',2);
        title(['LTME 2015 daily OBS nh4']);
        xlabel('time(month)','fontsize',13)
        ylabel('nh4 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[eom_15(1:end)]);
        xlim([1 length(nh4_b_15)])
        ylim([0 10])
        set(gca,'xticklabel',1:12,'fontsize',10);
saveas(gcf,[num2str(i) '-st_LTME_2015_daily_OBS_nh4.png']);
close all;
end

for i = 1:9
fig = figure; hold on;
        plot(nh4_b_16(i,:),'b*','linew',2);
        title(['LTME 2016 daily OBS nh4']);
        xlabel('time(month)','fontsize',13)
        ylabel('nh4 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[eom_16(1:end)]);
        xlim([1 length(nh4_b_16)])
        ylim([0 10])
        set(gca,'xticklabel',1:12,'fontsize',10);
saveas(gcf,[num2str(i) '-st_LTME_2016_daily_OBS_nh4.png']);
close all;
end

for i = 1:9
fig = figure; hold on;
        plot(nh4_b_17(i,:),'b*','linew',2);
        title(['LTME 2017 daily OBS nh4']);
        xlabel('time(month)','fontsize',13)
        ylabel('nh4 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[eom_17(1:end)]);
        xlim([1 length(nh4_b_17)])
        ylim([0 10])
        set(gca,'xticklabel',1:12,'fontsize',10);
saveas(gcf,[num2str(i) '-st_LTME_2017_daily_OBS_nh4.png']);
close all;
end

%% no3
for i = 1:9
fig = figure; hold on;
        plot(no3_b_13(i,:),'b*','linew',2);
        title(['LTME 2013 daily OBS no3']);
        xlabel('time(month)','fontsize',13)
        ylabel('no3 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[eom_13(1:end)]);
        xlim([1 length(no3_b_13)])
%         ylim([0 100])
        set(gca,'xticklabel',1:12,'fontsize',10);
saveas(gcf,[num2str(i) '-st_LTME_2013_daily_OBS_no3.png']);
% close all;
end

for i = 1:9
fig = figure; hold on;
        plot(no3_b_15(i,:),'b*','linew',2);
        title(['LTME 2015 daily OBS no3']);
        xlabel('time(month)','fontsize',13)
        ylabel('no3 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[eom_15(1:end)]);
        xlim([1 length(no3_b_15)])
%         ylim([0 10])
        set(gca,'xticklabel',1:12,'fontsize',10);
saveas(gcf,[num2str(i) '-st_LTME_2015_daily_OBS_no3.png']);
close all;
end

for i = 1:9
fig = figure; hold on;
        plot(no3_b_16(i,:),'b*','linew',2);
        title(['LTME 2016 daily OBS no3']);
        xlabel('time(month)','fontsize',13)
        ylabel('no3 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[eom_16(1:end)]);
        xlim([1 length(no3_b_16)])
%         ylim([0 10])
        set(gca,'xticklabel',1:12,'fontsize',10);
saveas(gcf,[num2str(i) '-st_LTME_2016_daily_OBS_no3.png']);
close all;
end

for i = 1:9
fig = figure; hold on;
        plot(no3_b_17(i,:),'b*','linew',2);
        title(['LTME 2017 daily OBS no3']);
        xlabel('time(month)','fontsize',13)
        ylabel('no3 (mmol N / m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[eom_17(1:end)]);
        xlim([1 length(no3_b_17)])
%         ylim([0 10])
        set(gca,'xticklabel',1:12,'fontsize',10);
saveas(gcf,[num2str(i) '-st_LTME_2017_daily_OBS_no3.png']);
close all;
end

%% chl
for i = 1:9
fig = figure; hold on;
        plot(chl_b_15(i,:),'b*','linew',2);
        title(['LTME 2015 daily OBS chl']);
        xlabel('time(month)','fontsize',13)
        ylabel('chl (mg/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[eom_15(1:end)]);
        xlim([1 length(chl_b_15)])
        ylim([0 10])
        set(gca,'xticklabel',1:12,'fontsize',10);
saveas(gcf,[num2str(i) '-st_LTME_2015_daily_OBS_chl.png']);
close all;
end

for i = 1:9
fig = figure; hold on;
        plot(chl_b_16(i,:),'b*','linew',2);
        title(['LTME 2016 daily OBS chl']);
        xlabel('time(month)','fontsize',13)
        ylabel('chl (mg/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        set(gca,'xtick',[eom_16(1:end)]);
        xlim([1 length(chl_b_16)])
        ylim([0 10])
        set(gca,'xticklabel',1:12,'fontsize',10);
saveas(gcf,[num2str(i) '-st_LTME_2016_daily_OBS_chl.png']);
close all;
end

