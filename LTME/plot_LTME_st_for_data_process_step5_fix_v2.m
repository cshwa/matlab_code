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

save('LTME_observation_data.mat')


close all; clear; clc
cd D:\장기생태\observation\관측\장기생태프로젝트\남해_DB
load LTME_observation_data.mat

% make eom_d
k=0
for i = 2013:2013
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end
t_tick_pre=sum(eom_d,2);

for i = 1:size(eom_d,1)
    for j = 1:size(eom_d,2)
        eom_d_each(i,j) = sum(eom_d(i,1:j));
    end
end


k=0
for i = 2016:2016
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end
t_tick_pre=sum(eom_d,2);

for i = 1:size(eom_d,1)
    for j = 1:size(eom_d,2)
        eom_d_each_lp(i,j) = sum(eom_d(i,1:j));
    end
end


t_s_16(:,59)=[]; t_b_16(:,59)=[]; s_s_16(:,59)=[]; s_b_16(:,59)=[]; 
no3_s_16(:,59)=[]; no3_b_16(:,59)=[]; nh4_s_16(:,59)=[]; nh4_b_16(:,59)=[];
chl_s_16(:,59)=[]; chl_b_16(:,59)=[];


% make monthly mean
clearvars ltme_*
ltme_nh4_b = NaN(9,12,4);
ltme_no3_b = NaN(9,12,4);
ltme_t_b = NaN(9,12,4);
ltme_s_b = NaN(9,12,4);
ltme_chl_b = NaN(9,12,4);
ltme_nh4_s = NaN(9,12,4);
ltme_no3_s = NaN(9,12,4);
ltme_t_s = NaN(9,12,4);
ltme_s_s = NaN(9,12,4);
ltme_chl_s = NaN(9,12,4);
for j = 1:9
    for i =1:12
    if i ==1 
        % 2013
        ltme_nh4_b(j,i,1)=nanmean(nh4_b_13(j,1:eom_d_each(1,i)),2);
        ltme_no3_b(j,i,1)=nanmean(no3_b_13(j,1:eom_d_each(1,i)),2);
        ltme_t_b(j,i,1)=nanmean(t_b_13(j,1:eom_d_each(1,i)),2);
        ltme_s_b(j,i,1)=nanmean(s_b_13(j,1:eom_d_each(1,i)),2);
        ltme_nh4_s(j,i,1)=nanmean(nh4_s_13(j,1:eom_d_each(1,i)),2);
        ltme_no3_s(j,i,1)=nanmean(no3_s_13(j,1:eom_d_each(1,i)),2);
        ltme_t_s(j,i,1)=nanmean(t_s_13(j,1:eom_d_each(1,i)),2);
        ltme_s_s(j,i,1)=nanmean(s_s_13(j,1:eom_d_each(1,i)),2);
        % 2015
        ltme_nh4_b(j,i,2)=nanmean(nh4_b_15(j,1:eom_d_each(1,i)),2);
        ltme_no3_b(j,i,2)=nanmean(no3_b_15(j,1:eom_d_each(1,i)),2);
        ltme_t_b(j,i,2)=nanmean(t_b_15(j,1:eom_d_each(1,i)),2);
        ltme_s_b(j,i,2)=nanmean(s_b_15(j,1:eom_d_each(1,i)),2);
        ltme_chl_b(j,i,2)=nanmean(chl_b_15(j,1:eom_d_each(1,i)),2);
        ltme_nh4_s(j,i,2)=nanmean(nh4_s_15(j,1:eom_d_each(1,i)),2);
        ltme_no3_s(j,i,2)=nanmean(no3_s_15(j,1:eom_d_each(1,i)),2);
        ltme_t_s(j,i,2)=nanmean(t_s_15(j,1:eom_d_each(1,i)),2);
        ltme_s_s(j,i,2)=nanmean(s_s_15(j,1:eom_d_each(1,i)),2);
        ltme_chl_s(j,i,2)=nanmean(chl_s_15(j,1:eom_d_each(1,i)),2);
        % 2016
        ltme_nh4_b(j,i,3)=nanmean(nh4_b_16(j,1:eom_d_each(1,i)),2);
        ltme_no3_b(j,i,3)=nanmean(no3_b_16(j,1:eom_d_each(1,i)),2);
        ltme_t_b(j,i,3)=nanmean(t_b_16(j,1:eom_d_each(1,i)),2);
        ltme_s_b(j,i,3)=nanmean(s_b_16(j,1:eom_d_each(1,i)),2);
        ltme_chl_b(j,i,3)=nanmean(chl_b_16(j,1:eom_d_each(1,i)),2);
        ltme_nh4_s(j,i,3)=nanmean(nh4_s_16(j,1:eom_d_each(1,i)),2);
        ltme_no3_s(j,i,3)=nanmean(no3_s_16(j,1:eom_d_each(1,i)),2);
        ltme_t_s(j,i,3)=nanmean(t_s_16(j,1:eom_d_each(1,i)),2);
        ltme_s_s(j,i,3)=nanmean(s_s_16(j,1:eom_d_each(1,i)),2);
        ltme_chl_s(j,i,3)=nanmean(chl_s_16(j,1:eom_d_each(1,i)),2);
        % 2017
        ltme_nh4_b(j,i,4)=nanmean(nh4_b_17(j,1:eom_d_each(1,i)),2);
        ltme_no3_b(j,i,4)=nanmean(no3_b_17(j,1:eom_d_each(1,i)),2);
        ltme_t_b(j,i,4)=nanmean(t_b_17(j,1:eom_d_each(1,i)),2);
        ltme_s_b(j,i,4)=nanmean(s_b_17(j,1:eom_d_each(1,i)),2);
        ltme_nh4_s(j,i,4)=nanmean(nh4_s_17(j,1:eom_d_each(1,i)),2);
        ltme_no3_s(j,i,4)=nanmean(no3_s_17(j,1:eom_d_each(1,i)),2);
        ltme_t_s(j,i,4)=nanmean(t_s_17(j,1:eom_d_each(1,i)),2);
        ltme_s_s(j,i,4)=nanmean(s_s_17(j,1:eom_d_each(1,i)),2);
    else
        % 2013
        ltme_nh4_b(j,i,1)=nanmean(nh4_b_13(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_no3_b(j,i,1)=nanmean(no3_b_13(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_t_b(j,i,1)=nanmean(t_b_13(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_s_b(j,i,1)=nanmean(s_b_13(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_nh4_s(j,i,1)=nanmean(nh4_s_13(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_no3_s(j,i,1)=nanmean(no3_s_13(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_t_s(j,i,1)=nanmean(t_s_13(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_s_s(j,i,1)=nanmean(s_s_13(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        % 2015
        ltme_nh4_b(j,i,2)=nanmean(nh4_b_15(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_no3_b(j,i,2)=nanmean(no3_b_15(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_t_b(j,i,2)=nanmean(t_b_15(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_s_b(j,i,2)=nanmean(s_b_15(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_chl_b(j,i,2)=nanmean(chl_b_15(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_nh4_s(j,i,2)=nanmean(nh4_s_15(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_no3_s(j,i,2)=nanmean(no3_s_15(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_t_s(j,i,2)=nanmean(t_s_15(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_s_s(j,i,2)=nanmean(s_s_15(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_chl_s(j,i,2)=nanmean(chl_s_15(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        % 2016
        ltme_nh4_b(j,i,3)=nanmean(nh4_b_16(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_no3_b(j,i,3)=nanmean(no3_b_16(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_t_b(j,i,3)=nanmean(t_b_16(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_s_b(j,i,3)=nanmean(s_b_16(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_chl_b(j,i,3)=nanmean(chl_b_16(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_nh4_s(j,i,3)=nanmean(nh4_s_16(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_no3_s(j,i,3)=nanmean(no3_s_16(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_t_s(j,i,3)=nanmean(t_s_16(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_s_s(j,i,3)=nanmean(s_s_16(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_chl_s(j,i,3)=nanmean(chl_s_16(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        % 2017
        ltme_nh4_b(j,i,4)=nanmean(nh4_b_17(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_no3_b(j,i,4)=nanmean(no3_b_17(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_t_b(j,i,4)=nanmean(t_b_17(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_s_b(j,i,4)=nanmean(s_b_17(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_nh4_s(j,i,4)=nanmean(nh4_s_17(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_no3_s(j,i,4)=nanmean(no3_s_17(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_t_s(j,i,4)=nanmean(t_s_17(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
        ltme_s_s(j,i,4)=nanmean(s_s_17(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)),2);
    end
    end
end

ltme_chl_s_clim =nanmean(ltme_chl_s,3); ltme_nh4_b_clim =nanmean(ltme_chl_b,3);
ltme_nh4_s_clim =nanmean(ltme_nh4_s,3); ltme_nh4_b_clim =nanmean(ltme_nh4_b,3);
ltme_no3_s_clim =nanmean(ltme_no3_s,3); ltme_no3_b_clim =nanmean(ltme_no3_b,3);
ltme_t_s_clim =nanmean(ltme_t_s,3); ltme_t_b_clim =nanmean(ltme_t_b,3);
ltme_s_s_clim =nanmean(ltme_s_s,3); ltme_s_b_clim =nanmean(ltme_s_b,3);

clearvars tem_*
tem_chl_s=ltme_chl_s(3:8,:,:); tem_chl_b=ltme_chl_b(3:8,:,:);
tem_nh4_s=ltme_nh4_s(3:8,:,:); tem_nh4_b=ltme_nh4_b(3:8,:,:);
tem_no3_s=ltme_no3_s(3:8,:,:); tem_no3_b=ltme_no3_b(3:8,:,:);
tem_t_s=ltme_t_s(3:8,:,:); tem_t_b=ltme_t_b(3:8,:,:);
tem_s_s=ltme_s_s(3:8,:,:); tem_s_b=ltme_s_b(3:8,:,:);
dim_all_c = size(tem_chl_s,1)*size(tem_chl_s,3)
dim_all = size(tem_t_s,1)*size(tem_t_s,3)

for i = 1:12 %month
tem_chl_s_m(:,i)=reshape(squeeze(tem_chl_s(:,i,:)),dim_all_c,1); tem_chl_b_m(:,i)=reshape(squeeze(tem_chl_b(:,i,:)),dim_all_c,1);
tem_nh4_s_m(:,i)=reshape(squeeze(tem_nh4_s(:,i,:)),dim_all,1); tem_nh4_b_m(:,i)=reshape(squeeze(tem_nh4_b(:,i,:)),dim_all,1);
tem_no3_s_m(:,i)=reshape(squeeze(tem_no3_s(:,i,:)),dim_all,1); tem_no3_b_m(:,i)=reshape(squeeze(tem_no3_b(:,i,:)),dim_all,1);
tem_t_s_m(:,i)=reshape(squeeze(tem_t_s(:,i,:)),dim_all,1); tem_t_b_m(:,i)=reshape(squeeze(tem_t_b(:,i,:)),dim_all,1);
tem_s_s_m(:,i)=reshape(squeeze(tem_s_s(:,i,:)),dim_all,1); tem_s_b_m(:,i)=reshape(squeeze(tem_s_b(:,i,:)),dim_all,1);
end

for i = 1:12 %month
ltme_chl_s_mean(i) = nanmean(tem_chl_s_m(:,i)); ltme_chl_b_mean(i) = nanmean(tem_chl_b_m(:,i));
ltme_nh4_s_mean(i) = nanmean(tem_nh4_s_m(:,i)); ltme_nh4_b_mean(i) = nanmean(tem_nh4_b_m(:,i));
ltme_no3_s_mean(i) = nanmean(tem_no3_s_m(:,i)); ltme_no3_b_mean(i) = nanmean(tem_no3_b_m(:,i));
ltme_t_s_mean(i) = nanmean(tem_t_s_m(:,i)); ltme_t_b_mean(i) = nanmean(tem_t_b_m(:,i));
ltme_s_s_mean(i) = nanmean(tem_s_s_m(:,i)); ltme_s_b_mean(i) = nanmean(tem_s_b_m(:,i));

ltme_chl_s_std(i) = nanstd(tem_chl_s_m(:,i)); ltme_chl_b_std(i) = nanstd(tem_chl_b_m(:,i));
ltme_nh4_s_std(i) = nanstd(tem_nh4_s_m(:,i)); ltme_nh4_b_std(i) = nanstd(tem_nh4_b_m(:,i));
ltme_no3_s_std(i) = nanstd(tem_no3_s_m(:,i)); ltme_no3_b_std(i) = nanstd(tem_no3_b_m(:,i));
ltme_t_s_std(i) = nanstd(tem_t_s_m(:,i)); ltme_t_b_std(i) = nanstd(tem_t_b_m(:,i));
ltme_s_s_std(i) = nanstd(tem_s_s_m(:,i)); ltme_s_b_std(i) = nanstd(tem_s_b_m(:,i));
end

%% make it day form
for i =1:12
if i == 1
    ltme_chl_s_mean_d(1:eom_d_each(1,i)) = ltme_chl_s_mean(i);
    ltme_chl_b_mean_d(1:eom_d_each(1,i)) = ltme_chl_b_mean(i);
    ltme_nh4_s_mean_d(1:eom_d_each(1,i)) = ltme_nh4_s_mean(i);
    ltme_nh4_b_mean_d(1:eom_d_each(1,i)) = ltme_nh4_b_mean(i);
    ltme_no3_s_mean_d(1:eom_d_each(1,i)) = ltme_no3_s_mean(i);
    ltme_no3_b_mean_d(1:eom_d_each(1,i)) = ltme_no3_b_mean(i);
    ltme_t_s_mean_d(1:eom_d_each(1,i)) = ltme_t_s_mean(i);
    ltme_t_b_mean_d(1:eom_d_each(1,i)) = ltme_t_b_mean(i);
    ltme_s_s_mean_d(1:eom_d_each(1,i)) = ltme_s_s_mean(i);
    ltme_s_b_mean_d(1:eom_d_each(1,i)) = ltme_s_b_mean(i);
    ltme_chl_s_std_d(1:eom_d_each(1,i)) = ltme_chl_s_std(i);
    ltme_chl_b_std_d(1:eom_d_each(1,i)) = ltme_chl_b_std(i);
    ltme_nh4_s_std_d(1:eom_d_each(1,i)) = ltme_nh4_s_std(i);
    ltme_nh4_b_std_d(1:eom_d_each(1,i)) = ltme_nh4_b_std(i);
    ltme_no3_s_std_d(1:eom_d_each(1,i)) = ltme_no3_s_std(i);
    ltme_no3_b_std_d(1:eom_d_each(1,i)) = ltme_no3_b_std(i);
    ltme_t_s_std_d(1:eom_d_each(1,i)) = ltme_t_s_std(i);
    ltme_t_b_std_d(1:eom_d_each(1,i)) = ltme_t_b_std(i);
    ltme_s_s_std_d(1:eom_d_each(1,i)) = ltme_s_s_std(i);
    ltme_s_b_std_d(1:eom_d_each(1,i)) = ltme_s_b_std(i);
else
    ltme_chl_s_mean_d(eom_d_each(1,i-1)+1:eom_d_each(1,i)) = ltme_chl_s_mean(i);
    ltme_chl_b_mean_d(eom_d_each(1,i-1)+1:eom_d_each(1,i)) = ltme_chl_b_mean(i);
    ltme_nh4_s_mean_d(eom_d_each(1,i-1)+1:eom_d_each(1,i)) = ltme_nh4_s_mean(i);
    ltme_nh4_b_mean_d(eom_d_each(1,i-1)+1:eom_d_each(1,i)) = ltme_nh4_b_mean(i);
    ltme_no3_s_mean_d(eom_d_each(1,i-1)+1:eom_d_each(1,i)) = ltme_no3_s_mean(i);
    ltme_no3_b_mean_d(eom_d_each(1,i-1)+1:eom_d_each(1,i)) = ltme_no3_b_mean(i);
    ltme_t_s_mean_d(eom_d_each(1,i-1)+1:eom_d_each(1,i)) = ltme_t_s_mean(i);
    ltme_t_b_mean_d(eom_d_each(1,i-1)+1:eom_d_each(1,i)) = ltme_t_b_mean(i);
    ltme_s_s_mean_d(eom_d_each(1,i-1)+1:eom_d_each(1,i)) = ltme_s_s_mean(i);
    ltme_s_b_mean_d(eom_d_each(1,i-1)+1:eom_d_each(1,i)) = ltme_s_b_mean(i);
    ltme_chl_s_std_d(eom_d_each(1,i-1)+1:eom_d_each(1,i)) = ltme_chl_s_std(i);
    ltme_chl_b_std_d(eom_d_each(1,i-1)+1:eom_d_each(1,i)) = ltme_chl_b_std(i);
    ltme_nh4_s_std_d(eom_d_each(1,i-1)+1:eom_d_each(1,i)) = ltme_nh4_s_std(i);
    ltme_nh4_b_std_d(eom_d_each(1,i-1)+1:eom_d_each(1,i)) = ltme_nh4_b_std(i);
    ltme_no3_s_std_d(eom_d_each(1,i-1)+1:eom_d_each(1,i)) = ltme_no3_s_std(i);
    ltme_no3_b_std_d(eom_d_each(1,i-1)+1:eom_d_each(1,i)) = ltme_no3_b_std(i);
    ltme_t_s_std_d(eom_d_each(1,i-1)+1:eom_d_each(1,i)) = ltme_t_s_std(i);
    ltme_t_b_std_d(eom_d_each(1,i-1)+1:eom_d_each(1,i)) = ltme_t_b_std(i);
    ltme_s_s_std_d(eom_d_each(1,i-1)+1:eom_d_each(1,i)) = ltme_s_s_std(i);
    ltme_s_b_std_d(eom_d_each(1,i-1)+1:eom_d_each(1,i)) = ltme_s_b_std(i);
end
end

save('LTME_observation_data_v2.mat','*_mean_d','*_std_d')



% right_here
%% make linear coeff. and save file

%% salt
color_pic = lines(size(regime_salt_yr,1));
marker_sty = {'o','+','x','^','>','h','p','s','d','.','*','v','<','o','+','x','^','>','h','p','s','d','.','*','v','<',...
    'o','+','x','^','>','h','p','s','d','.','*','v','<','o','+','x','^','>','h','p','s','d','.','*','v','<','o','+','x','^','>','h','p','s','d','.','*','v','<'};
xp = 1:22;
j=0
figure; hold on;
for i = 1:size(regime_salt_yr,1)
  clearvars reg_data_salt xp_w_salt pf_w_salt
reg_data_salt = regime_salt_yr(i,:);
if isnan(reg_data_salt(1)) == 0 
    j = j+1;
    xp_w_salt = find(isnan(reg_data_salt)==0);
    pf_w_salt = polyfit(xp_w_salt,reg_data_salt(xp_w_salt),1);
    yp_w_salt(i,:) = polyval(pf_w_salt,xp);
    scatter(1:22,regime_salt_yr(i,:),marker_sty{i},'MarkerEdgeColor',color_pic(i,:));
    plot(1:22, yp_w_salt(i,:),'color',color_pic(i,:));
    coeff_salt(i,:) = pf_w_salt;
    temp_case(j) = i;
end
hold on
end
xlabel('time(year)','fontsize',13)
ylabel('salt (mg/m^3)','fontsize',13)
set(gca,'xtick',[1:2:22]);
set(gca,'xlim',[1 22]);
set(gca,'xticklabel',1997:2:2018);
title('KOEM-표층염분 연평균','fontsize',13)
grid on
set(gca,'fontsize',13)
% ylim([32 35])


%% no3
j=0
figure; hold on;
for i = 1:size(regime_no3_yr,1)
  clearvars reg_data_no3 xp_w_no3 pf_w_no3
reg_data_no3 = regime_no3_yr(i,:);
if isnan(reg_data_no3(1)) == 0 
    j = j+1;
    xp_w_no3 = find(isnan(reg_data_no3)==0);
    pf_w_no3 = polyfit(xp_w_no3,reg_data_no3(xp_w_no3),1);
    yp_w_no3(i,:) = polyval(pf_w_no3,xp);
    scatter(1:22,regime_no3_yr(i,:),marker_sty{i},'MarkerEdgeColor',color_pic(i,:));
    plot(1:22, yp_w_no3(i,:),'color',color_pic(i,:));
    coeff_no3(i,:) = pf_w_no3;
    temp_case(j) = i;
end
hold on
end
xlabel('time(year)','fontsize',13)
ylabel('NO3 (ug/L)','fontsize',13)
set(gca,'xtick',[1:2:22]);
set(gca,'xlim',[1 22]);
set(gca,'xticklabel',1997:2:2018);
title('KOEM관측-표층질산염 연평균','fontsize',13)
grid on
set(gca,'fontsize',13)
% ylim([32 35])
% legend('205-01','205-02','205-03','205-04','205-05')

%temp
j=0
figure; hold on;
for i = 1:size(regime_temp_yr,1)
  clearvars reg_data_temp xp_w_temp pf_w_temp
reg_data_temp = regime_temp_yr(i,:);
if isnan(reg_data_temp(1)) == 0 
    j = j+1;
    xp_w_temp = find(isnan(reg_data_temp)==0);
    pf_w_temp = polyfit(xp_w_temp,reg_data_temp(xp_w_temp),1);
    yp_w_temp(i,:) = polyval(pf_w_temp,xp);
    scatter(1:22,regime_temp_yr(i,:),marker_sty{i},'MarkerEdgeColor',color_pic(i,:));
    plot(1:22, yp_w_temp(i,:),'color',color_pic(i,:));
    coeff_temp(i,:) = pf_w_temp;
    temp_case(j) = i;
end
hold on
end
xlabel('time(year)','fontsize',13)
ylabel('temperature (^oC)','fontsize',13)
set(gca,'xtick',[1:2:22]);
set(gca,'xlim',[1 22]);
set(gca,'xticklabel',1997:2:2018);
title('KOEM관측-표층수온 연평균','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([13 20])

save('kodc_just_mean_yearly_linear_trend.mat','-v7.3');

filename=['KOEM_chl_bot_' num2str(i) 'mth']; 

print('-dpng',filename); 


close all

    % save('KOEM_name_tag.mat','name_tag');