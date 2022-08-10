close all; clear; clc;   % -v3

% input the sigma
% sig = 3; %% sigma
sig = 2; %% sigma

% % station name tag
% % (for preventing the miss picking using tag from 1 col's name on excel) 

% port
for i = 1:3
name_tag_1{i} = ['여수신항 H' num2str(i)] 
end

name_tag_2{1} = ['광양항 H' num2str(1)] 

name_tag_3{1} = ['삼천포항 H' num2str(1)] 

% estuary
for i = 1:5
name_tag_4{i} = ['가막만 ' num2str(i,'%02d')] 
end

for i = 1:25
name_tag_5{i} = ['섬진강하구 ' num2str(i,'%02d')] 
end

for i = 1:2
name_tag_6{i} = ['진주만 ' num2str(i,'%02d')] 
end

%coastal 
for i = 1:28
name_tag_7{i} = ['대한해협연안 ' num2str(i,'%02d')] 
end

% combining the tag and outter point excluding
name_tag = name_tag_1'; 
name_tag{end+1:end+length(name_tag_2)} = name_tag_2; 
name_tag{end+1:end+length(name_tag_3)} = name_tag_3; 
size_tag = length(name_tag);

for i = 1:length(name_tag_4)
name_tag{size_tag+i} = name_tag_4{i}; 
end
size_tag = length(name_tag);

for i = 1:length(name_tag_5)
name_tag{size_tag+i} = name_tag_5{i}; 
end
size_tag = length(name_tag);

for i = 1:length(name_tag_6)
name_tag{size_tag+i} = name_tag_6{i}; 
end
size_tag = length(name_tag);

% %  when skip the outer point
% % for i = 1:length(name_tag_7)-18  % 연안 3, 6, 13:28 has to be remove (18)  [it's located at out of domain] 
% %  new_i = [1:28]; new_i(3)=[]; new_i(6)=[]; new_i(1)=[];
% %  name_tag{size_tag+i} = name_tag_7{new_i(i)};
% % end

for i = 1:length(name_tag_7) % 연안 3, 6, 13:28 has to be remove (18)  [it's located at out of domain] 
 name_tag{size_tag+i} = name_tag_7{i};
end

%% pick the row on the excel which has same name with tag
[raw_p txt_p]=xlsread('해양환경측정망(항만측정망).xls','sheet','');
txt_matc_p = txt_p(3:end,1); % name list
txt_date_p = txt_p(3:end,3); % date list
temp_sur_p = txt_p(3:end,4); 
temp_bot_p = txt_p(3:end,5); 
salt_sur_p = txt_p(3:end,6); 
salt_bot_p = txt_p(3:end,7);
do_sur_p = txt_p(3:end,10); 
do_bot_p = txt_p(3:end,11);
nh4_sur_p = txt_p(3:end,14); 
nh4_bot_p = txt_p(3:end,15);
no3_sur_p = txt_p(3:end,18); 
no3_bot_p = txt_p(3:end,19);
po4_sur_p = txt_p(3:end,24); 
po4_bot_p = txt_p(3:end,25);
chl_sur_p = txt_p(3:end,32); 
chl_bot_p = txt_p(3:end,33);

[raw_es txt_es]=xlsread('해양환경측정망(하천영향및반폐쇄성해역환경측정망).xls','sheet','');
txt_matc_es = txt_es(3:end,1);
txt_date_es = txt_es(3:end,3);
temp_sur_es = txt_es(3:end,4); 
temp_bot_es = txt_es(3:end,5); 
salt_sur_es = txt_es(3:end,6); 
salt_bot_es = txt_es(3:end,7); 
do_sur_es = txt_es(3:end,10); 
do_bot_es = txt_es(3:end,11);
nh4_sur_es = txt_es(3:end,14); 
nh4_bot_es = txt_es(3:end,15);
no3_sur_es = txt_es(3:end,18); 
no3_bot_es = txt_es(3:end,19);
po4_sur_es = txt_es(3:end,24); 
po4_bot_es = txt_es(3:end,25);
chl_sur_es = txt_es(3:end,32); 
chl_bot_es = txt_es(3:end,33);

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
po4_sur_co = txt_co(3:end,24); 
po4_bot_co = txt_co(3:end,25);
chl_sur_co = txt_co(3:end,32); 
chl_bot_co = txt_co(3:end,33);

merge_txt = [txt_matc_p; txt_matc_es; txt_matc_co;]; % name list
merge_date = [txt_date_p; txt_date_es; txt_date_co;]; % date list
merge_temp_sur = [temp_sur_p; temp_sur_es; temp_sur_co;];
merge_temp_bot = [temp_bot_p; temp_bot_es; temp_bot_co;];
merge_salt_sur = [salt_sur_p; salt_sur_es; salt_sur_co;];
merge_salt_bot = [salt_bot_p; salt_bot_es; salt_bot_co;];
merge_po4_sur = [po4_sur_p; po4_sur_es; po4_sur_co;];
merge_po4_bot = [po4_bot_p; po4_bot_es; po4_bot_co;];
merge_do_sur = [do_sur_p; do_sur_es; do_sur_co;];
merge_do_bot = [do_bot_p; do_bot_es; do_bot_co;];
merge_nh4_sur = [nh4_sur_p; nh4_sur_es; nh4_sur_co;];
merge_nh4_bot = [nh4_bot_p; nh4_bot_es; nh4_bot_co;];
merge_no3_sur = [no3_sur_p; no3_sur_es; no3_sur_co;];
merge_no3_bot = [no3_bot_p; no3_bot_es; no3_bot_co;];
merge_chl_sur = [chl_sur_p; chl_sur_es; chl_sur_co;];
merge_chl_bot = [chl_bot_p; chl_bot_es; chl_bot_co;];

merge_data_txt = [txt_p(3:end,4:end); txt_es(3:end,4:end); txt_co(3:end,4:end);];

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


%temp
clearvars temp
for i = 1:length(indx)
    temp = merge_temp_sur(indx{i});
    for j = 1:size(indx_date,2) %mth
        temp_sur_clim(i,j) = mean(str2num(char(temp(indx_date{i,j}))));
    end
end
temp_sur_clim(:,end+1) = NaN; % 2018.12.

clearvars temp
for i = 1:length(indx)
    temp = merge_temp_bot(indx{i});
    for j = 1:size(indx_date,2) %mth
        temp_bot_clim(i,j) = mean(str2num(char(temp(indx_date{i,j}))));
    end
end
temp_bot_clim(:,end+1) = NaN;

%salt
clearvars temp
for i = 1:length(indx)
    temp = merge_salt_sur(indx{i});
    for j = 1:size(indx_date,2) %mth
        salt_sur_clim(i,j) = mean(str2num(char(temp(indx_date{i,j}))));
    end
end
salt_sur_clim(:,end+1) = NaN; % 2018.12.

clearvars temp
for i = 1:length(indx)
    temp = merge_salt_bot(indx{i});
    for j = 1:size(indx_date,2) %mth
        salt_bot_clim(i,j) = mean(str2num(char(temp(indx_date{i,j}))));
    end
end
salt_bot_clim(:,end+1) = NaN; % 2018.12.


%do
clearvars temp
for i = 1:length(indx)
    temp = merge_do_sur(indx{i});
    for j = 1:size(indx_date,2) %mth
        do_sur_clim(i,j) = mean(str2num(char(temp(indx_date{i,j}))));
    end
end
do_sur_clim(:,end+1) = NaN; % 2018.12.

clearvars temp
for i = 1:length(indx)
    temp = merge_do_bot(indx{i});
    for j = 1:size(indx_date,2) %mth
        do_bot_clim(i,j) = mean(str2num(char(temp(indx_date{i,j}))));
    end
end
do_bot_clim(:,end+1) = NaN; % 2018.12.

%nh4
clearvars temp
for i = 1:length(indx)
    temp = merge_nh4_sur(indx{i});
    for j = 1:size(indx_date,2) %mth
        nh4_sur_clim(i,j) = mean(str2num(char(temp(indx_date{i,j}))));
    end
end
nh4_sur_clim(:,end+1) = NaN; % 2018.12.

clearvars temp
for i = 1:length(indx)
    temp = merge_nh4_bot(indx{i});
    for j = 1:size(indx_date,2) %mth
        nh4_bot_clim(i,j) = mean(str2num(char(temp(indx_date{i,j}))));
    end
end
nh4_bot_clim(:,end+1) = NaN; % 2018.12.

%no3
clearvars temp
for i = 1:length(indx)
    temp = merge_no3_sur(indx{i});
    for j = 1:size(indx_date,2) %mth
        no3_sur_clim(i,j) = mean(str2num(char(temp(indx_date{i,j}))));
    end
end
no3_sur_clim(:,end+1) = NaN; % 2018.12.

clearvars temp
for i = 1:length(indx)
    temp = merge_no3_bot(indx{i});
    for j = 1:size(indx_date,2) %mth
        no3_bot_clim(i,j) = mean(str2num(char(temp(indx_date{i,j}))));
    end
end
no3_bot_clim(:,end+1) = NaN; % 2018.12.

%chl
clearvars temp
for i = 1:length(indx)
    temp = merge_chl_sur(indx{i});
    for j = 1:size(indx_date,2) %mth
        chl_sur_clim(i,j) = mean(str2num(char(temp(indx_date{i,j}))));
    end
end
chl_sur_clim(:,end+1) = NaN; % 2018.12.

clearvars temp
for i = 1:length(indx)
    temp = merge_chl_bot(indx{i});
    for j = 1:size(indx_date,2) %mth
        chl_bot_clim(i,j) = mean(str2num(char(temp(indx_date{i,j}))));
    end
end
chl_bot_clim(:,end+1) = NaN; % 2018.12.

%po4
clearvars temp
for i = 1:length(indx)
    temp = merge_po4_sur(indx{i});
    for j = 1:size(indx_date,2) %mth
        po4_sur_clim(i,j) = mean(str2num(char(temp(indx_date{i,j}))));
    end
end
po4_sur_clim(:,end+1) = NaN; % 2018.12.

clearvars temp
for i = 1:length(indx)
    temp = merge_po4_bot(indx{i});
    for j = 1:size(indx_date,2) %mth
        po4_bot_clim(i,j) = mean(str2num(char(temp(indx_date{i,j}))));
    end
end
po4_bot_clim(:,end+1) = NaN; % 2018.12.


lon=ncread('grid_sumjin_v1970_fix_3m.nc','lon_rho');
lat=ncread('grid_sumjin_v1970_fix_3m.nc','lat_rho');
mask=ncread('grid_sumjin_v1970_fix_3m.nc','mask_rho');
load KOEM_st_info_(decimal_deg).mat
lat_koem = [lat_1; lat_2; lat_3;];
lon_koem = [lon_1; lon_2; lon_3;];


% spatial region index
sp_gy = [4,22,28,29,30,32,33,34,35];
for i  = 1:length(name_tag)
    % make it NaN for all st. except gy.
    if sum(i ~= sp_gy) == 9
        do_sur_clim(i,:) = NaN;
        no3_sur_clim(i,:) = NaN;
        temp_sur_clim(i,:) = NaN;
        salt_sur_clim(i,:) = NaN;
        chl_sur_clim(i,:) = NaN;
        nh4_sur_clim(i,:) = NaN;
        
        po4_sur_clim(i,:) = NaN;
        po4_bot_clim(i,:) = NaN;
        
        do_bot_clim(i,:) = NaN;
        no3_bot_clim(i,:) = NaN;
        temp_bot_clim(i,:) = NaN;
        salt_bot_clim(i,:) = NaN;
        chl_bot_clim(i,:) = NaN;
        nh4_bot_clim(i,:) = NaN;
    end
end

%% extract over 3sig
% for i = 1:length(name_tag)
    clearvars regime_*
    
    
    % sur
    clearvars idx1 tempo_data
    tempo_data= do_sur_clim; 
    idx1 = find(isnan(tempo_data) == 0);
    regime_do=tempo_data;
%     regime_do(find(tempo_data > mean(tempo_data(idx1)) + sig*std(tempo_data(idx1))))=NaN;
%     regime_do(find(tempo_data < mean(tempo_data(idx1)) - sig*std(tempo_data(idx1))))=NaN;
    mtx_regime_do = regime_do;
    regm_do = nanmean(regime_do)

    clearvars idx1 tempo_data
    tempo_data= no3_sur_clim; 
    idx1 = find(isnan(tempo_data) == 0);
    regime_no3=tempo_data;
%     regime_no3(find(tempo_data > mean(tempo_data(idx1)) + sig*std(tempo_data(idx1))))=NaN;
%     regime_no3(find(tempo_data < mean(tempo_data(idx1)) - sig*std(tempo_data(idx1))))=NaN;
    mtx_regime_no3 = regime_no3;
    regm_no3 = nanmean(regime_no3)


    clearvars idx1 tempo_data
    tempo_data= temp_sur_clim; 
    idx1 = find(isnan(tempo_data) == 0);
    regime_temp=tempo_data;
%     regime_temp(find(tempo_data > mean(tempo_data(idx1)) + sig*std(tempo_data(idx1))))=NaN;
%     regime_temp(find(tempo_data < mean(tempo_data(idx1)) - sig*std(tempo_data(idx1))))=NaN;
    mtx_regime_temp = regime_temp;
    regm_temp = nanmean(regime_temp)


    clearvars idx1 tempo_data
    tempo_data= salt_sur_clim;
    idx1 = find(isnan(tempo_data) == 0);
    regime_salt=tempo_data;
%     regime_salt(find(tempo_data > mean(tempo_data(idx1)) + sig*std(tempo_data(idx1))))=NaN;
%     regime_salt(find(tempo_data < mean(tempo_data(idx1)) - sig*std(tempo_data(idx1))))=NaN;
    mtx_regime_salt = regime_salt;
    regm_salt = nanmean(regime_salt)

    clearvars idx1 tempo_data
    tempo_data= po4_sur_clim;
    idx1 = find(isnan(tempo_data) == 0);
    regime_po4=tempo_data;
%     regime_po4(find(tempo_data > mean(tempo_data(idx1)) + sig*std(tempo_data(idx1))))=NaN;
%     regime_po4(find(tempo_data < mean(tempo_data(idx1)) - sig*std(tempo_data(idx1))))=NaN;
    mtx_regime_po4 = regime_po4;
    regm_po4 = nanmean(regime_po4)
    
    clearvars idx1 tempo_data
    tempo_data= nh4_sur_clim;
    idx1 = find(isnan(tempo_data) == 0);
    regime_nh4=tempo_data;
%     regime_nh4(find(tempo_data > mean(tempo_data(idx1)) + sig*std(tempo_data(idx1))))=NaN;
%     regime_nh4(find(tempo_data < mean(tempo_data(idx1)) - sig*std(tempo_data(idx1))))=NaN;
    mtx_regime_nh4 = regime_nh4;
    regm_nh4 = nanmean(regime_nh4)
    
    clearvars idx1 tempo_data
    tempo_data= chl_sur_clim;
    idx1 = find(isnan(tempo_data) == 0);
    regime_chl=tempo_data;
%     regime_chl(find(tempo_data > mean(tempo_data(idx1)) + sig*std(tempo_data(idx1))))=NaN;
%     regime_chl(find(tempo_data < mean(tempo_data(idx1)) - sig*std(tempo_data(idx1))))=NaN;
    mtx_regime_chl = regime_chl;
    regm_chl = nanmean(regime_chl)
    
    
    %bot
    clearvars idx1 tempo_data
    tempo_data= do_bot_clim;
    idx1 = find(isnan(tempo_data) == 0);
    regime_do_b=tempo_data;
%     regime_do_b(find(tempo_data > mean(tempo_data(idx1)) + sig*std(tempo_data(idx1))))=NaN;
%     regime_do_b(find(tempo_data < mean(tempo_data(idx1)) - sig*std(tempo_data(idx1))))=NaN;
    mtx_regime_do_b = regime_do_b;
    regm_do_b = nanmean(regime_do_b)

    clearvars idx1 tempo_data
    tempo_data= no3_bot_clim;
    idx1 = find(isnan(tempo_data) == 0);
    regime_no3_b=tempo_data;
%     regime_no3_b(find(tempo_data > mean(tempo_data(idx1)) + sig*std(tempo_data(idx1))))=NaN;
%     regime_no3_b(find(tempo_data < mean(tempo_data(idx1)) - sig*std(tempo_data(idx1))))=NaN;
    mtx_regime_no3_b = regime_no3_b;
    regm_no3_b = nanmean(regime_no3_b)


    clearvars idx1 tempo_data
    tempo_data= temp_bot_clim;
    idx1 = find(isnan(tempo_data) == 0);
    regime_temp_b=tempo_data;
%     regime_temp_b(find(tempo_data > mean(tempo_data(idx1)) + sig*std(tempo_data(idx1))))=NaN;
%     regime_temp_b(find(tempo_data < mean(tempo_data(idx1)) - sig*std(tempo_data(idx1))))=NaN;
    mtx_regime_temp_b = regime_temp_b;
    regm_temp_b = nanmean(regime_temp)


    clearvars idx1 tempo_data
    tempo_data= salt_bot_clim;
    idx1 = find(isnan(tempo_data) == 0);
    regime_salt_b=tempo_data;
%     regime_salt_b(find(tempo_data > mean(tempo_data(idx1)) + sig*std(tempo_data(idx1))))=NaN;
%     regime_salt_b(find(tempo_data < mean(tempo_data(idx1)) - sig*std(tempo_data(idx1))))=NaN;
    mtx_regime_salt_b = regime_salt_b;
    regm_salt_b = nanmean(regime_salt)
  
    clearvars idx1 tempo_data
    tempo_data= po4_bot_clim;
    idx1 = find(isnan(tempo_data) == 0);
    regime_po4_b=tempo_data;
%     regime_po4_b(find(tempo_data > mean(tempo_data(idx1)) + sig*std(tempo_data(idx1))))=NaN;
%     regime_po4_b(find(tempo_data < mean(tempo_data(idx1)) - sig*std(tempo_data(idx1))))=NaN;
    mtx_regime_po4_b = regime_po4_b;
    regm_po4_b = nanmean(regime_po4)
    
    clearvars idx1 tempo_data
    tempo_data= nh4_bot_clim;
    idx1 = find(isnan(tempo_data) == 0);
    regime_nh4_b=tempo_data;
%     regime_nh4_b(find(tempo_data > mean(tempo_data(idx1)) + sig*std(tempo_data(idx1))))=NaN;
%     regime_nh4_b(find(tempo_data < mean(tempo_data(idx1)) - sig*std(tempo_data(idx1))))=NaN;
    mtx_regime_nh4_b = regime_nh4_b;
    regm_nh4_b = nanmean(regime_nh4)

    clearvars idx1 tempo_data
    tempo_data= chl_bot_clim;
    idx1 = find(isnan(tempo_data) == 0);
    regime_chl_b=tempo_data;
%     regime_chl_b(find(tempo_data > mean(tempo_data(idx1)) + sig*std(tempo_data(idx1))))=NaN;
%     regime_chl_b(find(tempo_data < mean(tempo_data(idx1)) - sig*std(tempo_data(idx1))))=NaN;
    mtx_regime_chl_b = regime_chl_b;
    regm_chl_b = nanmean(regime_chl)
    
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
% end

clearvars regime_*
regime_do = mtx_regime_do;
regime_no3 = mtx_regime_no3;
regime_temp = mtx_regime_temp;
regime_salt = mtx_regime_salt;
regime_nh4 = mtx_regime_nh4;
regime_chl = mtx_regime_chl;

regime_do_b = mtx_regime_do_b;
regime_no3_b = mtx_regime_no3_b;
regime_temp_b = mtx_regime_temp_b;
regime_salt_b = mtx_regime_salt_b;
regime_nh4_b = mtx_regime_nh4_b;
regime_chl_b = mtx_regime_chl_b;

regime_po4 = mtx_regime_po4;
regime_po4_b = mtx_regime_po4_b;

save('koem_timeseires_monthly_gy_only_raw.mat');
% save(['koem_timeseires_monthly_gy_only_',num2str(sig),'sig.mat']);
return

for i = 1:length(name_tag)
            clearvars tempo_*
        tempo_do= regime_do(i,:);
        tempo_no3= regime_no3(i,:);
        tempo_temp= regime_temp(i,:);
        tempo_salt= regime_salt(i,:);
        tempo_nh4= regime_nh4(i,:);
        tempo_do_b= regime_do_b(i,:);
        tempo_no3_b= regime_no3_b(i,:);
        tempo_temp_b= regime_temp_b(i,:);
        tempo_salt_b= regime_salt_b(i,:);
        tempo_nh4_b= regime_nh4_b(i,:);
        
    for j = 1:length(1997:2018)
        %make it yearly mean      
        regime_do_yr(i,j) = nanmean(tempo_do((j-1)*12+1:12*j)); 
        regime_no3_yr(i,j) = nanmean(tempo_no3((j-1)*12+1:12*j)); 
        regime_temp_yr(i,j) = nanmean(tempo_temp((j-1)*12+1:12*j)); 
        regime_salt_yr(i,j) = nanmean(tempo_salt((j-1)*12+1:12*j)); 
        regime_nh4_yr(i,j) = nanmean(tempo_nh4((j-1)*12+1:12*j)); 

        regime_do_b_yr(i,j) = nanmean(tempo_do_b((j-1)*12+1:12*j)); 
        regime_no3_b_yr(i,j) = nanmean(tempo_no3_b((j-1)*12+1:12*j)); 
        regime_temp_b_yr(i,j) = nanmean(tempo_temp_b((j-1)*12+1:12*j)); 
        regime_salt_b_yr(i,j) = nanmean(tempo_salt_b((j-1)*12+1:12*j)); 
        regime_nh4_b_yr(i,j) = nanmean(tempo_nh4_b((j-1)*12+1:12*j)); 
    end
end

% return
% = regime_do_yr;
nonan_case_n = find(isnan(regime_no3_yr(:,1))==0);
nonan_case_t = find(isnan(regime_temp_yr(:,1))==0);
nonan_case_s = find(isnan(regime_salt_yr(:,1))==0);
sp_mean_no3_yr = mean(regime_no3_yr(nonan_case_n,:),1);
sp_mean_temp_yr = mean(regime_temp_yr(nonan_case_t,:),1);
sp_mean_salt_yr = mean(regime_salt_yr(nonan_case_s,:),1);
% sp_mean_salt_yr = mean(regime_salt_yr(nonan_case_s,:),1);
% regime_temp_yr(i,j)
% regime_salt_yr(i,j)
% regime_nh4_yr(i,j)

save('koem_timeseires_yearly_gy_only_raw.mat');
% save('koem_timeseires_yearly_gy_only_2sig.mat');
return

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