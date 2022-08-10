close all; clear; clc;   % -v3
% % station name tag
% % (for preventing the miss picking using tag from 1 col's name on excel) 
%coastal
% j=[7;];
% for i = 1:1
% name_tag_7{i} = ['대한해협연안 ' num2str(j(i),'%02d')] 
% end
P_MW = 30.973762;
PO4_MW =94.971482;
N_MW = 14.006720;
NO3_MW = 62.005010;
NH4_MW = 18.038508;


j=[3;];
for i = 1:1
name_tag_7{i} = ['남해' num2str(j(i),'%01d')] 
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

[raw_co txt_co]=xlsread('C:\Users\user\Desktop\GY_메타\해양환경측정망 연안_남해3.xlsx','sheet0','');
txt_matc_co = txt_co(2:end,2);
txt_date_co = txt_co(2:37,5);
for i = 37:length(txt_co)
    txt_date_co{i-1} = {[char(txt_co(i,3)),'-',char(txt_co(i,4))]};
end
temp_sur_co = txt_co(2:end,9); 
temp_bot_co = txt_co(2:end,10); 
salt_sur_co = txt_co(2:end,11); 
salt_bot_co = txt_co(2:end,12);
do_sur_co = txt_co(2:end,15); 
do_bot_co = txt_co(2:end,16);
nh4_sur_co = txt_co(2:end,19); 
nh4_bot_co = txt_co(2:end,20);
no3_sur_co = txt_co(2:end,23); 
no3_bot_co = txt_co(2:end,24);
po4_sur_co = txt_co(2:end,29); 
po4_bot_co = txt_co(2:end,30);
chl_sur_co = txt_co(2:end,37); 
chl_bot_co = txt_co(2:end,38);

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
merge_po4_sur = [po4_sur_co;];
merge_po4_bot = [po4_bot_co;];
merge_chl_sur = [chl_sur_co;];
merge_chl_bot = [chl_bot_co;];

merge_data_txt = [txt_co(2:end,9:end);];


%merge
for i = 1:length(name_tag)
   if  sum(strcmp(name_tag{i}, merge_txt)) ~= 0
       indx{i} = find([strcmp(name_tag{i}, merge_txt)] == 1)     
   end
end

%% make date to be 'yymm' form
for i = 1:length(merge_date)
    clearvars temp
temp = char(merge_date{i});
if size(temp,2) ~= 7
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

%% make 1997 to 2020 'yymm' form
k=0
for i = 1997:2020
    for j = 1:12
        k=k+1;
        ref_date{k,1} = [num2str(i) '-' num2str(j,'%02d')];
    end
end

%% make 1997 to 2020 'yymm' form
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

%po4
clearvars temp
for i = 1:length(indx)
    temp = merge_po4_sur(indx{i});
    for j = 1:size(indx_date,2) %mth
        po4_sur_clim(i,j) = mean(str2num(char(temp(indx_date{i,j}))));
    end
end

clearvars temp
for i = 1:length(indx)
    temp = merge_po4_bot(indx{i});
    for j = 1:size(indx_date,2) %mth
        po4_bot_clim(i,j) = mean(str2num(char(temp(indx_date{i,j}))));
    end
end


%% cut 3 regime
cut_1 = 96 % 2004-12
cut_2 = 228 % 2015-12

chl_sur_1 = chl_sur_clim(1:cut_1);
% chl_sur_2 = chl_sur_clim(cut_1+1:cut_2);
% chl_sur_3 = chl_sur_clim(cut_2+1:end);

chl_bot_1 = chl_bot_clim(1:cut_1);
% chl_bot_2 = chl_bot_clim(cut_1+1:cut_2);
% chl_bot_3 = chl_bot_clim(cut_2+1:end);

po4_sur_1 = po4_sur_clim(1:cut_1);
% po4_sur_2 = po4_sur_clim(cut_1+1:cut_2);
% po4_sur_3 = po4_sur_clim(cut_2+1:end);

po4_bot_1 = po4_bot_clim(1:cut_1);
% po4_bot_2 = po4_bot_clim(cut_1+1:cut_2);
% po4_bot_3 = po4_bot_clim(cut_2+1:end);


no3_sur_1 =  no3_sur_clim(1:cut_1);
% no3_sur_2 =  no3_sur_clim(cut_1+1:cut_2);
% no3_sur_3 =  no3_sur_clim(cut_2+1:end);

no3_bot_1 =  no3_bot_clim(1:cut_1);
% no3_bot_2 =  no3_bot_clim(cut_1+1:cut_2);
% no3_bot_3 =  no3_bot_clim(cut_2+1:end);

nh4_sur_1 =  nh4_sur_clim(1:cut_1);
% nh4_sur_2 =  nh4_sur_clim(cut_1+1:cut_2);
% nh4_sur_3 =  nh4_sur_clim(cut_2+1:end);

nh4_bot_1 =  nh4_bot_clim(1:cut_1);
% nh4_bot_2 =  nh4_bot_clim(cut_1+1:cut_2);
% nh4_bot_3 =  nh4_bot_clim(cut_2+1:end);

do_sur_1 =  do_sur_clim(1:cut_1);
% do_sur_2 =  do_sur_clim(cut_1+1:cut_2);
% do_sur_3 =  do_sur_clim(cut_2+1:end);

do_bot_1 =  do_bot_clim(1:cut_1);
% do_bot_2 =  do_bot_clim(cut_1+1:cut_2);
% do_bot_3 =  do_bot_clim(cut_2+1:end);

salt_sur_1 =  salt_sur_clim(1:cut_1);
% salt_sur_2 =  salt_sur_clim(cut_1+1:cut_2);
% salt_sur_3 =  salt_sur_clim(cut_2+1:end);

salt_bot_1 =  salt_bot_clim(1:cut_1);
% salt_bot_2 =  salt_bot_clim(cut_1+1:cut_2);
% salt_bot_3 =  salt_bot_clim(cut_2+1:end);

temp_sur_1 =  temp_sur_clim(1:cut_1);
% temp_sur_2 =  temp_sur_clim(cut_1+1:cut_2);
% temp_sur_3 =  temp_sur_clim(cut_2+1:end);

temp_bot_1 =  temp_bot_clim(1:cut_1);
% temp_bot_2 =  temp_bot_clim(cut_1+1:cut_2);
% temp_bot_3 =  temp_bot_clim(cut_2+1:end);


%% extract sigma
sig = 3; 

clearvars tempd
tempd=chl_sur_1;
tempd(tempd > nanmean(chl_sur_1) + sig*nanstd(chl_sur_1)) =NaN;
tempd(tempd < nanmean(chl_sur_1) - sig*nanstd(chl_sur_1)) =NaN;
chl_sur_1_s = tempd;



clearvars tempd
tempd=chl_bot_1;
tempd(tempd > nanmean(chl_bot_1) + sig*nanstd(chl_bot_1)) =NaN;
tempd(tempd < nanmean(chl_bot_1) - sig*nanstd(chl_bot_1)) =NaN;
chl_bot_1_s = tempd;


clearvars tempd
tempd=nh4_sur_1;
tempd(tempd > nanmean(nh4_sur_1) + sig*nanstd(nh4_sur_1)) =NaN;
tempd(tempd < nanmean(nh4_sur_1) - sig*nanstd(nh4_sur_1)) =NaN;
nh4_sur_1_s = tempd;



clearvars tempd
tempd=nh4_bot_1;
tempd(tempd > nanmean(nh4_bot_1) + sig*nanstd(nh4_bot_1)) =NaN;
tempd(tempd < nanmean(nh4_bot_1) - sig*nanstd(nh4_bot_1)) =NaN;
nh4_bot_1_s = tempd;


clearvars tempd
tempd=no3_sur_1;
tempd(tempd > nanmean(no3_sur_1) + sig*nanstd(no3_sur_1)) =NaN;
tempd(tempd < nanmean(no3_sur_1) - sig*nanstd(no3_sur_1)) =NaN;
no3_sur_1_s = tempd;


clearvars tempd
tempd=no3_bot_1;
tempd(tempd > nanmean(no3_bot_1) + sig*nanstd(no3_bot_1)) =NaN;
tempd(tempd < nanmean(no3_bot_1) - sig*nanstd(no3_bot_1)) =NaN;
no3_bot_1_s = tempd;


clearvars tempd
tempd=po4_sur_1;
tempd(tempd > nanmean(po4_sur_1) + sig*nanstd(po4_sur_1)) =NaN;
tempd(tempd < nanmean(po4_sur_1) - sig*nanstd(po4_sur_1)) =NaN;
po4_sur_1_s = tempd;


clearvars tempd
tempd=po4_bot_1;
tempd(tempd > nanmean(po4_bot_1) + sig*nanstd(po4_bot_1)) =NaN;
tempd(tempd < nanmean(po4_bot_1) - sig*nanstd(po4_bot_1)) =NaN;
po4_bot_1_s = tempd;

clearvars tempd
tempd=do_sur_1;
tempd(tempd > nanmean(do_sur_1) + sig*nanstd(do_sur_1)) =NaN;
tempd(tempd < nanmean(do_sur_1) - sig*nanstd(do_sur_1)) =NaN;
do_sur_1_s = tempd;


clearvars tempd
tempd=do_bot_1;
tempd(tempd > nanmean(do_bot_1) + sig*nanstd(do_bot_1)) =NaN;
tempd(tempd < nanmean(do_bot_1) - sig*nanstd(do_bot_1)) =NaN;
do_bot_1_s = tempd;


clearvars tempd
tempd=temp_sur_1;
tempd(tempd > nanmean(temp_sur_1) + sig*nanstd(temp_sur_1)) =NaN;
tempd(tempd < nanmean(temp_sur_1) - sig*nanstd(temp_sur_1)) =NaN;
temp_sur_1_s = tempd;


clearvars tempd
tempd=temp_bot_1;
tempd(tempd > nanmean(temp_bot_1) + sig*nanstd(temp_bot_1)) =NaN;
tempd(tempd < nanmean(temp_bot_1) - sig*nanstd(temp_bot_1)) =NaN;
temp_bot_1_s = tempd;


clearvars tempd
tempd=salt_sur_1;
tempd(tempd > nanmean(salt_sur_1) + sig*nanstd(salt_sur_1)) =NaN;
tempd(tempd < nanmean(salt_sur_1) - sig*nanstd(salt_sur_1)) =NaN;
salt_sur_1_s = tempd;


clearvars tempd
tempd=salt_bot_1;
tempd(tempd > nanmean(salt_bot_1) + sig*nanstd(salt_bot_1)) =NaN;
tempd(tempd < nanmean(salt_bot_1) - sig*nanstd(salt_bot_1)) =NaN;
salt_bot_1_s = tempd;


% make monthly climate
for i = 1:12 %month   
mon_clim_sur_chl_1(i) = nanmean(chl_sur_1_s(i:12:end));

mon_clim_bot_chl_1(i) = nanmean(chl_bot_1_s(i:12:end));

mon_clim_sur_nh4_1(i) = nanmean(nh4_sur_1_s(i:12:end));
 
mon_clim_bot_nh4_1(i) = nanmean(nh4_bot_1_s(i:12:end));

mon_clim_sur_no3_1(i) = nanmean(no3_sur_1_s(i:12:end));

mon_clim_bot_no3_1(i) = nanmean(no3_bot_1_s(i:12:end));

mon_clim_sur_po4_1(i) = nanmean(po4_sur_1_s(i:12:end));

mon_clim_bot_po4_1(i) = nanmean(po4_bot_1_s(i:12:end));

mon_clim_sur_do_1(i) = nanmean(do_sur_1_s(i:12:end));

mon_clim_bot_do_1(i) = nanmean(do_bot_1_s(i:12:end));

mon_clim_sur_temp_1(i) = nanmean(temp_sur_1_s(i:12:end));

mon_clim_bot_temp_1(i) = nanmean(temp_bot_1_s(i:12:end));

mon_clim_sur_salt_1(i) = nanmean(salt_sur_1_s(i:12:end));

mon_clim_bot_salt_1(i) = nanmean(salt_bot_1_s(i:12:end));

end
%%
mon_clim_sur_chl_1 = [mon_clim_sur_chl_1, mon_clim_sur_chl_1, mon_clim_sur_chl_1]

mon_clim_bot_chl_1 = [mon_clim_bot_chl_1, mon_clim_bot_chl_1, mon_clim_bot_chl_1]

mon_clim_sur_no3_1 = [mon_clim_sur_no3_1, mon_clim_sur_no3_1, mon_clim_sur_no3_1]

mon_clim_bot_no3_1 = [mon_clim_bot_no3_1, mon_clim_bot_no3_1, mon_clim_bot_no3_1]

mon_clim_sur_po4_1 = [mon_clim_sur_po4_1, mon_clim_sur_po4_1, mon_clim_sur_po4_1]

mon_clim_bot_po4_1 = [mon_clim_bot_po4_1, mon_clim_bot_po4_1, mon_clim_bot_po4_1]

mon_clim_sur_do_1 = [mon_clim_sur_do_1, mon_clim_sur_do_1, mon_clim_sur_do_1]

mon_clim_bot_do_1 = [mon_clim_bot_do_1, mon_clim_bot_do_1, mon_clim_bot_do_1]

mon_clim_sur_nh4_1 = [mon_clim_sur_nh4_1, mon_clim_sur_nh4_1, mon_clim_sur_nh4_1]

mon_clim_bot_nh4_1 = [mon_clim_bot_nh4_1, mon_clim_bot_nh4_1, mon_clim_bot_nh4_1]

mon_clim_sur_temp_1 = [mon_clim_sur_temp_1, mon_clim_sur_temp_1, mon_clim_sur_temp_1]

mon_clim_bot_temp_1 = [mon_clim_bot_temp_1, mon_clim_bot_temp_1, mon_clim_bot_temp_1]

mon_clim_sur_salt_1 = [mon_clim_sur_salt_1, mon_clim_sur_salt_1, mon_clim_sur_salt_1]

mon_clim_bot_salt_1 = [mon_clim_bot_salt_1, mon_clim_bot_salt_1, mon_clim_bot_salt_1]

%%

t=1:length(mon_clim_sur_chl_1);
mon_clim_sur_chl_1(isnan(mon_clim_sur_chl_1)) = interp1(t(~isnan(mon_clim_sur_chl_1)),mon_clim_sur_chl_1(~isnan(mon_clim_sur_chl_1)),t(isnan(mon_clim_sur_chl_1)));

if sum(~isnan(mon_clim_bot_chl_1)) == 0
    mon_clim_bot_chl_1=mon_clim_sur_chl_1;
elseif sum(~isnan(mon_clim_bot_chl_1)) > 0
    t=1:length(mon_clim_bot_chl_1);
    mon_clim_bot_chl_1(isnan(mon_clim_bot_chl_1)) = interp1(t(~isnan(mon_clim_bot_chl_1)),mon_clim_bot_chl_1(~isnan(mon_clim_bot_chl_1)),t(isnan(mon_clim_bot_chl_1)));
end

t=1:length(mon_clim_sur_no3_1);
mon_clim_sur_no3_1(isnan(mon_clim_sur_no3_1)) = interp1(t(~isnan(mon_clim_sur_no3_1)),mon_clim_sur_no3_1(~isnan(mon_clim_sur_no3_1)),t(isnan(mon_clim_sur_no3_1)));

t=1:length(mon_clim_bot_no3_1);
mon_clim_bot_no3_1(isnan(mon_clim_bot_no3_1)) = interp1(t(~isnan(mon_clim_bot_no3_1)),mon_clim_bot_no3_1(~isnan(mon_clim_bot_no3_1)),t(isnan(mon_clim_bot_no3_1)));

t=1:length(mon_clim_sur_po4_1);
mon_clim_sur_po4_1(isnan(mon_clim_sur_po4_1)) = interp1(t(~isnan(mon_clim_sur_po4_1)),mon_clim_sur_po4_1(~isnan(mon_clim_sur_po4_1)),t(isnan(mon_clim_sur_po4_1)));

t=1:length(mon_clim_bot_po4_1);
mon_clim_bot_po4_1(isnan(mon_clim_bot_po4_1)) = interp1(t(~isnan(mon_clim_bot_po4_1)),mon_clim_bot_po4_1(~isnan(mon_clim_bot_po4_1)),t(isnan(mon_clim_bot_po4_1)));

t=1:length(mon_clim_sur_nh4_1);
mon_clim_sur_nh4_1(isnan(mon_clim_sur_nh4_1)) = interp1(t(~isnan(mon_clim_sur_nh4_1)),mon_clim_sur_nh4_1(~isnan(mon_clim_sur_nh4_1)),t(isnan(mon_clim_sur_nh4_1)));

t=1:length(mon_clim_bot_nh4_1);
mon_clim_bot_nh4_1(isnan(mon_clim_bot_nh4_1)) = interp1(t(~isnan(mon_clim_bot_nh4_1)),mon_clim_bot_nh4_1(~isnan(mon_clim_bot_nh4_1)),t(isnan(mon_clim_bot_nh4_1)));

t=1:length(mon_clim_sur_do_1);
mon_clim_sur_do_1(isnan(mon_clim_sur_do_1)) = interp1(t(~isnan(mon_clim_sur_do_1)),mon_clim_sur_do_1(~isnan(mon_clim_sur_do_1)),t(isnan(mon_clim_sur_do_1)));

t=1:length(mon_clim_bot_do_1);
mon_clim_bot_do_1(isnan(mon_clim_bot_do_1)) = interp1(t(~isnan(mon_clim_bot_do_1)),mon_clim_bot_do_1(~isnan(mon_clim_bot_do_1)),t(isnan(mon_clim_bot_do_1)));

t=1:length(mon_clim_sur_temp_1);
mon_clim_sur_temp_1(isnan(mon_clim_sur_temp_1)) = interp1(t(~isnan(mon_clim_sur_temp_1)),mon_clim_sur_temp_1(~isnan(mon_clim_sur_temp_1)),t(isnan(mon_clim_sur_temp_1)));

t=1:length(mon_clim_bot_temp_1);
mon_clim_bot_temp_1(isnan(mon_clim_bot_temp_1)) = interp1(t(~isnan(mon_clim_bot_temp_1)),mon_clim_bot_temp_1(~isnan(mon_clim_bot_temp_1)),t(isnan(mon_clim_bot_temp_1)));

t=1:length(mon_clim_sur_salt_1);
mon_clim_sur_salt_1(isnan(mon_clim_sur_salt_1)) = interp1(t(~isnan(mon_clim_sur_salt_1)),mon_clim_sur_salt_1(~isnan(mon_clim_sur_salt_1)),t(isnan(mon_clim_sur_salt_1)));

t=1:length(mon_clim_bot_salt_1);
mon_clim_bot_salt_1(isnan(mon_clim_bot_salt_1)) = interp1(t(~isnan(mon_clim_bot_salt_1)),mon_clim_bot_salt_1(~isnan(mon_clim_bot_salt_1)),t(isnan(mon_clim_bot_salt_1)));

mon_clim_sur_chl_1_in= mon_clim_sur_chl_1(13:24); 
mon_clim_bot_chl_1_in= mon_clim_bot_chl_1(13:24); 

mon_clim_sur_no3_1_in= mon_clim_sur_no3_1(13:24); 
mon_clim_bot_no3_1_in= mon_clim_bot_no3_1(13:24); 

mon_clim_sur_nh4_1_in= mon_clim_sur_nh4_1(13:24); 
mon_clim_bot_nh4_1_in= mon_clim_bot_nh4_1(13:24);

mon_clim_sur_do_1_in= mon_clim_sur_do_1(13:24); 
mon_clim_bot_do_1_in= mon_clim_bot_do_1(13:24); 

mon_clim_sur_temp_1_in= mon_clim_sur_temp_1(13:24); 
mon_clim_bot_temp_1_in= mon_clim_bot_temp_1(13:24); 

mon_clim_sur_salt_1_in= mon_clim_sur_salt_1(13:24); 
mon_clim_bot_salt_1_in= mon_clim_bot_salt_1(13:24); 

save('koem_input_fix_1997to2004_3sig(3regime)_v7.mat');


