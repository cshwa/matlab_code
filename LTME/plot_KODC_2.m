close all; clear; clc;   % -v3

% % station name tag
% % (for preventing the miss picking using tag from 1 col's name on excel) 

% port
name_tag_1{1} = [num2str(0,'%02d'),'번'] 
name_tag_1{2} = [num2str(16,'%02d'),'번'] 
name_tag_1{3} = [num2str(17,'%02d'),'번'] 

% combining the tag and outter point excluding
% name_tag = name_tag_1'; 
% name_tag{end+1:end+length(name_tag_2)} = name_tag_2; 
% name_tag{end+1:end+length(name_tag_3)} = name_tag_3; 
% size_tag = length(name_tag);


%% pick the row on the excel which has same name with tag
[raw_p txt_p]=xlsread('정선해양조사_400line_from80.xls','sheet','');
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
       no3_p(i) = str2num(char(txt_cut(i,14)))/1000;
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
for i = 1:3
clearvars temp temp_ymd temp_ymd_c
temp = char(date_st{i});
temp_ymd=temp(:,1:end-9);

for j = 1:length(temp_ymd)
    temp_ymd_c{j} = temp_ymd(j,:);
end
date_ymd{i,1} = temp_ymd_c;
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
        ref_yymm(m,:)=[num2str(i+1979) '-' num2str(n,'%02d')];
    for j = 1:eom_d(i,n)
        l=l+1; % index for 1~31(or 30 or 29)
        k=k+1; % index for 1~366
        ref_yymmdd{k}=[num2str(i+1979) '-' num2str(n,'%02d') '-'  num2str(j,'%02d')];
    end
    end
end



% matched date 'yymm' form
for j = 1:3 % st. axis
    for i = 1:length(ref_yymmdd) % date axis
       if  sum(strcmp(ref_yymmdd{i}, date_ymd{j})) ~= 0
           indx_date{j,i} = find([strcmp(ref_yymmdd{i}, date_ymd{j})] == 1);     
       end
    end
end

% matched date 'yymm' form on surf & bottom
for j = 1:3 % st. axis
    for i = 1:length(ref_yymmdd) % date axis
       if  sum(strcmp(ref_yymmdd{i}, date_ymd{j})) ~= 0
           clearvars indx_date
           indx_date = find([strcmp(ref_yymmdd{i}, date_ymd{j})] == 1); 
           indx_date_s{j,i} = indx_date(find(depth_kodc{j}(1,indx_date) == 0));
           indx_date_b{j,i} = indx_date(find(depth_kodc{j}(1,indx_date) == max(depth_kodc{j}(1,indx_date))));
       end
    end 
end

return

% clearvars indx_date_*
% for j = 1:3 % st. axis
%     for i = 1:length(ref_yymmdd) % date axis
%        if  sum(strcmp(ref_yymmdd{i}, date_ymd{j})) ~= 0
%            clearvars indx_date
%            indx_date = find([strcmp(ref_yymmdd{i}, date_ymd{j})] == 1);
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
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_p(indx_date_s{i,j}))) ~= 0     
        temp_sur(i,j) = mean(temp_p(indx_date_s{i,j}));
    elseif sum(size(temp_p(indx_date_s{i,j}))) == 0 
        temp_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_s,1)
    for j = 1:size(indx_date_s,2)
    if sum(size(salt_p(indx_date_s{i,j}))) ~= 0     
        salt_sur(i,j) = mean(salt_p(indx_date_s{i,j}));
    elseif sum(size(salt_p(indx_date_s{i,j}))) == 0 
        salt_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_s,1)
    for j = 1:size(indx_date_s,2)
    if sum(size(do_p(indx_date_s{i,j}))) ~= 0     
        do_sur(i,j) = mean(do_p(indx_date_s{i,j}));
    elseif sum(size(do_p(indx_date_s{i,j}))) == 0 
        do_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_s,1)
    for j = 1:size(indx_date_s,2)
    if sum(size(no3_p(indx_date_s{i,j}))) ~= 0     
        no3_sur(i,j) = mean(no3_p(indx_date_s{i,j}));
    elseif sum(size(no3_p(indx_date_s{i,j}))) == 0 
        no3_sur(i,j) = NaN;
    end
    end
end


%bot
for i = 1:size(indx_date_b,1)
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_p(indx_date_b{i,j}))) ~= 0     
        temp_bot(i,j) = mean(temp_p(indx_date_b{i,j}));
    elseif sum(size(temp_p(indx_date_b{i,j}))) == 0 
        temp_bot(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    for j = 1:size(indx_date_b,2)
    if sum(size(salt_p(indx_date_b{i,j}))) ~= 0     
        salt_bot(i,j) = mean(salt_p(indx_date_b{i,j}));
    elseif sum(size(salt_p(indx_date_b{i,j}))) == 0 
        salt_bot(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    for j = 1:size(indx_date_b,2)
    if sum(size(do_p(indx_date_b{i,j}))) ~= 0     
        do_bot(i,j) = mean(do_p(indx_date_b{i,j}));
    elseif sum(size(do_p(indx_date_b{i,j}))) == 0 
        do_bot(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    for j = 1:size(indx_date_b,2)
    if sum(size(no3_p(indx_date_b{i,j}))) ~= 0     
        no3_bot(i,j) = mean(no3_p(indx_date_b{i,j}));
    elseif sum(size(no3_p(indx_date_b{i,j}))) == 0 
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


clearvars temp
for i = 1:length(indx)
    temp = merge_temp_bot(indx{i});
    for j = 1:size(indx_date_mm,2) %mth
        temp_bot_clim(i,j) = mean(str2num(char(temp(indx_date_mm{i,j}))));
    end
end

%salt
clearvars temp
for i = 1:length(indx)
    temp = merge_salt_sur(indx{i});
    for j = 1:size(indx_date_mm,2) %mth
        salt_sur_clim(i,j) = mean(str2num(char(temp(indx_date_mm{i,j}))));
    end
end

clearvars temp
for i = 1:length(indx)
    temp = merge_salt_bot(indx{i});
    for j = 1:size(indx_date_mm,2) %mth
        salt_bot_clim(i,j) = mean(str2num(char(temp(indx_date_mm{i,j}))));
    end
end

%do
clearvars temp
for i = 1:length(indx)
    temp = merge_do_sur(indx{i});
    for j = 1:size(indx_date_mm,2) %mth
        do_sur_clim(i,j) = mean(str2num(char(temp(indx_date_mm{i,j}))));
    end
end

clearvars temp
for i = 1:length(indx)
    temp = merge_do_bot(indx{i});
    for j = 1:size(indx_date_mm,2) %mth
        do_bot_clim(i,j) = mean(str2num(char(temp(indx_date_mm{i,j}))));
    end
end

%nh4
clearvars temp
for i = 1:length(indx)
    temp = merge_nh4_sur(indx{i});
    for j = 1:size(indx_date_mm,2) %mth
        nh4_sur_clim(i,j) = mean(str2num(char(temp(indx_date_mm{i,j}))));
    end
end

clearvars temp
for i = 1:length(indx)
    temp = merge_nh4_bot(indx{i});
    for j = 1:size(indx_date_mm,2) %mth
        nh4_bot_clim(i,j) = mean(str2num(char(temp(indx_date_mm{i,j}))));
    end
end

%no3
clearvars temp
for i = 1:length(indx)
    temp = merge_no3_sur(indx{i});
    for j = 1:size(indx_date_mm,2) %mth
        no3_sur_clim(i,j) = mean(str2num(char(temp(indx_date_mm{i,j}))));
    end
end

clearvars temp
for i = 1:length(indx)
    temp = merge_no3_bot(indx{i});
    for j = 1:size(indx_date_mm,2) %mth
        no3_bot_clim(i,j) = mean(str2num(char(temp(indx_date_mm{i,j}))));
    end
end

%chl
clearvars temp
for i = 1:length(indx)
    temp = merge_chl_sur(indx{i});
    for j = 1:size(indx_date_mm,2) %mth
        chl_sur_clim(i,j) = mean(str2num(char(temp(indx_date_mm{i,j}))));
    end
end

clearvars temp
for i = 1:length(indx)
    temp = merge_chl_bot(indx{i});
    for j = 1:size(indx_date_mm,2) %mth
        chl_bot_clim(i,j) = mean(str2num(char(temp(indx_date_mm{i,j}))));
    end
end


lon=ncread('grid_sumjin_v1970_fix_3m.nc','lon_rho');
lat=ncread('grid_sumjin_v1970_fix_3m.nc','lat_rho');
mask=ncread('grid_sumjin_v1970_fix_3m.nc','mask_rho');
load KOEM_st_info_(decimal_deg).mat
lat_koem = [lat_1; lat_2; lat_3;];
lon_koem = [lon_1; lon_2; lon_3;];

%% interp

%temp
clearvars temp
for i = 1:size(temp_sur_clim,2)
    temp = griddata(lon_koem,lat_koem,temp_sur_clim(:,i),lon,lat);
    temp_sur_clim_g(:,:,i) = temp;
    clearvars temp
end

clearvars temp
for i = 1:size(temp_bot_clim,2)
    temp = griddata(lon_koem,lat_koem,temp_bot_clim(:,i),lon,lat);
    temp_bot_clim_g(:,:,i) = temp;
    clearvars temp
end

%salt
clearvars temp
for i = 1:size(salt_sur_clim,2)
    temp = griddata(lon_koem,lat_koem,salt_sur_clim(:,i),lon,lat);
    salt_sur_clim_g(:,:,i) = temp;
    clearvars temp
end

clearvars temp
for i = 1:size(salt_bot_clim,2)
    temp = griddata(lon_koem,lat_koem,salt_bot_clim(:,i),lon,lat);
    salt_bot_clim_g(:,:,i) = temp;
    clearvars temp
end

%do
clearvars temp
for i = 1:size(do_sur_clim,2)
    temp = griddata(lon_koem,lat_koem,do_sur_clim(:,i),lon,lat);
    do_sur_clim_g(:,:,i) = temp;
    clearvars temp
end

clearvars temp
for i = 1:size(do_bot_clim,2)
    temp = griddata(lon_koem,lat_koem,do_bot_clim(:,i),lon,lat);
    do_bot_clim_g(:,:,i) = temp;
    clearvars temp
end

%nh4
clearvars temp
for i = 1:size(nh4_sur_clim,2)
    temp = griddata(lon_koem,lat_koem,nh4_sur_clim(:,i),lon,lat);
    nh4_sur_clim_g(:,:,i) = temp;
    clearvars temp
end

clearvars temp
for i = 1:size(nh4_bot_clim,2)
    temp = griddata(lon_koem,lat_koem,nh4_bot_clim(:,i),lon,lat);
    nh4_bot_clim_g(:,:,i) = temp;
    clearvars temp
end

%no3
clearvars temp
for i = 1:size(no3_sur_clim,2)
    temp = griddata(lon_koem,lat_koem,no3_sur_clim(:,i),lon,lat);
    no3_sur_clim_g(:,:,i) = temp;
    clearvars temp
end

clearvars temp
for i = 1:size(no3_bot_clim,2)
    temp = griddata(lon_koem,lat_koem,no3_bot_clim(:,i),lon,lat);
    no3_bot_clim_g(:,:,i) = temp;
    clearvars temp
end

%chl
clearvars temp
for i = 1:size(chl_sur_clim,2)
    temp = griddata(lon_koem,lat_koem,chl_sur_clim(:,i),lon,lat);
    chl_sur_clim_g(:,:,i) = temp;
    clearvars temp
end

clearvars temp
for i = 1:size(chl_bot_clim,2)
    temp = griddata(lon_koem,lat_koem,chl_bot_clim(:,i),lon,lat);
    chl_bot_clim_g(:,:,i) = temp;
    clearvars temp
end


%% check
% for i = 1: 11; figure; pcolor(lon,lat,temp_sur_clim_g(:,:,i).*(mask./mask)); shading interp; end
% for i = 1: 11; figure; pcolor(lon,lat,salt_sur_clim_g(:,:,i).*(mask./mask)); shading interp; end
% for i = 1: 11; figure; pcolor(lon,lat,do_sur_clim_g(:,:,i).*(mask./mask)); shading interp; end
% for i = 1: 11; figure; pcolor(lon,lat,nh4_sur_clim_g(:,:,i).*(mask./mask)); shading interp; end
% for i = 1: 11; figure; pcolor(lon,lat,no3_sur_clim_g(:,:,i).*(mask./mask)); shading interp; end
% for i = 1: 11; figure; pcolor(lon,lat,chl_sur_clim_g(:,:,i).*(mask./mask)); shading interp; end
% 
% for i = 1: 11; figure; pcolor(lon,lat,temp_bot_clim_g(:,:,i).*(mask./mask)); shading interp; end
% for i = 1: 11; figure; pcolor(lon,lat,salt_bot_clim_g(:,:,i).*(mask./mask)); shading interp; end
% for i = 1: 11; figure; pcolor(lon,lat,do_bot_clim_g(:,:,i).*(mask./mask)); shading interp; end
% for i = 1: 11; figure; pcolor(lon,lat,nh4_bot_clim_g(:,:,i).*(mask./mask)); shading interp; end
% for i = 1: 11; figure; pcolor(lon,lat,no3_bot_clim_g(:,:,i).*(mask./mask)); shading interp; end
% for i = 1: 11; figure; pcolor(lon,lat,chl_bot_clim_g(:,:,i).*(mask./mask)); shading interp; end

%% filling nan grid on the ROMS

%temp
clearvars temp
for i = 1:size(temp_sur_clim,2)
    temp=squeeze(temp_sur_clim_g(:,:,i));
    if length(temp(~isnan(temp))) ~= 0  %% there somthing in  the grid
        temp(isnan(temp)) =  griddata(lon(~isnan(temp)),lat(~isnan(temp)),temp(~isnan(temp)),lon(isnan(temp)),lat(isnan(temp)),'nearest');
        temp_sur_clim_g_f(:,:,i) = temp;
    elseif length(temp(~isnan(temp))) == 0 %% there's nothing in the grid
        temp_sur_clim_g_f(:,:,i) = NaN(length(lon(:,1)),length(lon(1,:)));
    end
    clearvars temp
end

clearvars temp
for i = 1:size(temp_bot_clim,2)
    temp=squeeze(temp_bot_clim_g(:,:,i));
    if length(temp(~isnan(temp))) ~= 0  %% there somthing in  the grid
        temp(isnan(temp)) =  griddata(lon(~isnan(temp)),lat(~isnan(temp)),temp(~isnan(temp)),lon(isnan(temp)),lat(isnan(temp)),'nearest');
        temp_bot_clim_g_f(:,:,i) = temp;
    elseif length(temp(~isnan(temp))) == 0 %% there's nothing in the grid
        temp_bot_clim_g_f(:,:,i) = NaN(length(lon(:,1)),length(lon(1,:)));
    end
    clearvars temp
end

%salt
clearvars temp
for i = 1:size(salt_sur_clim,2)
    temp=squeeze(salt_sur_clim_g(:,:,i));
    if length(temp(~isnan(temp))) ~= 0  %% there somthing in  the grid
        temp(isnan(temp)) =  griddata(lon(~isnan(temp)),lat(~isnan(temp)),temp(~isnan(temp)),lon(isnan(temp)),lat(isnan(temp)),'nearest');
        salt_sur_clim_g_f(:,:,i) = temp;
    elseif length(temp(~isnan(temp))) == 0 %% there's nothing in the grid
        salt_sur_clim_g_f(:,:,i) = NaN(length(lon(:,1)),length(lon(1,:)));
    end
    clearvars temp
end

clearvars temp
for i = 1:size(salt_bot_clim,2)
    temp=squeeze(salt_bot_clim_g(:,:,i));
    if length(temp(~isnan(temp))) ~= 0  %% there somthing in  the grid
        temp(isnan(temp)) =  griddata(lon(~isnan(temp)),lat(~isnan(temp)),temp(~isnan(temp)),lon(isnan(temp)),lat(isnan(temp)),'nearest');
        salt_bot_clim_g_f(:,:,i) = temp;
    elseif length(temp(~isnan(temp))) == 0 %% there's nothing in the grid
        salt_bot_clim_g_f(:,:,i) = NaN(length(lon(:,1)),length(lon(1,:)));
    end
    clearvars temp
end


%do
clearvars temp
for i = 1:size(do_sur_clim,2)
    temp=squeeze(do_sur_clim_g(:,:,i));
    if length(temp(~isnan(temp))) ~= 0  %% there somthing in  the grid
        temp(isnan(temp)) =  griddata(lon(~isnan(temp)),lat(~isnan(temp)),temp(~isnan(temp)),lon(isnan(temp)),lat(isnan(temp)),'nearest');
        do_sur_clim_g_f(:,:,i) = temp;
    elseif length(temp(~isnan(temp))) == 0 %% there's nothing in the grid
        do_sur_clim_g_f(:,:,i) = NaN(length(lon(:,1)),length(lon(1,:)));
    end
    clearvars temp
end

clearvars temp
for i = 1:size(do_bot_clim,2)
    temp=squeeze(do_bot_clim_g(:,:,i));
    if length(temp(~isnan(temp))) ~= 0  %% there somthing in  the grid
        temp(isnan(temp)) =  griddata(lon(~isnan(temp)),lat(~isnan(temp)),temp(~isnan(temp)),lon(isnan(temp)),lat(isnan(temp)),'nearest');
        do_bot_clim_g_f(:,:,i) = temp;
    elseif length(temp(~isnan(temp))) == 0 %% there's nothing in the grid
        do_bot_clim_g_f(:,:,i) = NaN(length(lon(:,1)),length(lon(1,:)));
    end
    clearvars temp
end


%nh4
clearvars temp
for i = 1:size(nh4_sur_clim,2)
    temp=squeeze(nh4_sur_clim_g(:,:,i));
    if length(temp(~isnan(temp))) ~= 0  %% there somthing in  the grid
        temp(isnan(temp)) =  griddata(lon(~isnan(temp)),lat(~isnan(temp)),temp(~isnan(temp)),lon(isnan(temp)),lat(isnan(temp)),'nearest');
        nh4_sur_clim_g_f(:,:,i) = temp;
    elseif length(temp(~isnan(temp))) == 0 %% there's nothing in the grid
        nh4_sur_clim_g_f(:,:,i) = NaN(length(lon(:,1)),length(lon(1,:)));
    end
    clearvars temp
end

clearvars temp
for i = 1:size(nh4_bot_clim,2)
    temp=squeeze(nh4_bot_clim_g(:,:,i));
    if length(temp(~isnan(temp))) ~= 0  %% there somthing in  the grid
        temp(isnan(temp)) =  griddata(lon(~isnan(temp)),lat(~isnan(temp)),temp(~isnan(temp)),lon(isnan(temp)),lat(isnan(temp)),'nearest');
        nh4_bot_clim_g_f(:,:,i) = temp;
    elseif length(temp(~isnan(temp))) == 0 %% there's nothing in the grid
        nh4_bot_clim_g_f(:,:,i) = NaN(length(lon(:,1)),length(lon(1,:)));
    end
    clearvars temp
end

%no3
clearvars temp
for i = 1:size(no3_sur_clim,2)
    temp=squeeze(no3_sur_clim_g(:,:,i));
    if length(temp(~isnan(temp))) ~= 0  %% there somthing in  the grid
        temp(isnan(temp)) =  griddata(lon(~isnan(temp)),lat(~isnan(temp)),temp(~isnan(temp)),lon(isnan(temp)),lat(isnan(temp)),'nearest');
        no3_sur_clim_g_f(:,:,i) = temp;
    elseif length(temp(~isnan(temp))) == 0 %% there's nothing in the grid
        no3_sur_clim_g_f(:,:,i) = NaN(length(lon(:,1)),length(lon(1,:)));
    end
    clearvars temp
end

clearvars temp
for i = 1:size(no3_bot_clim,2)
    temp=squeeze(no3_bot_clim_g(:,:,i));
    if length(temp(~isnan(temp))) ~= 0  %% there somthing in  the grid
        temp(isnan(temp)) =  griddata(lon(~isnan(temp)),lat(~isnan(temp)),temp(~isnan(temp)),lon(isnan(temp)),lat(isnan(temp)),'nearest');
        no3_bot_clim_g_f(:,:,i) = temp;
    elseif length(temp(~isnan(temp))) == 0 %% there's nothing in the grid
        no3_bot_clim_g_f(:,:,i) = NaN(length(lon(:,1)),length(lon(1,:)));
    end
    clearvars temp
end

%chl
clearvars temp
for i = 1:size(chl_sur_clim,2)
    temp=squeeze(chl_sur_clim_g(:,:,i));
    if length(temp(~isnan(temp))) ~= 0  %% there somthing in  the grid
        temp(isnan(temp)) =  griddata(lon(~isnan(temp)),lat(~isnan(temp)),temp(~isnan(temp)),lon(isnan(temp)),lat(isnan(temp)),'nearest');
        chl_sur_clim_g_f(:,:,i) = temp;
    elseif length(temp(~isnan(temp))) == 0 %% there's nothing in the grid
        chl_sur_clim_g_f(:,:,i) = NaN(length(lon(:,1)),length(lon(1,:)));
    end
    clearvars temp
end

clearvars temp
for i = 1:size(chl_bot_clim,2)
    temp=squeeze(chl_bot_clim_g(:,:,i));
    if length(temp(~isnan(temp))) ~= 0  %% there somthing in  the grid
        temp(isnan(temp)) =  griddata(lon(~isnan(temp)),lat(~isnan(temp)),temp(~isnan(temp)),lon(isnan(temp)),lat(isnan(temp)),'nearest');
        chl_bot_clim_g_f(:,:,i) = temp;
    elseif length(temp(~isnan(temp))) == 0 %% there's nothing in the grid
        chl_bot_clim_g_f(:,:,i) = NaN(length(lon(:,1)),length(lon(1,:)));
    end
    clearvars temp
end

save('koem_climate.mat', '*_clim_g_f');

clearvars temp
for i = 1:size(lon,1)
    for j = 1:size(lon,2)
            % temp
            temp_t_s=squeeze(temp_sur_clim_g_f(i,j,:)); % pick time series of one point
            temp_t_b=squeeze(temp_bot_clim_g_f(i,j,:)); % pick time series of one point
            
            % salt
            temp_s_s=squeeze(salt_sur_clim_g_f(i,j,:)); % pick time series of one point
            temp_s_b=squeeze(salt_bot_clim_g_f(i,j,:)); % pick time series of one point
            
            % do
            temp_d_s=squeeze(do_sur_clim_g_f(i,j,:)); % pick time series of one point
            temp_d_b=squeeze(do_bot_clim_g_f(i,j,:)); % pick time series of one point
            
            % nh4
            temp_nh_s=squeeze(nh4_sur_clim_g_f(i,j,:)); % pick time series of one point
            temp_nh_b=squeeze(nh4_bot_clim_g_f(i,j,:)); % pick time series of one point
            
            % no3
            temp_no_s=squeeze(no3_sur_clim_g_f(i,j,:)); % pick time series of one point
            temp_no_b=squeeze(no3_bot_clim_g_f(i,j,:)); % pick time series of one point
            
            % chl
            temp_ch_s=squeeze(chl_sur_clim_g_f(i,j,:)); % pick time series of one point
            temp_ch_b=squeeze(chl_bot_clim_g_f(i,j,:)); % pick time series of one point

     
     temp_t_s = [temp_t_s(end); NaN; temp_t_s; NaN; NaN; temp_t_s(2);];
     temp_t_b = [temp_t_b(end); NaN; temp_t_b; NaN; NaN; temp_t_b(2);];
     
     temp_s_s = [temp_s_s(end); NaN; temp_s_s; NaN; NaN; temp_s_s(2);];
     temp_s_b = [temp_s_b(end); NaN; temp_s_b; NaN; NaN; temp_s_b(2);];
     
     temp_d_s = [temp_d_s(end); NaN; temp_d_s; NaN; NaN; temp_d_s(2);];
     temp_d_b = [temp_d_b(end); NaN; temp_d_b; NaN; NaN; temp_d_b(2);];
     
     temp_nh_s = [temp_nh_s(end); NaN; temp_nh_s; NaN; NaN; temp_nh_s(2);];
     temp_nh_b = [temp_nh_b(end); NaN; temp_nh_b; NaN; NaN; temp_nh_b(2);];
     
     temp_no_s = [temp_no_s(end); NaN; temp_no_s; NaN; NaN; temp_no_s(2);];
     temp_no_b = [temp_no_b(end); NaN; temp_no_b; NaN; NaN; temp_no_b(2);];
     
     temp_ch_s = [temp_ch_s(end); NaN; temp_ch_s; NaN; NaN; temp_ch_s(2);];
     temp_ch_b = [temp_ch_b(end); NaN; temp_ch_b; NaN; NaN; temp_ch_b(2);];
     
     % filling NaNs on the time order
     t = 1:length(temp_t_s); %general time
     temp_t_s(isnan(temp_t_s)) = interp1(t(~isnan(temp_t_s)), temp_t_s(~isnan(temp_t_s)), t(isnan(temp_t_s)) ); 
     temp_t_b(isnan(temp_t_b)) = interp1(t(~isnan(temp_t_b)), temp_t_b(~isnan(temp_t_b)), t(isnan(temp_t_b)) ); 
     temp_t_s(1:2)=[]; temp_t_s(end-1:end)=[]; % 11:12, 1:2 mth removing
     temp_t_b(1:2)=[]; temp_t_b(end-1:end)=[]; % 11:12, 1:2 mth removing
     
     temp_s_s(isnan(temp_s_s)) = interp1(t(~isnan(temp_s_s)), temp_s_s(~isnan(temp_s_s)), t(isnan(temp_s_s)) ); 
     temp_s_b(isnan(temp_s_b)) = interp1(t(~isnan(temp_s_b)), temp_s_b(~isnan(temp_s_b)), t(isnan(temp_s_b)) );
     temp_s_s(1:2)=[]; temp_s_s(end-1:end)=[]; % 11:12, 1:2 mth removing
     temp_s_b(1:2)=[]; temp_s_b(end-1:end)=[]; % 11:12, 1:2 mth removing
     
     temp_d_s(isnan(temp_d_s)) = interp1(t(~isnan(temp_d_s)), temp_d_s(~isnan(temp_d_s)), t(isnan(temp_d_s)) ); 
     temp_d_b(isnan(temp_d_b)) = interp1(t(~isnan(temp_d_b)), temp_d_b(~isnan(temp_d_b)), t(isnan(temp_d_b)) );
     temp_d_s(1:2)=[]; temp_d_s(end-1:end)=[]; % 11:12, 1:2 mth removing
     temp_d_b(1:2)=[]; temp_d_b(end-1:end)=[]; % 11:12, 1:2 mth removing
     
     temp_nh_s(isnan(temp_nh_s)) = interp1(t(~isnan(temp_nh_s)), temp_nh_s(~isnan(temp_nh_s)), t(isnan(temp_nh_s)) ); 
     temp_nh_b(isnan(temp_nh_b)) = interp1(t(~isnan(temp_nh_b)), temp_nh_b(~isnan(temp_nh_b)), t(isnan(temp_nh_b)) );
     temp_nh_s(1:2)=[]; temp_nh_s(end-1:end)=[]; % 11:12, 1:2 mth removing
     temp_nh_b(1:2)=[]; temp_nh_b(end-1:end)=[]; % 11:12, 1:2 mth removing
     
     temp_no_s(isnan(temp_no_s)) = interp1(t(~isnan(temp_no_s)), temp_no_s(~isnan(temp_no_s)), t(isnan(temp_no_s)) ); 
     temp_no_b(isnan(temp_no_b)) = interp1(t(~isnan(temp_no_b)), temp_no_b(~isnan(temp_no_b)), t(isnan(temp_no_b)) );
     temp_no_s(1:2)=[]; temp_no_s(end-1:end)=[]; % 11:12, 1:2 mth removing
     temp_no_b(1:2)=[]; temp_no_b(end-1:end)=[]; % 11:12, 1:2 mth removing
     
     temp_ch_s(isnan(temp_ch_s)) = interp1(t(~isnan(temp_ch_s)), temp_ch_s(~isnan(temp_ch_s)), t(isnan(temp_ch_s)) ); 
     temp_ch_b(isnan(temp_ch_b)) = interp1(t(~isnan(temp_ch_b)), temp_ch_b(~isnan(temp_ch_b)), t(isnan(temp_ch_b)) );
     temp_ch_s(1:2)=[]; temp_ch_s(end-1:end)=[]; % 11:12, 1:2 mth removing
     temp_ch_b(1:2)=[]; temp_ch_b(end-1:end)=[]; % 11:12, 1:2 mth removing
     
      % interp on the layer order  

     prepretemp = [NaN(20,12)];
     prepretemp(20,1:12) = temp_t_s'; prepretemp(1,1:12) = temp_t_b';
     
     prepresalt = [NaN(20,12)];
     prepresalt(20,1:12) = temp_s_s'; prepresalt(1,1:12) = temp_s_b';
     
     prepredo = [NaN(20,12)];
     prepredo(20,1:12) = temp_d_s'; prepredo(1,1:12) = temp_d_b';
     
     preprenh4 = [NaN(20,12)];
     preprenh4(20,1:12) = temp_nh_s'; preprenh4(1,1:12) = temp_nh_b';
     
     prepreno3 = [NaN(20,12)];
     prepreno3(20,1:12) = temp_no_s'; prepreno3(1,1:12) = temp_no_b';
     
     preprech = [NaN(20,12)];
     preprech(20,1:12) = temp_ch_s'; preprech(1,1:12) = temp_ch_b';
     
          for k = 1:12 %time
              temp_t = squeeze(prepretemp(:,k));
              temp_s = squeeze(prepresalt(:,k));
              temp_do =  squeeze(prepredo(:,k));
              temp_nh4 =  squeeze(preprenh4(:,k));
              temp_no3 =  squeeze(prepreno3(:,k));
              temp_ch =  squeeze(preprech(:,k));
              
              g_t = 1:length(temp_t); %general layer num.
              temp_t(isnan(temp_t))=interp1(g_t(~isnan(temp_t)), temp_t(~isnan(temp_t)), g_t(isnan(temp_t)) );
              temp_s(isnan(temp_s))=interp1(g_t(~isnan(temp_s)), temp_s(~isnan(temp_s)), g_t(isnan(temp_s)) );
              temp_do(isnan(temp_do))=interp1(g_t(~isnan(temp_do)), temp_do(~isnan(temp_do)), g_t(isnan(temp_do)) );
              temp_nh4(isnan(temp_nh4))=interp1(g_t(~isnan(temp_nh4)), temp_nh4(~isnan(temp_nh4)), g_t(isnan(temp_nh4)) );
              temp_no3(isnan(temp_no3))=interp1(g_t(~isnan(temp_no3)), temp_no3(~isnan(temp_no3)), g_t(isnan(temp_no3)) );
              temp_ch(isnan(temp_ch))=interp1(g_t(~isnan(temp_ch)), temp_ch(~isnan(temp_ch)), g_t(isnan(temp_ch)) );
              
              pretemp(:,k) = temp_t;
              presalt(:,k) = temp_s;
              predo(:,k) = temp_do;
              preno4(:,k) = temp_nh4;
              preno3(:,k) = temp_no3;
              prech(:,k) = temp_ch;
              clearvars temp_t temp_s temp_do temp_nh4 temp_no3 temp_ch g_t
          end
               
          % completet
          temperature(i,j,:,:) = pretemp;
          salt(i,j,:,:) = presalt;
          do(i,j,:,:) = predo;
          nh4(i,j,:,:) = preno4;
          no3(i,j,:,:) = preno3;
          chla(i,j,:,:) = prech;
          
          clearvars temp_t* temp_s_* temp_d* temp_n* temp_ch*  pre* 
    end
end

return

save('KOEM_interped_climate.mat'); 

%% figure

%% temp
temp_range=[4:1:28];
for i = 1:12
figure%('position',[550 550 550 550],'PaperUnits','inches','PaperPosition',[0 0 7.5 3.5]);
% set(gca,'Position',[0.2 0.17 0.73 0.65]);
hold on;
% [cs,h]=contour(lon,lat,squeeze(temp_sur_clim_g(:,:,i)),temp_range,'linecolor',[.01 .01 .01],'linewidth',1);
% clabel(cs,h,'fontsize',17,'labelspacing',700,'fontweight','bold');
pcolor(lon,lat,squeeze(temperature(:,:,20,i)).*(mask./mask));
% colormap summer(20);xlim([0 0.3])

H1=colorbar;
set(get(H1,'title'),'string','^oC','fontsize',18,'fontweight','bold');
% set(H1,'position',[.2242 .1915 .25 .04])
caxis([temp_range(1) temp_range(end)]);
xlim([min(min(lon)) max(max(lon))])
ylim([min(min(lat)) max(max(lat))])
shading flat;

ax(1)=gca;
set(gca,'box','on','linewidth',1.5,'fontsize',20)
set(get(ax(1),'xlabel'),'string','Lon','fontsize',20,'fontweight','bold');
set(get(ax(1),'ylabel'),'string','Lat','fontsize',20,'fontweight','bold');
set(gca,'fontweight','bold')

filename=['KOEM_Temp_sur_' num2str(i) 'mth']; 

print('-dpng',filename);
end

temp_range=[4:1:28];
for i = 1:12
figure;
hold on;
pcolor(lon,lat,squeeze(temperature(:,:,1,i)).*(mask./mask));

H1=colorbar;
set(get(H1,'title'),'string','^oC','fontsize',18,'fontweight','bold');
caxis([temp_range(1) temp_range(end)]);
xlim([min(min(lon)) max(max(lon))])
ylim([min(min(lat)) max(max(lat))])
shading flat;

ax(1)=gca;
set(gca,'box','on','linewidth',1.5,'fontsize',20)
set(get(ax(1),'xlabel'),'string','Lon','fontsize',20,'fontweight','bold');
set(get(ax(1),'ylabel'),'string','Lat','fontsize',20,'fontweight','bold');
set(gca,'fontweight','bold')

filename=['KOEM_Temp_bot_' num2str(i) 'mth']; 

print('-dpng',filename);
end

%% salt
salt_range=[2:1:35];
for i = 1:12
figure;
hold on;
pcolor(lon,lat,squeeze(salt(:,:,20,i)).*(mask./mask));

H1=colorbar;
set(get(H1,'title'),'string','PSU','fontsize',18,'fontweight','bold');
% set(H1,'position',[.2242 .1915 .25 .04])
caxis([salt_range(1) salt_range(end)]);
xlim([min(min(lon)) max(max(lon))])
ylim([min(min(lat)) max(max(lat))])
shading flat;

ax(1)=gca;
set(gca,'box','on','linewidth',1.5,'fontsize',20)
set(get(ax(1),'xlabel'),'string','Lon','fontsize',20,'fontweight','bold');
set(get(ax(1),'ylabel'),'string','Lat','fontsize',20,'fontweight','bold');
set(gca,'fontweight','bold')

filename=['KOEM_salt_sur_' num2str(i) 'mth']; 

print('-dpng',filename); 
end

salt_range=[2:1:35];
for i = 1:size(salt_bot_clim,2)
figure;
hold on;
pcolor(lon,lat,squeeze(salt(:,:,1,i)).*(mask./mask));

H1=colorbar;
set(get(H1,'title'),'string','PSU','fontsize',18,'fontweight','bold');
% set(H1,'position',[.2242 .1915 .25 .04])
caxis([salt_range(1) salt_range(end)]);
xlim([min(min(lon)) max(max(lon))])
ylim([min(min(lat)) max(max(lat))])
shading flat;

ax(1)=gca;
set(gca,'box','on','linewidth',1.5,'fontsize',20)
set(get(ax(1),'xlabel'),'string','Lon','fontsize',20,'fontweight','bold');
set(get(ax(1),'ylabel'),'string','Lat','fontsize',20,'fontweight','bold');
set(gca,'fontweight','bold')

filename=['KOEM_salt_bot_' num2str(i) 'mth']; 

print('-dpng',filename); 
end

%% DO
do_range=[4:1:14];
for i = 1:12
figure;
hold on;
pcolor(lon,lat,squeeze(do(:,:,20,i)).*(mask./mask));

H1=colorbar;
set(get(H1,'title'),'string','DO-mg/L','fontsize',18,'fontweight','bold');
% set(H1,'position',[.2242 .1915 .25 .04])
caxis([do_range(1) do_range(end)]);
xlim([min(min(lon)) max(max(lon))])
ylim([min(min(lat)) max(max(lat))])
shading flat;

ax(1)=gca;
set(gca,'box','on','linewidth',1.5,'fontsize',20)
set(get(ax(1),'xlabel'),'string','Lon','fontsize',20,'fontweight','bold');
set(get(ax(1),'ylabel'),'string','Lat','fontsize',20,'fontweight','bold');
set(gca,'fontweight','bold')

filename=['KOEM_do_sur_' num2str(i) 'mth']; 

print('-dpng',filename); 
end

do_range=[4:1:14];
for i = 1:12
figure;
hold on;
pcolor(lon,lat,squeeze(do(:,:,1,i)).*(mask./mask));

H1=colorbar;
set(get(H1,'title'),'string','DO-mg/L','fontsize',18,'fontweight','bold');
% set(H1,'position',[.2242 .1915 .25 .04])
caxis([do_range(1) do_range(end)]);
xlim([min(min(lon)) max(max(lon))])
ylim([min(min(lat)) max(max(lat))])
shading flat;

ax(1)=gca;
set(gca,'box','on','linewidth',1.5,'fontsize',20)
set(get(ax(1),'xlabel'),'string','Lon','fontsize',20,'fontweight','bold');
set(get(ax(1),'ylabel'),'string','Lat','fontsize',20,'fontweight','bold');
set(gca,'fontweight','bold')

filename=['KOEM_do_bot_' num2str(i) 'mth']; 

print('-dpng',filename); 
end

%% NH4
% nh4_range=[0:1:135];
% nh4_range=[0:1:90];
nh4_range=[0:1:205];
for i = 1:12
figure;
hold on;
pcolor(lon,lat,squeeze(nh4(:,:,20,i)).*(mask./mask));

H1=colorbar;
set(get(H1,'title'),'string','nh4-ug/L','fontsize',18,'fontweight','bold');
% set(H1,'position',[.2242 .1915 .25 .04])
caxis([nh4_range(1) nh4_range(end)]);
xlim([min(min(lon)) max(max(lon))])
ylim([min(min(lat)) max(max(lat))])
shading flat;

ax(1)=gca;
set(gca,'box','on','linewidth',1.5,'fontsize',20)
set(get(ax(1),'xlabel'),'string','Lon','fontsize',20,'fontweight','bold');
set(get(ax(1),'ylabel'),'string','Lat','fontsize',20,'fontweight','bold');
set(gca,'fontweight','bold')

filename=['KOEM_nh4_sur_' num2str(i) 'mth']; 

print('-dpng',filename); 
end

nh4_range=[0:1:205];
for i = 1:12
figure;
hold on;
pcolor(lon,lat,squeeze(nh4(:,:,1,i)).*(mask./mask));

H1=colorbar;
set(get(H1,'title'),'string','nh4-ug/L','fontsize',18,'fontweight','bold');
% set(H1,'position',[.2242 .1915 .25 .04])
caxis([nh4_range(1) nh4_range(end)]);
xlim([min(min(lon)) max(max(lon))])
ylim([min(min(lat)) max(max(lat))])
shading flat;

ax(1)=gca;
set(gca,'box','on','linewidth',1.5,'fontsize',20)
set(get(ax(1),'xlabel'),'string','Lon','fontsize',20,'fontweight','bold');
set(get(ax(1),'ylabel'),'string','Lat','fontsize',20,'fontweight','bold');
set(gca,'fontweight','bold')

filename=['KOEM_nh4_bot_' num2str(i) 'mth']; 

print('-dpng',filename); 
end

%% NO3
no3_range=[0:1:650];
for i = 1:12
figure;
hold on;
pcolor(lon,lat,squeeze(no3(:,:,20,i)).*(mask./mask));

H1=colorbar;
set(get(H1,'title'),'string','no3-ug/L','fontsize',18,'fontweight','bold');
% set(H1,'position',[.2242 .1915 .25 .04])
caxis([no3_range(1) no3_range(end)]);
xlim([min(min(lon)) max(max(lon))])
ylim([min(min(lat)) max(max(lat))])
shading flat;

ax(1)=gca;
set(gca,'box','on','linewidth',1.5,'fontsize',20)
set(get(ax(1),'xlabel'),'string','Lon','fontsize',20,'fontweight','bold');
set(get(ax(1),'ylabel'),'string','Lat','fontsize',20,'fontweight','bold');
set(gca,'fontweight','bold')

filename=['KOEM_no3_sur_' num2str(i) 'mth']; 

print('-dpng',filename); 
end

no3_range=[0:1:650];
for i = 1:12
figure;
hold on;
pcolor(lon,lat,squeeze(no3(:,:,1,i)).*(mask./mask));
% colormap summer(20);xlim([0 0.3])

H1=colorbar;
set(get(H1,'title'),'string','no3-ug/L','fontsize',18,'fontweight','bold');
% set(H1,'position',[.2242 .1915 .25 .04])
caxis([no3_range(1) no3_range(end)]);
xlim([min(min(lon)) max(max(lon))])
ylim([min(min(lat)) max(max(lat))])
shading flat;

ax(1)=gca;
set(gca,'box','on','linewidth',1.5,'fontsize',20)
set(get(ax(1),'xlabel'),'string','Lon','fontsize',20,'fontweight','bold');
set(get(ax(1),'ylabel'),'string','Lat','fontsize',20,'fontweight','bold');
set(gca,'fontweight','bold')

filename=['KOEM_no3_bot_' num2str(i) 'mth']; 

print('-dpng',filename); 
end


%% chla
chl_range=[0:1:18];
for i = 1:12
figure;
hold on;
pcolor(lon,lat,squeeze(chla(:,:,20,i)).*(mask./mask));

H1=colorbar;
set(get(H1,'title'),'string','chl-ug/L','fontsize',18,'fontweight','bold');
% set(H1,'position',[.2242 .1915 .25 .04])
caxis([chl_range(1) chl_range(end)]);
xlim([min(min(lon)) max(max(lon))])
ylim([min(min(lat)) max(max(lat))])
shading flat;

ax(1)=gca;
set(gca,'box','on','linewidth',1.5,'fontsize',20)
set(get(ax(1),'xlabel'),'string','Lon','fontsize',20,'fontweight','bold');
set(get(ax(1),'ylabel'),'string','Lat','fontsize',20,'fontweight','bold');
set(gca,'fontweight','bold')

filename=['KOEM_chl_sur_' num2str(i) 'mth']; 

print('-dpng',filename); 
end

chl_range=[0:1:18];
for i = 1:12
figure;
hold on;
pcolor(lon,lat,squeeze(chla(:,:,1,i)).*(mask./mask));

H1=colorbar;
set(get(H1,'title'),'string','chl-ug/L','fontsize',18,'fontweight','bold');
% set(H1,'position',[.2242 .1915 .25 .04])
caxis([chl_range(1) chl_range(end)]);
xlim([min(min(lon)) max(max(lon))])
ylim([min(min(lat)) max(max(lat))])
shading flat;

ax(1)=gca;
set(gca,'box','on','linewidth',1.5,'fontsize',20)
set(get(ax(1),'xlabel'),'string','Lon','fontsize',20,'fontweight','bold');
set(get(ax(1),'ylabel'),'string','Lat','fontsize',20,'fontweight','bold');
set(gca,'fontweight','bold')

filename=['KOEM_chl_bot_' num2str(i) 'mth']; 

print('-dpng',filename); 
end

close all

    % save('KOEM_name_tag.mat','name_tag');