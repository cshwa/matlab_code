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
for i = 1:3
clearvars temp temp_ymd temp_ymd_c
temp = char(date_st{i});
temp_ymd=temp(:,1:end-12);

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
        ref_yymm{m}=[num2str(i+1979) '-' num2str(n,'%02d')];
    for j = 1:eom_d(i,n)
        l=l+1; % index for 1~31(or 30 or 29)
        k=k+1; % index for 1~366
        ref_yymmdd{k}=[num2str(i+1979) '-' num2str(n,'%02d') '-'  num2str(j,'%02d')];
    end
    end
end



% matched date 'yymm' form on surf & bottom
for j = 1:3 % st. axis
    for i = 1:length(ref_yymm) % date axis
       if  sum(strcmp(ref_yymm{i}, date_ymd{j})) ~= 0
           clearvars indx_date
           indx_date = find([strcmp(ref_yymm{i}, date_ymd{j})] == 1); 
           indx_date_s{j,i} = indx_date(find(depth_kodc{j}(1,indx_date) == 0));
           indx_date_b{j,i} = indx_date(find(depth_kodc{j}(1,indx_date) == max(depth_kodc{j}(1,indx_date))));
       end
    end 
end

indx_date_s = cell(3,length(ref_yymm));
indx_date_b = cell(3,length(ref_yymm));
for j = 1:3 % st. axis
    for i = 1:length(ref_yymm) % date axis
       if  sum(strcmp(ref_yymm{i}, date_ymd{j})) ~= 0
           clearvars indx_date
           indx_date = find([strcmp(ref_yymm{i}, date_ymd{j})] == 1); 
           indx_date_s{j,i} = indx_date(find(depth_kodc{j}(1,indx_date) == 0));
           indx_date_b{j,i} = indx_date(find(depth_kodc{j}(1,indx_date) == max(depth_kodc{j}(1,indx_date))));
%        elseif sum(strcmp(ref_yymm{i}, date_ymd{j})) == 0
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
return
save('KODC_data_monthly.mat','*_bot','*_sur','ref_yymm', 'date_ymd');

%% interp >> to do regime shift test

do_sur = do_sur(2,:)';
do_bot = do_bot(2,:)';
no3_sur =  no3_sur(2,:)'.* 14;
no3_bot = no3_bot(2,:)' .* 14;
temp_sur = temp_sur(2,:)';
temp_bot = temp_bot(2,:)';
salt_sur = salt_sur(2,:)';
salt_bot = salt_bot(2,:)';

t=1:length(1980:2019)*12;
%      temp_t_s(isnan(temp_t_s)) = interp1(t(~isnan(temp_t_s)), temp_t_s(~isnan(temp_t_s)), t(isnan(temp_t_s)) ); 

temp_sur(isnan(temp_sur)) = interp1( t(~isnan(temp_sur)), temp_sur(~isnan(temp_sur)), t(isnan(temp_sur)) ); 
temp_bot(isnan(temp_bot)) = interp1( t(~isnan(temp_bot)), temp_bot(~isnan(temp_bot)), t(isnan(temp_bot)) ); 

salt_sur(isnan(salt_sur)) = interp1( t(~isnan(salt_sur)), salt_sur(~isnan(salt_sur)), t(isnan(salt_sur)) ); 
salt_bot(isnan(salt_bot)) = interp1( t(~isnan(salt_bot)), salt_bot(~isnan(salt_bot)), t(isnan(salt_bot)) ); 

do_sur(isnan(do_sur)) = interp1( t(~isnan(do_sur)), do_sur(~isnan(do_sur)), t(isnan(do_sur)) ); 
do_bot(isnan(do_bot)) = interp1( t(~isnan(do_bot)), do_bot(~isnan(do_bot)), t(isnan(do_bot)) ); 

no3_sur(isnan(no3_sur)) = interp1( t(~isnan(no3_sur)), no3_sur(~isnan(no3_sur)), t(isnan(no3_sur)) ); 
no3_bot(isnan(no3_bot)) = interp1( t(~isnan(no3_bot)), no3_bot(~isnan(no3_bot)), t(isnan(no3_bot)) ); 


