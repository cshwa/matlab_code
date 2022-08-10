
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  RCP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc; 
cd D:\장기생태\Dynamic\06_river_physics
[rawa1 txta1]=xlsread('여수_AWS_1990to2019.xls','yeosu_1990to1999','');
[rawa2 txta2]=xlsread('여수_AWS_1990to2019.xls','yeosu_2000to2009','');
[rawa3 txta3]=xlsread('여수_AWS_1990to2019.xls','yeosu_2010to2019','');
[rawa4 txta4]=xlsread('여수_AWS_2020.xlsx','여수_AWS_2020','');
rcp_air = load('J:\장기생태_2021\Dynamic\result\CMIP5\tair_test66_2006-2100.mat'); 

cd D:\장기생태\Dynamic\06_river_physics\환경과학원
[raw txt]=xlsread('수질측정지점_악양.xls','검색결과','');
[raw2 txt2]=xlsread('악양_수질_일반측정망_201901-201912.xls','수질(일반측정망)','');
[raw3 txt3]=xlsread('악양_수질_일반측정망_202001-202012.xls','수질(일반측정망)','');

dash_c = '-';

n_2018 =char(txt(2:end,5)); % 2018 date
n_2019 =char(txt2(3:end,3)); % 2018 date
n_2020 =char(txt3(3:end,3)); % 2018 date

n_2018_ymd = [n_2018(:,1:4) repmat(dash_c,length(n_2018(:,1:4)),1) n_2018(:,6:7) repmat(dash_c,length(n_2018(:,1:4)),1) n_2018(:,9:10)];
n_2019_ymd = [n_2019(:,1:4) repmat(dash_c,length(n_2019(:,1:4)),1) n_2019(:,6:7) repmat(dash_c,length(n_2019(:,1:4)),1) n_2019(:,9:10)];
n_2020_ymd = [n_2020(:,1:4) repmat(dash_c,length(n_2020(:,1:4)),1) n_2020(:,6:7) repmat(dash_c,length(n_2020(:,1:4)),1) n_2020(:,9:10)];

r_txt_ud = flipud(txt);
% r_date_txt=[char(r_txt_ud(1:end-2,4)) repmat(dash_c,length(r_txt_ud(1:end-2,4)),1) char(r_txt_ud(1:end-2,5)) ...
%     repmat(dash_c,length(r_txt_ud(1:end-2,4)),1) char(r_txt_ud(1:end-2,6))];
r_date_txt = [flipud(n_2018_ymd(1:end-12,:)); n_2019_ymd; n_2020_ymd;]; % concaternate (merge)

r_temp_txt=[r_txt_ud(13:end-1,13)];
r_temp=str2num(char(r_temp_txt));
r_temp = [r_temp; raw2(:,14); raw3(:,14);]; % concaternate (merge)

temp_air=[rawa1(:,6); rawa2(:,6); rawa3(:,6); rawa4(:,7);];
temp_air_rcp = reshape(rcp_air.var_mer,1,(size(rcp_air.var_mer,1)*size(rcp_air.var_mer,2)));
temp_date=[txta1(2:end,2); txta2(2:end,2); txta3(2:end,2); txta4(2:end,3);];

%% make 2006~2100 year date vector

% make 365 mm-dd
for i = 1:12
  eom_d(i) = eomday(1981,i); % 1981 is no leap-yr
end

clearvars temp_date_noleap_c temp_date_noleap
zi = 0
for i=2006:2100
    for j=1:12
        for k = 1:eom_d(j)
            zi = zi+1;
            temp_date_noleap(zi,:) = [num2str(i),'-',num2str(j,'%02d'),'-',num2str(k,'%02d')];
            temp_date_noleap_c{zi,1} = [num2str(i),'-',num2str(j,'%02d'),'-',num2str(k,'%02d')];
        end
    end
end


%% yearly mean plot
char_temp_date=char(temp_date_noleap);
for  i = 1:length(temp_date)
    temp_date_air_yy{i,1} = char_temp_date(i,1:4); % delete year from air temp date
end

% pick matched date from air temp date
clearvars indx_air_yy
k=0;
for i = 1990:2020
    k=k+1;
       indx_air_yy{k} = find([strcmp(num2str(i), temp_date_air_yy)] == 1)     
end

for i = 1:length(indx_air_yy)
temp_present_air(i)=mean(temp_air(indx_air_yy{i}));
end


figure; hold on;
plot(length(1990:2006):length(1990:2100),mean(rcp_air.var_mer,1));
plot(1:length(1990:2020),temp_present_air);
xticks(1:5:length(1990:2100));
xticklabels(1990:5:2100);
xtickangle(45);
xlim([1 length(1990:2100)]);
grid on; ylabel('temp (^oC)')
legend('RCP85-IPSL-MR','AWS-Yeosu')

%% raw data plot
for i = 1:length(indx_air_yy)
yy_st(i)=indx_air_yy{i}(1)
end

for i = length(indx_air_yy)+1:length(1990:2100)
yy_st(i) = yy_st(i-1) + 365;
end

figure; hold on;
plot(1:length(temp_air),temp_air,'r','linew',2);
plot(indx_air_yy{17}(1):indx_air_yy{17}(1)+(size(rcp_air.var_mer,1) * (size(rcp_air.var_mer,2)))-1,...
    reshape(rcp_air.var_mer,(size(rcp_air.var_mer,1) * (size(rcp_air.var_mer,2))),1),'b');
xticks(yy_st(1:5:end));
xticklabels(1990:5:2100);
xtickangle(45);
xlim([1 yy_st(end)+365]);
grid on; ylabel('temp (^oC)')
legend('RCP85-IPSL-MR','AWS-Yeosu')


for i = 1:length(r_date_txt)
    r_date_txt_c{i,1} = r_date_txt(i,6:end); % delete year
end

clearvars char_temp_date temp_date_air
char_temp_date=char(temp_date_noleap_c); %air_temp date
for  i = 1:length(temp_date_noleap_c)
    temp_date_air{i,1} = char_temp_date(i,6:end); % delete year from air temp date
end

% % make 366 mm-dd
% for i = 1:12
%   eom_d(i) = eomday(1980,i); % 1980 is leap-yr
% end


k=0
for i = 1:12
    l=0
    for j = 1:eom_d(i)
        l=l+1; % index for 1~31(or 30 or 29)
        k=k+1; % index for 1~366
        mth_d_txt(k,:)=[num2str(i,'%02d') '-'  num2str(l,'%02d')]
    end
end

% make it to cell-array
for i = 1:length(mth_d_txt)
    mth_d_txt_c{i,1} = mth_d_txt(i,:); % delete year
end

% pick matched date from water temp date
for i = 1:length(mth_d_txt_c)
       indx{i} = find([strcmp(mth_d_txt_c{i}, r_date_txt_c)] == 1)     
end

% make quasi_climate
for i = 1:length(mth_d_txt_c)
    if size(indx{i},1) == 0
        w_temp(i) = NaN;   
    else
        w_temp(i) = mean(r_temp(indx{i}));
    end
end

% pick matched date from air temp date
for i = 1:length(mth_d_txt_c)
       indx_air{i} = find([strcmp(mth_d_txt_c{i}, temp_date_air)] == 1)     
end

% make climate from air temp
for i = 1:length(mth_d_txt_c)
    if size(indx_air{i},1) == 0
        a_temp(i) = NaN;   
    else
        a_temp(i) = nanmean(temp_air_rcp(indx_air{i}));
    end
end


pf = polyfit(1:length(a_temp),a_temp,7);
xp = 1:length(a_temp);
yp = polyval(pf,xp);

xp_w = find(isnan(w_temp)==0);
pf_w = polyfit(xp_w,w_temp(xp_w),7);
yp_w = polyval(pf_w,xp);

% water raw - mean_structure (climate)
% rr_temp=r_temp;
for i = 1:length(mth_d_txt_c)
    if size(indx{i},1) == 0
     disp('didnt exist the value')  
    else
        rr_temp(indx{i}) =  r_temp(indx{i}) - yp_w(i);
    end
end

% air raw - daily mean climate
temp_air2 = temp_air_rcp;
for i = 1:length(mth_d_txt_c)
    if size(indx_air{i},1) == 0
    disp('didnt exist the value')  
    else
        temp_air2(indx_air{i}) =  temp_air_rcp(indx_air{i}) - yp(i);
    end
end

% make water date into yy mm dd form.
for i = 1:length(r_date_txt)
    r_date_txt_yymmdd{i,1} = r_date_txt(i,:);
end

% find(isnan(temp_air)==1)
% interp_t = 1:length(temp_air);
% temp_air(isnan(temp_air))=interp1(interp_t(~isnan(temp_air)),temp_air(~isnan(temp_air)), interp_t(isnan(temp_air))); % filling middle nan;

% find matched date from yy mm dd form.
j=0
for i = 1:length(r_date_txt_yymmdd)
   if  sum(strcmp(r_date_txt_yymmdd{i}, temp_date_noleap_c)) ~= 0
       j=j+1;
       indx_yymmdd(j) = find([strcmp(r_date_txt_yymmdd{i}, temp_date_noleap_c)] == 1);    
   else
       disp(i)
   end
end

%% scatter plot
clearvars temp_air_re b1 yCalc1 X b YCalc2
temp_air_re = temp_air2(indx_yymmdd);

%regression
figure; plot(temp_air_re); hold on; plot(rr_temp,'r');hold off;
return
%slope y = b1*x
b1 = temp_air_re\rr_temp % x\y for getting slop
yCalc1 = b1*temp_air_re';
% Slope & Intercept y = b0 + b1*x
X = [ones(length(temp_air_re),1) temp_air_re']; %b_0 b_1 
b = X\rr_temp'
yCalc2 = X*b;

% num of samples
clearvars size_in
for i = 1:length(indx)
    size_in(i) = size(indx{i},1);
end

figure;
plot(size_in); hold on;
xlabel('Days','fontsize',13)
ylabel('number of samples','fontsize',13)
title('number of samples (2008~2020)','fontsize',13)
grid on;
xlim([1 366])
legend('number of samples')
set(gca,'fontsize',15)


figure;
scatter(temp_air_re,rr_temp,'.');
hold on
% plot(temp_air_re,yCalc1,'r')
xlabel('air temp. (^oC)','fontsize',13)
ylabel('water temp. (^oC)','fontsize',13)
title('Linear Regression Relation Between air & water temp.')
grid on
plot(temp_air_re,yCalc2,'--','color','r')
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
gtext(['Corr. = ',num2str(corr(temp_air_re',rr_temp'),2)],'Color','k','FontSize',16)
set(gca,'fontsize',13)

figure;
plot(w_temp,'o'); hold on;
plot(a_temp,'r','linew',1.5);
plot(yp,'g','linew',1.5);
plot(yp_w,'c','linew',1.5);
xlabel('Days','fontsize',13)
ylabel('temp. (^oC)','fontsize',13)
title('compare daily climate air & water temp.','fontsize',13)
grid on;
xlim([1 366])
legend('sumjin(agyang)','IPSL-rcp85-air')
set(gca,'fontsize',15)

save('agyang_regression_RCP85_IPSL_MR_extract_climate_2100_high_airT.mat'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc; 
cd D:\장기생태\Dynamic\06_river_physics
[rawa1 txta1]=xlsread('여수_AWS_1990to2019.xls','yeosu_1990to1999','');
[rawa2 txta2]=xlsread('여수_AWS_1990to2019.xls','yeosu_2000to2009','');
[rawa3 txta3]=xlsread('여수_AWS_1990to2019.xls','yeosu_2010to2019','');
[rawa4 txta4]=xlsread('여수_AWS_2020.xlsx','여수_AWS_2020','');

cd D:\장기생태\Dynamic\06_river_physics\환경과학원
[raw txt]=xlsread('수질측정지점_하동.xls','검색결과','');
[raw2 txt2]=xlsread('하동_수질_일반측정망_201901-201912.xls','수질(일반측정망)','');
[raw3 txt3]=xlsread('하동_수질_일반측정망_202001-202012.xls','수질(일반측정망)','');

dash_c = '-';

n_2018 =char(txt(2:end,5)); % 2018 date
n_2019 =char(txt2(3:end,3)); % 2018 date
n_2020 =char(txt3(3:end,3)); % 2018 date

n_2018_ymd = [n_2018(:,1:4) repmat(dash_c,length(n_2018(:,1:4)),1) n_2018(:,6:7) repmat(dash_c,length(n_2018(:,1:4)),1) n_2018(:,9:10)];
n_2019_ymd = [n_2019(:,1:4) repmat(dash_c,length(n_2019(:,1:4)),1) n_2019(:,6:7) repmat(dash_c,length(n_2019(:,1:4)),1) n_2019(:,9:10)];
n_2020_ymd = [n_2020(:,1:4) repmat(dash_c,length(n_2020(:,1:4)),1) n_2020(:,6:7) repmat(dash_c,length(n_2020(:,1:4)),1) n_2020(:,9:10)];

r_txt_ud = flipud(txt);
% r_date_txt=[char(r_txt_ud(1:end-2,4)) repmat(dash_c,length(r_txt_ud(1:end-2,4)),1) char(r_txt_ud(1:end-2,5)) ...
%     repmat(dash_c,length(r_txt_ud(1:end-2,4)),1) char(r_txt_ud(1:end-2,6))];
r_date_txt = [flipud(n_2018_ymd(1:end-12,:)); n_2019_ymd; n_2020_ymd;]; % concaternate (merge)

r_temp_txt=[r_txt_ud(13:end-1,13)];
r_temp=str2num(char(r_temp_txt));
r_temp = [r_temp; raw2(:,14); raw3(:,14);]; % concaternate (merge)

temp_air=[rawa1(:,6); rawa2(:,6); rawa3(:,6); rawa4(:,7);];
temp_date=[txta1(2:end,2); txta2(2:end,2); txta3(2:end,2); txta4(2:end,3);];

for i = 1:length(r_date_txt)
    r_date_txt_c{i,1} = r_date_txt(i,6:end); % delete year
end

char_temp_date=char(temp_date);
for  i = 1:length(temp_date)
    temp_date_air{i,1} = char_temp_date(i,6:end); % delete year from air temp date
end

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
        mth_d_txt(k,:)=[num2str(i,'%02d') '-'  num2str(l,'%02d')]
    end
end

% make it to cell-array
for i = 1:length(mth_d_txt)
    mth_d_txt_c{i,1} = mth_d_txt(i,:); % delete year
end

% pick matched date from water temp date
for i = 1:length(mth_d_txt_c)
       indx{i} = find([strcmp(mth_d_txt_c{i}, r_date_txt_c)] == 1)     
end

% make quasi_climate
for i = 1:length(mth_d_txt_c)
    if size(indx{i},1) == 0
        w_temp(i) = NaN;   
    else
        w_temp(i) = mean(r_temp(indx{i}));
    end
end

% pick matched date from air temp date
for i = 1:length(mth_d_txt_c)
       indx_air{i} = find([strcmp(mth_d_txt_c{i}, temp_date_air)] == 1)     
end

% make climate from air temp
for i = 1:length(mth_d_txt_c)
    if size(indx_air{i},1) == 0
        a_temp(i) = NaN;   
    else
        a_temp(i) = nanmean(temp_air(indx_air{i}));
    end
end


pf = polyfit(1:length(a_temp),a_temp,7);
xp = 1:length(a_temp);
yp = polyval(pf,xp);

xp_w = find(isnan(w_temp)==0);
pf_w = polyfit(xp_w,w_temp(xp_w),7);
yp_w = polyval(pf_w,xp);

% water raw - mean_structure (climate)
% rr_temp=r_temp;
for i = 1:length(mth_d_txt_c)
    if size(indx{i},1) == 0
     disp('didnt exist the value')  
    else
        rr_temp(indx{i}) =  r_temp(indx{i}) - yp_w(i);
    end
end

% air raw - daily mean climate
temp_air2 = temp_air;
for i = 1:length(mth_d_txt_c)
    if size(indx_air{i},1) == 0
    disp('didnt exist the value')  
    else
        temp_air2(indx_air{i}) =  temp_air(indx_air{i}) - yp(i);
    end
end

% make water date into yy mm dd form.
for i = 1:length(r_date_txt)
    r_date_txt_yymmdd{i,1} = r_date_txt(i,:);
end

% find(isnan(temp_air)==1)
% interp_t = 1:length(temp_air);
% temp_air(isnan(temp_air))=interp1(interp_t(~isnan(temp_air)),temp_air(~isnan(temp_air)), interp_t(isnan(temp_air))); % filling middle nan;

% find matched date from yy mm dd form.
j=0
for i = 1:length(r_date_txt_yymmdd)
   if  sum(strcmp(r_date_txt_yymmdd{i}, temp_date)) ~= 0
       j=j+1;
       indx_yymmdd(j) = find([strcmp(r_date_txt_yymmdd{i}, temp_date)] == 1);    
   else
       disp(i)
   end
end

%% scatter plot
temp_air_re = temp_air2(indx_yymmdd);

%regression
figure; plot(temp_air_re); hold on; plot(rr_temp,'r');hold off;
return
%slope y = b1*x
b1 = temp_air_re\rr_temp' % x\y for getting slop
yCalc1 = b1*temp_air_re;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(temp_air_re),1) temp_air_re]; %b_0 b_1 
b = X\rr_temp'
yCalc2 = X*b;

% num of samples
clearvars size_in
for i = 1:length(indx)
    size_in(i) = size(indx{i},1);
end

figure;
plot(size_in); hold on;
xlabel('Days','fontsize',13)
ylabel('number of samples','fontsize',13)
title('number of samples (1994~2019)','fontsize',13)
grid on;
xlim([1 366])
legend('number of samples')
set(gca,'fontsize',15)


figure;
scatter(temp_air_re,rr_temp,'.');
hold on
% plot(temp_air_re,yCalc1,'r')
xlabel('air temp. (^oC)','fontsize',13)
ylabel('water temp. (^oC)','fontsize',13)
title('Linear Regression Relation Between air & water temp.')
grid on
plot(temp_air_re,yCalc2,'--','color','r')
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
gtext(['Corr. = ',num2str(corr(temp_air_re,rr_temp'),2)],'Color','k','FontSize',16)
set(gca,'fontsize',13)

figure;
plot(w_temp,'o'); hold on;
plot(a_temp,'r','linew',1.5);
plot(yp,'g','linew',1.5);
plot(yp_w,'c','linew',1.5);
xlabel('Days','fontsize',13)
ylabel('temp. (^oC)','fontsize',13)
title('compare daily climate air & water temp.','fontsize',13)
grid on;
xlim([1 366])
legend('sumjin(hadong)','yeosu-air')
set(gca,'fontsize',15)

save('hadong_regression_yeosu_extract_climate_2021_high_airT.mat'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc; 
cd D:\장기생태\Dynamic\06_river_physics
[rawa1 txta1]=xlsread('여수_AWS_1990to2019.xls','yeosu_1990to1999','');
[rawa2 txta2]=xlsread('여수_AWS_1990to2019.xls','yeosu_2000to2009','');
[rawa3 txta3]=xlsread('여수_AWS_1990to2019.xls','yeosu_2010to2019','');
[rawa4 txta4]=xlsread('여수_AWS_2020.xlsx','여수_AWS_2020','');

cd D:\장기생태\Dynamic\06_river\환경과학원
[raw txt]=xlsread('수질측정지점_구례_1989-2018.xls','검색결과','');
[raw2 txt2]=xlsread('수질_일반측정망_구례_201901-201912.xls','수질(일반측정망)','');
[raw3 txt3]=xlsread('수질_일반측정망_구례_202001-202012.xls','수질(일반측정망)','');

dash_c = '-';

n_2018 =char(txt(2:end,5)); % 2018 date
n_2019 =char(txt2(3:end,3)); % 2018 date
n_2020 =char(txt3(3:end,3)); % 2018 date

n_2018_ymd = [n_2018(:,1:4) repmat(dash_c,length(n_2018(:,1:4)),1) n_2018(:,6:7) repmat(dash_c,length(n_2018(:,1:4)),1) n_2018(:,9:10)];
n_2019_ymd = [n_2019(:,1:4) repmat(dash_c,length(n_2019(:,1:4)),1) n_2019(:,6:7) repmat(dash_c,length(n_2019(:,1:4)),1) n_2019(:,9:10)];
n_2020_ymd = [n_2020(:,1:4) repmat(dash_c,length(n_2020(:,1:4)),1) n_2020(:,6:7) repmat(dash_c,length(n_2020(:,1:4)),1) n_2020(:,9:10)];

r_txt_ud = flipud(txt);
% r_date_txt=[char(r_txt_ud(1:end-2,4)) repmat(dash_c,length(r_txt_ud(1:end-2,4)),1) char(r_txt_ud(1:end-2,5)) ...
%     repmat(dash_c,length(r_txt_ud(1:end-2,4)),1) char(r_txt_ud(1:end-2,6))];
r_date_txt = [flipud(n_2018_ymd(1:end-12,:)); n_2019_ymd; n_2020_ymd;]; % concaternate (merge)

r_temp_txt=[r_txt_ud(13:end-1,13)];
r_temp=str2num(char(r_temp_txt));
r_temp = [r_temp; raw2(:,14); raw3(:,14);]; % concaternate (merge)

temp_air=[rawa1(:,6); rawa2(:,6); rawa3(:,6); rawa4(:,7);];
temp_date=[txta1(2:end,2); txta2(2:end,2); txta3(2:end,2); txta4(2:end,3);];

for i = 1:length(r_date_txt)
    r_date_txt_c{i,1} = r_date_txt(i,6:end); % delete year
end

char_temp_date=char(temp_date);
for  i = 1:length(temp_date)
    temp_date_air{i,1} = char_temp_date(i,6:end); % delete year from air temp date
end

%% make 366 mm-dd
for i = 1:12
  eom_d(i) = eomday(1980,i); % 1980 is leap-yr
end

k=0
for i = 1:12
    l=0
    for j = 1:eom_d(i)
        l=l+1; % index for 1~31(or 30 or 29)
        k=k+1; % index for 1~366
        mth_d_txt(k,:)=[num2str(i,'%02d') '-'  num2str(l,'%02d')]
    end
end

% make it to cell-array
for i = 1:length(mth_d_txt)
    mth_d_txt_c{i,1} = mth_d_txt(i,:); % delete year
end

% pick matched date from water temp date
for i = 1:length(mth_d_txt_c)
       indx{i} = find([strcmp(mth_d_txt_c{i}, r_date_txt_c)] == 1)     
end

% make quasi_climate
for i = 1:length(mth_d_txt_c)
    if size(indx{i},1) == 0
        w_temp(i) = NaN;   
    else
        w_temp(i) = mean(r_temp(indx{i}));
    end
end

% pick matched date from air temp date
for i = 1:length(mth_d_txt_c)
       indx_air{i} = find([strcmp(mth_d_txt_c{i}, temp_date_air)] == 1)     
end

% make climate from air temp
for i = 1:length(mth_d_txt_c)
    if size(indx_air{i},1) == 0
        a_temp(i) = NaN;   
    else
        a_temp(i) = nanmean(temp_air(indx_air{i}));
    end
end


pf = polyfit(1:length(a_temp),a_temp,7);
xp = 1:length(a_temp);
yp = polyval(pf,xp);

xp_w = find(isnan(w_temp)==0);
pf_w = polyfit(xp_w,w_temp(xp_w),7);
yp_w = polyval(pf_w,xp);

% water raw - mean_structure (climate)
% rr_temp=r_temp;
for i = 1:length(mth_d_txt_c)
    if size(indx{i},1) == 0
     disp('didnt exist the value')  
    else
        rr_temp(indx{i}) =  r_temp(indx{i}) - yp_w(i);
    end
end

% air raw - daily mean climate
temp_air2 = temp_air;
for i = 1:length(mth_d_txt_c)
    if size(indx_air{i},1) == 0
    disp('didnt exist the value')  
    else
        temp_air2(indx_air{i}) =  temp_air(indx_air{i}) - yp(i);
    end
end

% make water date into yy mm dd form.
for i = 1:length(r_date_txt)
    r_date_txt_yymmdd{i,1} = r_date_txt(i,:);
end

% find(isnan(temp_air)==1)
% interp_t = 1:length(temp_air);
% temp_air(isnan(temp_air))=interp1(interp_t(~isnan(temp_air)),temp_air(~isnan(temp_air)), interp_t(isnan(temp_air))); % filling middle nan;

% find matched date from yy mm dd form.
j=0
for i = 1:length(r_date_txt_yymmdd)
   if  sum(strcmp(r_date_txt_yymmdd{i}, temp_date)) ~= 0
       j=j+1;
       indx_yymmdd(j) = find([strcmp(r_date_txt_yymmdd{i}, temp_date)] == 1);    
   else
       disp(i)
   end
end

%% scatter plot
temp_air_re = temp_air2(indx_yymmdd);

%regression
figure; plot(temp_air_re); hold on; plot(rr_temp,'r');hold off;
return
%slope y = b1*x
b1 = temp_air_re\rr_temp' % x\y for getting slop
yCalc1 = b1*temp_air_re;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(temp_air_re),1) temp_air_re]; %b_0 b_1 
b = X\rr_temp'
yCalc2 = X*b;

% num of samples
clearvars size_in
for i = 1:length(indx)
    size_in(i) = size(indx{i},1);
end

figure;
plot(size_in); hold on;
xlabel('Days','fontsize',13)
ylabel('number of samples','fontsize',13)
title('number of samples (1994~2019)','fontsize',13)
grid on;
xlim([1 366])
legend('number of samples')
set(gca,'fontsize',15)


figure;
scatter(temp_air_re,rr_temp,'.');
hold on
% plot(temp_air_re,yCalc1,'r')
xlabel('air temp. (^oC)','fontsize',13)
ylabel('water temp. (^oC)','fontsize',13)
title('Linear Regression Relation Between air & water temp.')
grid on
plot(temp_air_re,yCalc2,'--','color','r')
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
gtext(['Corr. = ',num2str(corr(temp_air_re,rr_temp'),2)],'Color','k','FontSize',16)
set(gca,'fontsize',13)

figure;
plot(w_temp,'o'); hold on;
plot(a_temp,'r','linew',1.5);
plot(yp,'g','linew',1.5);
plot(yp_w,'c','linew',1.5);
xlabel('Days','fontsize',13)
ylabel('temp. (^oC)','fontsize',13)
title('compare daily climate air & water temp.','fontsize',13)
grid on;
xlim([1 366])
legend('sumjin(Gure)','yeosu-air')
set(gca,'fontsize',15)

save('songjung_regression_yeosu_extract_climate_2021_high_airT.mat'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc; 
cd C:\Users\user\Desktop\장기생태
[raw txt]=xlsread('섬진강_수온_염분_DO.xls','sheet','');
[raw1 txt1]=xlsread('남해_AWS_1990to2018.xls','namhe_1990to1999','');
[raw2 txt2]=xlsread('남해_AWS_1990to2018.xls','namhe_2000to2009','');
[raw3 txt3]=xlsread('남해_AWS_1990to2018.xls','namhe_2010to2018','');
dash_c = '-';
r_txt_ud = flipud(txt);
r_date_txt=[char(r_txt_ud(1:end-2,4)) repmat(dash_c,length(r_txt_ud(1:end-2,4)),1) char(r_txt_ud(1:end-2,5)) ...
    repmat(dash_c,length(r_txt_ud(1:end-2,4)),1) char(r_txt_ud(1:end-2,6))];

r_temp_txt=[r_txt_ud(1:end-2,8)];
r_temp=str2num(char(r_temp_txt));
temp_air=[raw1(:,3); raw2(:,3); raw3(:,3);];
temp_date=[txt1(2:end,2); txt2(2:end,2); txt3(2:end,2);];

for i = 1:length(r_date_txt)
    r_date_txt_c{i,1} = r_date_txt(i,6:end); % delete year
end

char_temp_date=char(temp_date);
for  i = 1:length(temp_date)
    temp_date_air{i,1} = char_temp_date(i,6:end); % delete year from air temp date
end

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
        mth_d_txt(k,:)=[num2str(i,'%02d') '-'  num2str(l,'%02d')]
    end
end

% make it to cell-array
for i = 1:length(mth_d_txt)
    mth_d_txt_c{i,1} = mth_d_txt(i,:); % delete year
end

% pick matched date from water temp date
for i = 1:length(mth_d_txt_c)
       indx{i} = find([strcmp(mth_d_txt_c{i}, r_date_txt_c)] == 1)     
end

% make quasi_climate
for i = 1:length(mth_d_txt_c)
    if size(indx{i},1) == 0
        w_temp(i) = NaN;   
    else
        w_temp(i) = mean(r_temp(indx{i}));
    end
end

% pick matched date from air temp date
for i = 1:length(mth_d_txt_c)
       indx_air{i} = find([strcmp(mth_d_txt_c{i}, temp_date_air)] == 1)     
end

% make climate from air temp
for i = 1:length(mth_d_txt_c)
    if size(indx_air{i},1) == 0
        a_temp(i) = NaN;   
    else
        a_temp(i) = nanmean(temp_air(indx_air{i}));
    end
end


pf = polyfit(1:length(a_temp),a_temp,7);
xp = 1:length(a_temp);
yp = polyval(pf,xp);

xp_w = find(isnan(w_temp)==0);
pf_w = polyfit(xp_w,w_temp(xp_w),7);
yp_w = polyval(pf_w,xp);

% water raw - mean_structure (climate)
rr_temp=r_temp;
for i = 1:length(mth_d_txt_c)
    if size(indx{i},1) == 0
     disp('didnt exist the value')  
    else
        rr_temp(indx{i}) =  r_temp(indx{i}) - yp_w(i);
    end
end

% air raw - daily mean climate
temp_air2 = temp_air;
for i = 1:length(mth_d_txt_c)
    if size(indx_air{i},1) == 0
    disp('didnt exist the value')  
    else
        temp_air2(indx_air{i}) =  temp_air(indx_air{i}) - yp(i);
    end
end

% make water date into yy mm dd form.
for i = 1:length(r_date_txt)
    r_date_txt_yymmdd{i,1} = r_date_txt(i,:);
end

% find(isnan(temp_air)==1)
% interp_t = 1:length(temp_air);
% temp_air(isnan(temp_air))=interp1(interp_t(~isnan(temp_air)),temp_air(~isnan(temp_air)), interp_t(isnan(temp_air))); % filling middle nan;

% find matched date from yy mm dd form.
j=0
for i = 1:length(r_date_txt_yymmdd)
   if  sum(strcmp(r_date_txt_yymmdd{i}, temp_date)) ~= 0
       j=j+1;
       indx_yymmdd(j) = find([strcmp(r_date_txt_yymmdd{i}, temp_date)] == 1)     
   end
end

%% scatter plot
temp_air_re = temp_air2(indx_yymmdd);

%regression
figure; plot(temp_air_re); hold on; plot(rr_temp,'r');hold off;
return
%slope y = b1*x
b1 = temp_air_re\rr_temp % x\y for getting slop
yCalc1 = b1*temp_air_re;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(temp_air_re),1) temp_air_re]; %b_0 b_1 
b = X\rr_temp
yCalc2 = X*b;

figure;
scatter(temp_air_re,rr_temp,'.');
hold on
% plot(temp_air_re,yCalc1,'r')
xlabel('air temp. (^oC)','fontsize',13)
ylabel('water temp. (^oC)','fontsize',13)
title('Linear Regression Relation Between air & water temp.')
grid on
plot(temp_air_re,yCalc2,'--','color','r')
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
gtext(['Corr. = ',num2str(corr(temp_air_re,rr_temp),2)],'Color','k','FontSize',16)
set(gca,'fontsize',13)

figure;
plot(w_temp,'o'); hold on;
plot(a_temp,'r','linew',1.5);
plot(yp,'g','linew',1.5);
plot(yp_w,'c','linew',1.5);
xlabel('Days','fontsize',13)
ylabel('temp. (^oC)','fontsize',13)
title('compare daily climate air & water temp.','fontsize',13)
grid on;
xlim([1 366])
legend('sumjin(hadong)','namhe-air')
set(gca,'fontsize',15)

corr(temp_air_re,rr_temp)
save('sumjin_regression_namhe_extract_climate.mat'); 






