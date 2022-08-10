close all; clear; clc; 
cd D:\장기생태\Dynamic\06_river_physics
[rawa1 txta1]=xlsread('진주_AWS_1990to2019.xls','jinju_1990to1999','');
[rawa2 txta2]=xlsread('진주_AWS_1990to2019.xls','jinju_2000to2009','');
[rawa3 txta3]=xlsread('진주_AWS_1990to2019.xls','jinju_2010to2019','');
[rawa4 txta4]=xlsread('진주_AWS_2020.xlsx','진주_AWS_2020','');

cd D:\장기생태\Dynamic\06_river_physics
[raw txt]=xlsread('남강댐_환경부_정보.xls','sheet1','');
[raw_n1 txt_n1]=xlsread('남강댐_2018.xls','수질(일반측정망)','');
[raw_n2 txt_n2]=xlsread('남강댐_2019.xls','수질(일반측정망)','');
[raw_n3 txt_n3]=xlsread('남강댐_2020.xls','수질(일반측정망)','');

dash_c = '-';

n_2018 =char(txt_n1(3:end,3)); % 2018 date
n_2019 =char(txt_n2(3:end,3)); % 2018 date
n_2020 =char(txt_n3(3:end,3)); % 2018 date
n_2018_ymd = [n_2018(22:end,1:4) repmat(dash_c,length(n_2018(22:end,1:4)),1) n_2018(22:end,6:7) repmat(dash_c,length(n_2018(22:end,1:4)),1) n_2018(22:end,9:10)];
n_2019_ymd = [n_2019(:,1:4) repmat(dash_c,length(n_2019(:,1:4)),1) n_2019(:,6:7) repmat(dash_c,length(n_2019(:,1:4)),1) n_2019(:,9:10)];
n_2020_ymd = [n_2020(:,1:4) repmat(dash_c,length(n_2020(:,1:4)),1) n_2020(:,6:7) repmat(dash_c,length(n_2020(:,1:4)),1) n_2020(:,9:10)];
 
r_txt_ud = flipud(txt);
r_date_txt=[char(r_txt_ud(1:end-2,4)) repmat(dash_c,length(r_txt_ud(1:end-2,4)),1) char(r_txt_ud(1:end-2,5)) ...
    repmat(dash_c,length(r_txt_ud(1:end-2,4)),1) char(r_txt_ud(1:end-2,6))];
r_date_txt = [r_date_txt; n_2018_ymd; n_2019_ymd; n_2020_ymd;]; % concaternate (merge)

r_temp_txt=[r_txt_ud(1:end-2,8)];
r_temp=str2num(char(r_temp_txt));
r_temp = [r_temp; raw_n1(22:end,14); raw_n2(:,14); raw_n3(:,14);]; % concaternate (merge)
temp_air=[rawa1(:,3); rawa2(:,3); rawa3(:,3); rawa4(:,4);];
temp_date=[txta1(2:end,2); txta2(2:end,2); txta3(2:end,2); txta4(2:end,3);]; %air temp date

for i = 1:length(r_date_txt)
    r_date_txt_c{i,1} = r_date_txt(i,6:end); % delete year
end

char_temp_date=char(temp_date);
for  i = 1:length(temp_date)
    temp_date_air{i,1} = char_temp_date(i,6:end); % delete year from air temp date
end

for i = 1:length(r_date_txt)
    r_date_txt_yymm{i,1} = r_date_txt(i,1:7); % delete day
end

% make reference yymm (only for 1mth)
j=0
for i = str2num(r_date_txt(1,1:4)) : 2019 
    j=j+1;
    axe_ref_yymm{j,1} = [num2str(i), '-','01'];
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

% pick matched date from water temp date to make time axis
for i = 1:length(axe_ref_yymm)
       indx_yymm{i} = find([strcmp(axe_ref_yymm{i}, r_date_txt_yymm)] == 1)     
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
legend('sumjin(hadong)','yeosu-air')
set(gca,'fontsize',15)



% figure; hold on;
%     yp=temp_storage;
%     yp_w=temp_storage2;
%     if leapyear(1979+i) == 0
%         yp(60) = [];
%         yp_w(60) = [];
%     end
% diff_clim = temp_air(indx{i}:indx{i+1}-1) - yp'; % - climate
% % ['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)]
% X = [ones(length(diff_clim),1) diff_clim]; %b_0 b_1
% pre_recon_w=X*b
% recon_w = pre_recon_w + yp_w'; % + climate
% plot(recon_w,':');
% plot(w_temp,'o');
% plot(a_temp,'r','linew',1.5);
% plot(yp,'g','linew',1.5);
% plot(yp_w,'c','linew',1.5);
% xlabel('Days','fontsize',13)
% ylabel('temp. (^oC)','fontsize',13)
% title('compare daily climate air & water temp.','fontsize',13)
% grid on;
% xlim([1 366])
% legend('gawha-water','jinju-air')
% set(gca,'fontsize',15)
save('gawha_regression_extract_climate_2020.mat'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc; 
cd D:\장기생태\Dynamic\06_river_physics
[rawa1 txta1]=xlsread('여수_AWS_1990to2019.xls','yeosu_1990to1999','');
[rawa2 txta2]=xlsread('여수_AWS_1990to2019.xls','yeosu_2000to2009','');
[rawa3 txta3]=xlsread('여수_AWS_1990to2019.xls','yeosu_2010to2019','');
[rawa4 txta4]=xlsread('여수_AWS_2020.xlsx','여수_AWS_2020','');

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

temp_air=[rawa1(:,3); rawa2(:,3); rawa3(:,3); rawa4(:,4);];
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

save('agyang_regression_yeosu_extract_climate_2020.mat'); 

