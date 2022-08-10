%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc; 
cd D:\장기생태\Dynamic\06_river_physics
% [raw1 txt1]=xlsread('1980s_yeosu_airT.csv','1980s_yeosu_airT',''); % 2: date, 3: mean temp
[raw1 txt1]=xlsread('진주_AWS_1990to2019.xls','jinju_1990to1999','');
[raw2 txt2]=xlsread('진주_AWS_1990to2019.xls','jinju_2000to2009','');
[raw3 txt3]=xlsread('진주_AWS_1990to2019.xls','jinju_2010to2019','');
[rawa4 txta4]=xlsread('진주_AWS_2020.xlsx','진주_AWS_2020','');
cd  D:\장기생태\Dynamic\06_river_physics\환경과학원
load('gawha_regression_extract_climate_2020_high_airT.mat','yp_w','b','r_temp','r_date_txt'); 
load('gawha_regression_extract_climate_2020.mat','yp'); 

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
    r_date_txt_c{i,1} = r_date_txt(i,:); % delete year
end

char_temp_date=char(temp_date);
for  i = 1:length(temp_date)
    temp_date_air{i,1} = char_temp_date(i,6:end); % delete year from air temp date
end

% find matched date
j=0
for i = 1:length(r_date_txt_c)
   if  sum(strcmp(r_date_txt_c{i}, temp_date_c)) ~= 0
       j=j+1;
       indx_put(j) = find([strcmp(r_date_txt_c{i}, temp_date_c)] == 1)     
   end
end
input_w = NaN(length(temp_date_c),1);
input_w(indx_put) = r_temp;

% make reference yymm (only for 1mth)
j=0
for i = str2num(temp_date(1,1:4)) : str2num(temp_date(end,1:4)) 
    j=j+1;
    axe_ref_mm{j,1} = [num2str(i), '-','01'];
end

% pick matched date from water temp date
ref_ddmm{1,1} = ['01-01']
j=0
for i = 1:length(temp_date_air)
        if [strcmp(temp_date_air{i}, ref_ddmm)] == 1
            j=j+1;
        indx{j} = i;
        end
end

temp_storage=yp;
temp_storage2=yp_w;
for i=1:length(indx)
    yp=temp_storage;
    yp_w=temp_storage2;
    if leapyear(1989+i) == 0
        yp(60) = [];
        yp_w(60) = [];
    end
    
if i == length(indx) 
    diff_clim = temp_air(indx{i}:end) - yp'; % - climate
else
    diff_clim = temp_air(indx{i}:indx{i+1}-1) - yp'; % - climate
end

% ['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)]
X = [ones(length(diff_clim),1) diff_clim]; %b_0 b_1
pre_recon_w=X*b
recon_w = pre_recon_w + yp_w'; % + climate
% plot(recon_w)

merg_recon_w_c{i} = recon_w;

if i == 1
    merg_recon_w = recon_w;
else
    merg_recon_w = [merg_recon_w; recon_w;]; 
end
end

tx = 1:length(temp_air);
for i = 1:length(indx); indx_num(i)=indx{i}; end


figure; hold on
% plot(temp_air,'r','linew',1.5);
plot(tx,merg_recon_w,'c','linew',1.5);
scatter(tx(indx_put),r_temp,'r','linew',1.5);
xlabel('Days','fontsize',13)
ylabel('temp. (^oC)','fontsize',13)
title('compare recons. water temp','fontsize',13)
grid on;
xlim([1 length(temp_air)]);
set(gca,'xtick',indx_num);
set(gca,'xticklabel',1990:2019);
legend({'recon. water','gahwa-water'})
set(gca,'fontsize',15)
print('-dpng',['recons_gawha_1990~2020']); hold off;
save('gawha_recons_water_temp_present_2020_high_airT_mixed_daily.mat');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc; 
cd D:\장기생태\Dynamic\06_river_physics
% [raw1 txt1]=xlsread('1980s_yeosu_airT.csv','1980s_yeosu_airT',''); % 2: date, 3: mean temp
[rawa1 txta1]=xlsread('여수_AWS_1990to2019.xls','yeosu_1990to1999','');
[rawa2 txta2]=xlsread('여수_AWS_1990to2019.xls','yeosu_2000to2009','');
[rawa3 txta3]=xlsread('여수_AWS_1990to2019.xls','yeosu_2010to2019','');
[rawa4 txta4]=xlsread('여수_AWS_2020.xlsx','여수_AWS_2020','');
% [rawa5 txta5]=xlsread('D:\RIST\2021_RIST_광양만\KMA_ASOS\여수_AWS_2021_08.xlsx','Sheet1','');
rcp_air = load('J:\장기생태_2021\Dynamic\result\CMIP5\tair_test66_2006-2100.mat'); 
cd  D:\장기생태\Dynamic\06_river_physics\환경과학원
% load('agyang_regression_yeosu_extract_climate_2021_high_airT.mat','yp','yp_w','b','r_temp','r_date_txt'); 
load('agyang_regression_RCP85_IPSL_MR_extract_climate_2100_high_airT.mat','yp',...
    'yp_w','b','r_temp','r_date_txt','temp_date_noleap_c','temp_date_noleap');  % to 2020.12
% load('agyang_regression_yeosu_extract_climate_2020.mat','yp'); % to 2021.05

re_var_mer=reshape(rcp_air.var_mer(:,length(2006:2020)+1:end),365*length(16:95),1);
%% make 2016:2020 clim to make delta 
varrcp_ref=reshape(rcp_air.var_mer(:,length(2006:2015)+1:length(2006:2020)),365*length(11:15),1);
clearvars  varrcp_ref_clim
for i = 1:365
varrcp_ref_clim(i,1)=mean(varrcp_ref(i:365:end));
end

clearvars delta_airT_pre
k=0;
for i = length(2006:2020)+1:10:length(2006:2100)
k=k+1;
delta_airT_pre(:,k)=squeeze(mean(rcp_air.var_mer(:,i:i+9),2)) - varrcp_ref_clim; %RCP 10yr mean - ref
end

% plot(delta_airT_pre)

re_var_tot=reshape(rcp_air.var_mer(:,:),365*length(1:95),1);
temp_air_aws=[rawa1(:,3); rawa2(:,3); rawa3(:,3); rawa4(:,4);]; % AWS 1990~2020
re_var_tot=reshape(rcp_air.var_mer(:,:),365*length(1:95),1);

dash_c = '-';
% temp_air=[rawa1(:,6); rawa2(:,6); rawa3(:,6); rawa4(:,7); ];
temp_air=[re_var_tot;];
temp_date_c=temp_date_noleap_c;
temp_date=char(temp_date_noleap_c);

% dash_c = '-';
% temp_air=[raw1(:,3); raw2(:,3); raw3(:,3);];
temp_date_c=[txta1(2:end,2); txta2(2:end,2); txta3(2:end,2); txta4(2:end,3);]; %air temp date
% temp_date=char(temp_date_c);
% {'2016-01-01'} == 9497
aws_3reg=temp_air_aws(9497:length(temp_date_c),1);
aws_3reg(60)=[];
% 365*4+1:length(aws_3reg)
aws_3reg( 365*4+60) = [];

clearvars aws_3reg_clim
for i = 1:365
    aws_3reg_clim(i,1)=mean(aws_3reg(i:365:end));
end

clearvars aws_3reg_delta_clim
k=1;
for i = 1:size(delta_airT_pre,2)
    k=k+1;
    if i == 1
        aws_3reg_delta_clim(:,1) = aws_3reg_clim; % 2016~2020
        aws_3reg_delta_clim(:,k) = aws_3reg_clim + delta_airT_pre(:,i);
    else 
       aws_3reg_delta_clim(:,k) = aws_3reg_clim + delta_airT_pre(:,i);
    end
end

for i = 1:length(r_date_txt)
    r_date_txt_c{i,1} = r_date_txt(i,:); % delete year
end

char_temp_date=char(temp_date);
for  i = 1:length(temp_date)
    temp_date_air{i,1} = char_temp_date(i,6:end); % delete year from air temp date
end


% find matched date
j=0
for i = 1:length(r_date_txt_c)
   if  sum(strcmp(r_date_txt_c{i}, temp_date_c)) ~= 0
       j=j+1;
       indx_put(j) = find([strcmp(r_date_txt_c{i}, temp_date_c)] == 1);     
   else
       disp(i)
   end
end
% input_w = NaN(length(temp_date_c),1);
input_w(indx_put) = r_temp;

% make reference yymm (only for jan.)
j=0
for i = str2num(temp_date(1,1:4)) : str2num(temp_date(end,1:4)) 
    j=j+1;
    axe_ref_mm{j,1} = [num2str(i), '-','01'];
end

% pick matched date from water temp date
ref_ddmm{1,1} = ['01-01']
j=0
for i = 1:length(temp_date_air)
        if [strcmp(temp_date_air{i}, ref_ddmm)] == 1
            j=j+1;
        indx{j} = i;
        end
end

temp_storage=yp;
temp_storage2=yp_w;
% for i=1:length(indx)
k=0; 
for i=1:length(indx)  % 2006~2100
    yp=temp_storage;
    yp_w=temp_storage2;
    k=k+1;
    
%     if i <= 31
%         if leapyear(1989+i) == 0
%             yp(60) = [];
%             yp_w(60) = [];
%         end
%     else %2021:2100 has no leap year on cmip5 RCP
%         yp(60) = [];
%         yp_w(60) = [];
%     end
    
    if i==length(indx)
        yp(length(temp_air(indx{i}:end))+1:end)=[];
        yp_w(length(temp_air(indx{i}:end))+1:end)=[];
    end
    
if i < length(indx)
    diff_clim = temp_air(indx{i}:indx{i+1}-1) - yp'; % - climate
elseif i == length(indx)
    diff_clim = temp_air(indx{i}:end) - yp'; % - climate
end

% ['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)]
X = [ones(length(diff_clim),1) diff_clim]; %b_0 b_1
pre_recon_w=X*b
recon_w = pre_recon_w + yp_w'; % + climate
% plot(recon_w)

merg_recon_w_c{i} = recon_w;

if k == 1
    merg_recon_w = recon_w;
else
    merg_recon_w = [merg_recon_w; recon_w;]; 
end
end

tx = 1:length(temp_air);
for i = 1:length(indx); indx_num(i)=indx{i}; end


figure; hold on
% plot(temp_air,'r','linew',1.5);
plot(tx(indx_put),r_temp,'r','linew',1.5);
plot(tx,merg_recon_w,'b','linew',1.5);
scatter(tx(indx_put),r_temp,'r','linew',1.5);
xlabel('Days','fontsize',13)
ylabel('temp. (^oC)','fontsize',13)
title('compare recons. water temp','fontsize',13)
grid on;
xlim([1 length(temp_air)]);
set(gca,'xtick',indx_num(1:5:end));
set(gca,'xticklabel',2006:5:2100);
legend({'sumjin-water(agyang)','recon. water'})
set(gca,'fontsize',15)
print('-dpng',['recons_sumjin_2006~2100']); hold off;
xlim([1 length(temp_air)]);
ylim([-5 30]); 
xtickangle(45)


%% yearly mean
% for i = 1:length(1990:2100)
for i = 1:length(2006:2100)
recon_temp_yr(i)=mean(merg_recon_w(365*(i-1)+1:365*i));
end
figure; hold on
% plot(temp_air,'r','linew',1.5);
% plot(tx(indx_put),r_temp,'r','linew',1.5);
plot(recon_temp_yr,'b','linew',1.5);
% scatter(tx(indx_put),r_temp,'r','linew',1.5);
xlabel('Years','fontsize',13)
ylabel('temp. (^oC)','fontsize',13)
title('compare recons. water temp','fontsize',13)
grid on;
xlim([1 length(temp_air)]);
set(gca,'xtick',1:5:length(2006:2100));
set(gca,'xticklabel',2006:5:2100);
legend({'sumjin-water(agyang)','recon. water'})
set(gca,'fontsize',15)
print('-dpng',['recons_sumjin_2021~2019']); hold off;
xlim([1 length(recon_temp_yr)]);
ylim([14 18]); 
xtickangle(45)

save('sumjin_recons_water_temp(agyang)_present_2100_high_airT_mixed_daily.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc; 
cd D:\장기생태\Dynamic\06_river_physics
% [raw1 txt1]=xlsread('1980s_yeosu_airT.csv','1980s_yeosu_airT',''); % 2: date, 3: mean temp
% [raw1 txt1]=xlsread('여수_AWS_1990to2019.xls','yeosu_1990to1999','');
% [raw2 txt2]=xlsread('여수_AWS_1990to2019.xls','yeosu_2000to2009','');
[raw3 txt3]=xlsread('여수_AWS_1990to2019.xls','yeosu_2010to2019','');
load('songjung_regression_yeosu_extract_climate_2019.mat','yp','yp_w','b','r_temp','r_date_txt'); 

dash_c = '-';
temp_air=[ raw3(:,3)];
temp_date_c=[ txt3(2:end,2) ]; %air temp date
temp_date=char(temp_date_c);


r_date_txt(1:305,:) =[];
r_temp(1:305,:) =[];
for i = 1:length(r_date_txt)
    r_date_txt_c{i,1} = r_date_txt(i,:); % delete year
end

char_temp_date=char(temp_date);
for  i = 1:length(temp_date)
    temp_date_air{i,1} = char_temp_date(i,6:end); % delete year from air temp date
end

% find matched date
j=0
for i = 1:length(r_date_txt_c)
   if  sum(strcmp(r_date_txt_c{i}, temp_date_c)) ~= 0
       j=j+1;
       indx_put(j) = find([strcmp(r_date_txt_c{i}, temp_date_c)] == 1)     
   end
end
input_w = NaN(length(temp_date_c),1);
input_w(indx_put) = r_temp;

% make reference yymm (only for 1mth)
j=0
for i = str2num(temp_date(1,1:4)) : str2num(temp_date(end,1:4)) 
    j=j+1;
    axe_ref_mm{j,1} = [num2str(i), '-','01'];
end

% pick matched date from water temp date
ref_ddmm{1,1} = ['01-01']
j=0
for i = 1:length(temp_date_air)
        if [strcmp(temp_date_air{i}, ref_ddmm)] == 1
            j=j+1;
        indx{j} = i;
        end
end

temp_storage=yp;
temp_storage2=yp_w;
for i=1:length(indx)
    yp=temp_storage;
    yp_w=temp_storage2;
    if leapyear(1989+i) == 0
        yp(60) = [];
        yp_w(60) = [];
    end
    
if i == length(indx) 
    diff_clim = temp_air(indx{i}:end) - yp'; % - climate
else
    diff_clim = temp_air(indx{i}:indx{i+1}-1) - yp'; % - climate
end

% ['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)]
X = [ones(length(diff_clim),1) diff_clim]; %b_0 b_1
pre_recon_w=X*b
recon_w = pre_recon_w + yp_w'; % + climate
% plot(recon_w)

merg_recon_w_c{i} = recon_w;

if i == 1
    merg_recon_w = recon_w;
else
    merg_recon_w = [merg_recon_w; recon_w;]; 
end
end

tx = 1:length(temp_air);
for i = 1:length(indx); indx_num(i)=indx{i}; end

figure; hold on
% plot(temp_air,'r','linew',1.5);
plot(tx,merg_recon_w,'b','linew',1.5);
plot(tx(indx_put),r_temp,'r','linew',1.5);
scatter(tx(indx_put),r_temp,'r','linew',1.5);
xlabel('Days','fontsize',13)
ylabel('temp. (^oC)','fontsize',13)
title('compare recons. water temp','fontsize',13)
grid on;
xlim([1 length(temp_air)]);
set(gca,'xtick',indx_num);
set(gca,'xticklabel',2010:2019);
legend({'recon. water','sumjin-water'})
set(gca,'fontsize',15)
print('-dpng',['recons_sumjin_2010~2019']); hold off;
save('sumjin_recons_water_temp_present_2010_2019.mat');




