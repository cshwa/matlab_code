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
load('agyang_regression_yeosu_extract_climate_2021_high_airT.mat','yp_w','b','r_temp','r_date_txt'); 
% load('agyang_regression_RCP85_IPSL_MR_extract_climate_2100_high_airT.mat','yp',...
%     'yp_w','b','r_temp','r_date_txt','temp_date_noleap_c','temp_date_noleap');  % to 2020.12
load('agyang_regression_RCP85_IPSL_MR_extract_climate_2100_high_airT.mat',...
  'temp_date_noleap_c','temp_date_noleap');  % to 2020.12
load('agyang_regression_yeosu_extract_climate_2020.mat','yp'); % to 2021.05

sjtemp=load('D:\장기생태\Dynamic\06_river_physics\환경과학원\sumjin_recons_water_temp(agyang)_present_2020_high_airT_mixed_daily.mat');
sumjin_re_w_c = sjtemp.merg_recon_w_c;

re_var_mer=reshape(rcp_air.var_mer(:,length(2006:2020)+1:end),365*length(16:95),1);
%% make 2016:2020 clim to make delta 
varrcp_ref=reshape(rcp_air.var_mer(:,length(2006:2015)+1:length(2006:2020)),365*length(11:15),1);
clearvars  varrcp_ref_clim
for i = 1:365
varrcp_ref_clim(i,1)=mean(varrcp_ref(i:365:end));
end

clearvars delta_airT_pre
k=0;
for i = length(2006:2020)+1:length(2006:2100)
k=k+1;
delta_airT_pre(:,k)=squeeze(rcp_air.var_mer(:,i)) - varrcp_ref_clim; %RCP 10yr mean - ref
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
%     if i == 1
%         aws_3reg_delta_clim(:,1) = aws_3reg_clim; % 2016~2020
%         aws_3reg_delta_clim(:,k) = aws_3reg_clim + delta_airT_pre(:,i);
%     else 
%        aws_3reg_delta_clim(:,k) = aws_3reg_clim + delta_airT_pre(:,i);
%     end
        aws_3reg_delta_clim(:,i) = aws_3reg_clim + delta_airT_pre(:,i);
end

k=0;
for i = 27:31
    k=k+1;
temp = sumjin_re_w_c{i};
if length(temp) ~=365
    temp(60)=[];
end
recon_w(:,k) = temp;
end

recon_w_3reg=mean(recon_w,2);

clearvars recon_3reg_delta_clim
k=1;
for i = 1:size(delta_airT_pre,2)
    k=k+1;
%     if i == 1
%         aws_3reg_delta_clim(:,1) = aws_3reg_clim; % 2016~2020
%         aws_3reg_delta_clim(:,k) = aws_3reg_clim + delta_airT_pre(:,i);
%     else 
%        aws_3reg_delta_clim(:,k) = aws_3reg_clim + delta_airT_pre(:,i);
%     end
        recon_3reg_delta_clim(:,i) = recon_w_3reg + delta_airT_pre(:,i);       
end

plot(mean(recon_3reg_delta_clim,1))


aws_3reg_delta=reshape(aws_3reg_delta_clim,size(aws_3reg_delta_clim,1)*size(aws_3reg_delta_clim,2),1);

% mid check 
total_3reg_rcp = [aws_3reg; aws_3reg_delta];
figure; hold on
plot(aws_3reg,'r')
plot(length(aws_3reg)+1:length(aws_3reg_delta)+length(aws_3reg),aws_3reg_delta,'k')
xticks(1:365*10:length(aws_3reg_delta)+length(aws_3reg));
xticklabels(2016:10:2100)
xlim([1 length(aws_3reg_delta)+length(aws_3reg)])
grid on; ylim([-10 35])
ylabel('air temp(^oC)');
set(gca,'fontsize',13);


for i=1:length(2016:2100)
    total_3reg_rcp_yy(i) =mean(total_3reg_rcp(365*(i-1)+1:365*i));
end

%% check delta aws air temp vs. delta recon. water
figure; hold on; 
plot(mean(recon_3reg_delta_clim,1));
plot(mean(aws_3reg_delta_clim,1),'r');
xticks(1:10:length(2016:2100));
xticklabels(2016:10:2100)
ylabel('temp(^oC)');
legend('delta-river','delta-aws');
set(gca,'fontsize',15);
xlim([1 length(2016:2100)])
grid on


figure; 
plot(total_3reg_rcp_yy);
xticks(1:10:length(2016:2100));
xticklabels(2016:10:2100)
ylabel('air temp(^oC)');
set(gca,'fontsize',15);
xlim([1 length(2016:2100)])
grid on


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

clearvars merg_recon_w_c recon_w pre_recon_w X diff_clim
for i = 1:365
    yp_re(i,1) = mean(total_3reg_rcp(i:365:end));
end
yp(60)=[];yp_w(60)=[];
temp_storage=yp;
temp_storage2=yp_w;
% for i=1:length(indx)
k=0; 
for i=1:length(2016:2100)  % 2006~2100
    yp=temp_storage;
    yp_w=temp_storage2;
    k=k+1;
    
%     if i==length(indx)
%         yp(length(temp_air(indx{i}:end))+1:end)=[];
%         yp_w(length(temp_air(indx{i}:end))+1:end)=[];
%     end
%     
%     diff_clim = total_3reg_rcp(365*(i-1)+1:365*i) - yp'; % - climate
    diff_clim = total_3reg_rcp(365*(i-1)+1:365*i) - yp_re'; % - climate

% ['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)]
X = [ones(length(diff_clim'),1) diff_clim']; %b_0 b_1
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

tx = 1:length(total_3reg_rcp);
for i = 1:length(indx); indx_num(i)=indx{i}; end

for i=1:length(2016:2100)
    merg_recon_w_yy(i) =mean(merg_recon_w(365*(i-1)+1:365*i));
end

figure; hold on
plot(merg_recon_w_yy);
plot(total_3reg_rcp_yy,'r');
xticks(1:10:length(2016:2100));
xticklabels(2016:10:2100)
ylabel('air temp(^oC)');
set(gca,'fontsize',15);
xlim([1 length(2016:2100)])
grid on


figure; hold on
% plot(temp_air,'r','linew',1.5);
% plot(tx(indx_put),r_temp,'r','linew',1.5);
plot(tx,merg_recon_w,'b','linew',1.5);
% scatter(tx(indx_put),r_temp,'r','linew',1.5);
xlabel('Days','fontsize',13)
ylabel('temp. (^oC)','fontsize',13)
title('recons. water temp','fontsize',13)
grid on;
xlim([1 length(total_3reg_rcp)]);
set(gca,'xtick',1:365*5:length(total_3reg_rcp));
set(gca,'xticklabel',2016:5:2100);
% legend({'sumjin-water(agyang)','recon. water'})
set(gca,'fontsize',15)
print('-dpng',['recons_sumjin_2006~2100']); hold off;
xlim([1 length(total_3reg_rcp)]);
ylim([-5 30]); 
xtickangle(45)
