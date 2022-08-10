% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  gahwa river %%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all; clear; clc; 
% cd C:\Users\user\Desktop\장기생태 
% % [raw1 txt1]=xlsread('1980s_yeosu_airT.csv','1980s_yeosu_airT',''); % 2: date, 3: mean temp
% [raw2 txt2]=xlsread('1980s_jinju_airT.csv','1980s_jinju_airT','');
% load('gawha_regression_extract_climate.mat','yp','yp_w','b'); 
% 
% dash_c = '-';
% for i =2:size(txt2,1); temp_date{i-1,1}=txt2{i,2}; end
% temp_date=char(temp_date);
% temp_air=raw2(:,3);
% 
% char_temp_date=char(temp_date);
% for  i = 1:length(temp_date)
%     temp_date_air{i,1} = char_temp_date(i,6:end); % delete year from air temp date
% end
% 
% % make reference yymm (only for 1mth)
% j=0
% for i = str2num(temp_date(1,1:4)) : str2num(temp_date(end,1:4)) 
%     j=j+1;
%     axe_ref_mm{j,1} = [num2str(i), '-','01'];
% end
% 
% % pick matched date from water temp date
% ref_ddmm{1,1} = ['01-01']
% j=0
% for i = 1:length(temp_date_air)
%         if [strcmp(temp_date_air{i}, ref_ddmm)] == 1
%             j=j+1;
%         indx{j} = i;
%         end
% end
% 
% temp_storage=yp;
% temp_storage2=yp_w;
% for i=1:length(indx)
%     yp=temp_storage;
%     yp_w=temp_storage2;
%     if leapyear(1979+i) == 0
%         yp(60) = [];
%         yp_w(60) = [];
%     end
% if i == length(indx) 
%     diff_clim = temp_air(indx{i}:end) - yp'; % - climate
% else
%     diff_clim = temp_air(indx{i}:indx{i+1}-1) - yp'; % - climate
% end
% % ['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)]
% X = [ones(length(diff_clim),1) diff_clim]; %b_0 b_1
% pre_recon_w=X*b
% recon_w = pre_recon_w + yp_w'; % + climate
% 
% figure; hold on
% if i == length(indx) 
%     plot(temp_air(indx{i}:end),'r','linew',1.5);
% else
%     plot(temp_air(indx{i}:indx{i+1}-1),'r','linew',1.5);
% end
% plot(yp,'g','linew',1.5);
% plot(yp_w,'b','linew',1.5);
% plot(recon_w,'c','linew',1.5);
% xlabel('Days','fontsize',13)
% ylabel('temp. (^oC)','fontsize',13)
% title('compare daily climate air & water temp.','fontsize',13)
% grid on;
% xlim([1 366])
% legend({'jinju-air','climate-air','climate-water','recon. water'},'position',[.52,.25,.1,.2])
% set(gca,'fontsize',15)
% print('-dpng',['recons_gawha_' num2str(i+1979)]); hold off;
% end
% 
% for i=1:length(indx)
%     yp=temp_storage;
%     yp_w=temp_storage2;
%     if leapyear(1979+i) == 0
%         yp(60) = [];
%         yp_w(60) = [];
%     end
%     
% if i == length(indx) 
%     diff_clim = temp_air(indx{i}:end) - yp'; % - climate
% else
%     diff_clim = temp_air(indx{i}:indx{i+1}-1) - yp'; % - climate
% end
% 
% % ['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)]
% X = [ones(length(diff_clim),1) diff_clim]; %b_0 b_1
% pre_recon_w=X*b
% recon_w = pre_recon_w + yp_w'; % + climate
% % plot(recon_w)
% 
% merg_recon_w_c{i} = recon_w;
% 
% if i == 1
%     merg_recon_w = recon_w;
% else
%     merg_recon_w = [merg_recon_w; recon_w;]; 
% end
% end
% save('gawha_recons_water_temp_1980s.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc; 
cd D:\장기생태\Dynamic\06_river_physics
% [raw1 txt1]=xlsread('1980s_yeosu_airT.csv','1980s_yeosu_airT',''); % 2: date, 3: mean temp
[raw1 txt1]=xlsread('진주_AWS_1990to2019.xls','jinju_1990to1999','');
[raw2 txt2]=xlsread('진주_AWS_1990to2019.xls','jinju_2000to2009','');
[raw3 txt3]=xlsread('진주_AWS_1990to2019.xls','jinju_2010to2019','');
load('gawha_regression_extract_climate_2019.mat','yp','yp_w','b','r_temp','r_date_txt'); 

dash_c = '-';
temp_air=[raw1(:,3); raw2(:,3); raw3(:,3);];
temp_date_c=[txt1(2:end,2); txt2(2:end,2); txt3(2:end,2);]; %air temp date
temp_date=char(temp_date_c);

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
print('-dpng',['recons_gawha_1990~2019']); hold off;
save('gawha_recons_water_temp_present_2019.mat');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  sumjin %%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all; clear; clc; 
% cd C:\Users\user\Desktop\장기생태 
% % [raw1 txt1]=xlsread('1980s_yeosu_airT.csv','1980s_yeosu_airT',''); % 2: date, 3: mean temp
% [raw2 txt2]=xlsread('1980s_yeosu_airT.csv','1980s_yeosu_airT','');
% load('songjung_regression_yeosu_extract_climate.mat','yp','yp_w','b'); 
% 
% dash_c = '-';
% for i =2:size(txt2,1); temp_date{i-1,1}=txt2{i,2}; end
% temp_date=char(temp_date);
% temp_air=raw2(:,3);
% 
% char_temp_date=char(temp_date);
% for  i = 1:length(temp_date)
%     temp_date_air{i,1} = char_temp_date(i,6:end); % delete year from air temp date
% end
% 
% % make reference yymm (only for 1mth)
% j=0
% for i = str2num(temp_date(1,1:4)) : str2num(temp_date(end,1:4)) 
%     j=j+1;
%     axe_ref_mm{j,1} = [num2str(i), '-','01'];
% end
% 
% % pick matched date from water temp date
% ref_ddmm{1,1} = ['01-01']
% j=0
% for i = 1:length(temp_date_air)
%         if [strcmp(temp_date_air{i}, ref_ddmm)] == 1
%             j=j+1;
%         indx{j} = i;
%         end
% end
% 
% temp_storage=yp;
% temp_storage2=yp_w;
% for i=1:length(indx)
%     yp=temp_storage;
%     yp_w=temp_storage2;
%     if leapyear(1979+i) == 0
%         yp(60) = [];
%         yp_w(60) = [];
%     end
% if i == length(indx) 
%     diff_clim = temp_air(indx{i}:end) - yp'; % - climate
% else
%     diff_clim = temp_air(indx{i}:indx{i+1}-1) - yp'; % - climate
% end
% % ['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)]
% X = [ones(length(diff_clim),1) diff_clim]; %b_0 b_1
% pre_recon_w=X*b
% recon_w = pre_recon_w + yp_w'; % + climate
% 
% figure; hold on
% if i == length(indx) 
%     plot(temp_air(indx{i}:end),'r','linew',1.5);
% else
%     plot(temp_air(indx{i}:indx{i+1}-1),'r','linew',1.5);
% end
% plot(yp,'g','linew',1.5);
% plot(yp_w,'b','linew',1.5);
% plot(recon_w,'c','linew',1.5);
% xlabel('Days','fontsize',13)
% ylabel('temp. (^oC)','fontsize',13)
% title('compare daily climate air & water temp.','fontsize',13)
% grid on;
% xlim([1 366])
% legend({'yeosu-air','climate-air','climate-water','recon. water'},'position',[.52,.25,.1,.2])
% set(gca,'fontsize',15)
% print('-dpng',['recons_sumjin_' num2str(i+1979)]); hold off;
% end
% 
% 
% for i=1:length(indx)
%     yp=temp_storage;
%     yp_w=temp_storage2;
%     if leapyear(1979+i) == 0
%         yp(60) = [];
%         yp_w(60) = [];
%     end
%     
% if i == length(indx) 
%     diff_clim = temp_air(indx{i}:end) - yp'; % - climate
% else
%     diff_clim = temp_air(indx{i}:indx{i+1}-1) - yp'; % - climate
% end
% 
% % ['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)]
% X = [ones(length(diff_clim),1) diff_clim]; %b_0 b_1
% pre_recon_w=X*b
% recon_w = pre_recon_w + yp_w'; % + climate
% % plot(recon_w)
% 
% merg_recon_w_c{i} = recon_w;
% 
% if i == 1
%     merg_recon_w = recon_w;
% else
%     merg_recon_w = [merg_recon_w; recon_w;]; 
% end
% end
% save('sumjin_recons_water_temp_1980s.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc; 
cd D:\장기생태\Dynamic\06_river_physics
% [raw1 txt1]=xlsread('1980s_yeosu_airT.csv','1980s_yeosu_airT',''); % 2: date, 3: mean temp
[raw1 txt1]=xlsread('여수_AWS_1990to2019.xls','yeosu_1990to1999','');
[raw2 txt2]=xlsread('여수_AWS_1990to2019.xls','yeosu_2000to2009','');
[raw3 txt3]=xlsread('여수_AWS_1990to2019.xls','yeosu_2010to2019','');
load('songjung_regression_yeosu_extract_climate_2019.mat','yp','yp_w','b','r_temp','r_date_txt'); 

dash_c = '-';
temp_air=[raw1(:,3); raw2(:,3); raw3(:,3);];
temp_date_c=[txt1(2:end,2); txt2(2:end,2); txt3(2:end,2);]; %air temp date
temp_date=char(temp_date_c);

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
plot(tx(indx_put),r_temp,'r','linew',1.5);
plot(tx,merg_recon_w,'b','linew',1.5);
scatter(tx(indx_put),r_temp,'r','linew',1.5);
xlabel('Days','fontsize',13)
ylabel('temp. (^oC)','fontsize',13)
title('compare recons. water temp','fontsize',13)
grid on;
xlim([1 length(temp_air)]);
set(gca,'xtick',indx_num);
set(gca,'xticklabel',1990:2019);
legend({'recon. water','sumjin-water'})
set(gca,'fontsize',15)
print('-dpng',['recons_sumjin_1990~2019']); hold off;
xlim([1466 length(temp_air)]);
ylim([-1 35]); xtickangle(45)
save('sumjin_recons_water_temp_present_2019.mat');


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




