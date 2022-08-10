 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  gahwa river %%%%%%%%%%
% close all; clear; clc; 
% cd C:\Users\user\Desktop\장기생태
% [raw1 txt1]=xlsread('진주_AWS_1990to2018.xls','jinju_1990to1999','');
% [raw2 txt2]=xlsread('진주_AWS_1990to2018.xls','jinju_2000to2009','');
% [raw3 txt3]=xlsread('진주_AWS_1990to2018.xls','jinju_2010to2018','');
% load gawha_river_temp_obs_00to18.mat
% 
% temp=[raw1(:,3); raw2(:,3); raw3(:,3);];
% 
% find(isnan(temp)==1)
% 
% interp_t = 1:length(temp);
% temp(isnan(temp))=interp1(interp_t(~isnan(temp)),temp(~isnan(temp)), interp_t(isnan(temp))); % filling middle nan;

close all; clear; clc; 
cd C:\Users\user\Desktop\장기생태
[raw txt]=xlsread('가화천_수온_염분_DO.xls','sheet','');
[raw1 txt1]=xlsread('진주_AWS_1990to2018.xls','jinju_1990to1999','');
[raw2 txt2]=xlsread('진주_AWS_1990to2018.xls','jinju_2000to2009','');
[raw3 txt3]=xlsread('진주_AWS_1990to2018.xls','jinju_2010to2018','');
dash_c = '-';
r_txt_ud = flipud(txt);
r_date_txt=[char(r_txt_ud(1:end-2,4)) repmat(dash_c,length(r_txt_ud(1:end-2,4)),1) char(r_txt_ud(1:end-2,5)) ...
    repmat(dash_c,length(r_txt_ud(1:end-2,4)),1) char(r_txt_ud(1:end-2,6))];

r_temp_txt=[r_txt_ud(1:end-2,8)];
r_temp=str2num(char(r_temp_txt));
temp_air=[raw1(:,3); raw2(:,3); raw3(:,3);];
temp_date=[txt1(2:end,2); txt2(2:end,2); txt3(2:end,2);];

for i = 1:length(r_date_txt)
    r_date_txt_c{i,1} = r_date_txt(i,:);
end

find(isnan(temp_air)==1)
interp_t = 1:length(temp_air);
temp_air(isnan(temp_air))=interp1(interp_t(~isnan(temp_air)),temp_air(~isnan(temp_air)), interp_t(isnan(temp_air))); % filling middle nan;

% find matched date
j=0
for i = 1:length(r_date_txt_c)
   if  sum(strcmp(r_date_txt_c{i}, temp_date)) ~= 0
       j=j+1;
       indx(j) = find([strcmp(r_date_txt_c{i}, temp_date)] == 1)     
   end
end

temp_date(indx)

%regression
temp_air_re = temp_air(indx);
figure; plot(temp_air_re); hold on; plot(r_temp,'r');hold off;

%slope y = b1*x
b1 = temp_air_re\r_temp % x\y for getting slop
yCalc1 = b1*temp_air_re;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(temp_air_re),1) temp_air_re]; %b_0 b_1 
b = X\r_temp
yCalc2 = X*b;

figure;
scatter(temp_air_re,r_temp);
hold on
% plot(temp_air_re,yCalc1,'r')
xlabel('air temp. (^oC)','fontsize',13)
ylabel('water temp. (^oC)','fontsize',13)
title('Linear Regression Relation Between air & water temp.')
grid on
plot(temp_air_re,yCalc2,'--','color','r')
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
set(gca,'fontsize',13)

save('gawha_regression.mat','temp_air_re', 'r_temp','b','yCalc2'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc; 
cd C:\Users\user\Desktop\장기생태
[raw txt]=xlsread('섬진강_수온_염분_DO.xls','sheet','');
[raw1 txt1]=xlsread('여수_AWS_1990to2018.xls','yeosu_1990to1999','');
[raw2 txt2]=xlsread('여수_AWS_1990to2018.xls','yeosu_2000to2009','');
[raw3 txt3]=xlsread('여수_AWS_1990to2018.xls','yeosu_2010to2018','');
dash_c = '-';
r_txt_ud = flipud(txt);
r_date_txt=[char(r_txt_ud(1:end-2,4)) repmat(dash_c,length(r_txt_ud(1:end-2,4)),1) char(r_txt_ud(1:end-2,5)) ...
    repmat(dash_c,length(r_txt_ud(1:end-2,4)),1) char(r_txt_ud(1:end-2,6))];

r_temp_txt=[r_txt_ud(1:end-2,8)];
r_temp=str2num(char(r_temp_txt));
temp_air=[raw1(:,3); raw2(:,3); raw3(:,3);];
temp_date=[txt1(2:end,2); txt2(2:end,2); txt3(2:end,2);];

for i = 1:length(r_date_txt)
    r_date_txt_c{i,1} = r_date_txt(i,:);
end

find(isnan(temp_air)==1)
interp_t = 1:length(temp_air);
temp_air(isnan(temp_air))=interp1(interp_t(~isnan(temp_air)),temp_air(~isnan(temp_air)), interp_t(isnan(temp_air))); % filling middle nan;

% find matched date
j=0
for i = 1:length(r_date_txt_c)
   if  sum(strcmp(r_date_txt_c{i}, temp_date)) ~= 0
       j=j+1;
       indx(j) = find([strcmp(r_date_txt_c{i}, temp_date)] == 1)     
   end
end

temp_date(indx)

%regression
temp_air_re = temp_air(indx);
figure; plot(temp_air_re); hold on; plot(r_temp,'r');hold off;

%slope y = b1*x
b1 = temp_air_re\r_temp % x\y for getting slop
yCalc1 = b1*temp_air_re;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(temp_air_re),1) temp_air_re]; %b_0 b_1 
b = X\r_temp
yCalc2 = X*b;

figure;
scatter(temp_air_re,r_temp);
hold on
% plot(temp_air_re,yCalc1,'r')
xlabel('air temp. (^oC)','fontsize',13)
ylabel('water temp. (^oC)','fontsize',13)
title('Linear Regression Relation Between air & water temp.')
grid on
plot(temp_air_re,yCalc2,'--','color','r')
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
set(gca,'fontsize',13)

save('songjung_regression_yeosu.mat','temp_air_re', 'r_temp','b','yCalc2'); 
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
    r_date_txt_c{i,1} = r_date_txt(i,:);
end

find(isnan(temp_air)==1)
interp_t = 1:length(temp_air);
temp_air(isnan(temp_air))=interp1(interp_t(~isnan(temp_air)),temp_air(~isnan(temp_air)), interp_t(isnan(temp_air))); % filling middle nan;

% find matched date
j=0
for i = 1:length(r_date_txt_c)
   if  sum(strcmp(r_date_txt_c{i}, temp_date)) ~= 0
       j=j+1;
       indx(j) = find([strcmp(r_date_txt_c{i}, temp_date)] == 1)     
   end
end

temp_date(indx)

%regression
temp_air_re = temp_air(indx);
figure; plot(temp_air_re); hold on; plot(r_temp,'r');hold off;

%slope y = b1*x
b1 = temp_air_re\r_temp % x\y for getting slop
yCalc1 = b1*temp_air_re;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(temp_air_re),1) temp_air_re]; %b_0 b_1 
b = X\r_temp
yCalc2 = X*b;

figure;
scatter(temp_air_re,r_temp);
hold on
% plot(temp_air_re,yCalc1,'r')
xlabel('air temp. (^oC)','fontsize',13)
ylabel('water temp. (^oC)','fontsize',13)
title('Linear Regression Relation Between air & water temp.')
grid on
plot(temp_air_re,yCalc2,'--','color','r')
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
set(gca,'fontsize',13)
% save('songjung_regression_yeosu.mat','temp_air_re',
%'r_temp','b','yCalc2');  don't need to save it (it's not the best fit
% set)










