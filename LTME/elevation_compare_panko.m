close all; clear; clc; 
cd C:\Users\user\Desktop\장기생태
[raw txt]=xlsread('가화천_수온_염분_DO.xls','sheet','');
[raw1 txt1]=xlsread('판문_수위_2001to2003.xls','판문_수위_2001to2003','');
[raw2 txt2]=xlsread('판문_수위_2004to2006.xls','판문_수위_2004to2006','');
[raw3 txt3]=xlsread('판문_수위_2007to2009.xls','판문_수위_2007to2009','');
[raw4 txt4]=xlsread('판문_수위_2010to2012.xls','판문_수위_2010to2012','');
[raw5 txt5]=xlsread('판문_수위_2013to2015.xls','판문_수위_2013to2015','');

% dash_c = '-';
% r_txt_ud = flipud(txt);
% r_date_txt=[char(r_txt_ud(1:end-2,4)) repmat(dash_c,length(r_txt_ud(1:end-2,4)),1) char(r_txt_ud(1:end-2,5)) ...
%     repmat(dash_c,length(r_txt_ud(1:end-2,4)),1) char(r_txt_ud(1:end-2,6))];

% r_temp_txt=[r_txt_ud(1:end-2,8)];
% r_temp=str2num(char(r_temp_txt));

temp_elev=[flipud(raw1(:,1)); flipud(raw2(:,1)); flipud(raw3(:,1)); flipud(raw4(:,1)); flipud(raw5(:,1));];  %temp_elev is present Namgang dem
temp_date=[flipud(txt1(3:end,1)); flipud(txt2(3:end,1)); flipud(txt3(3:end,1)); flipud(txt4(3:end,1)); flipud(txt5(3:end,1)); ];

temp_elev(find(temp_elev == -999)) =NaN; %-999 is NaN 
temp_elev(find(temp_elev >= 50)) =NaN;


% interp_t = 1:length(temp_elev);
% temp_elev(isnan(temp_elev))=interp1(interp_t(~isnan(temp_elev)),temp_elev(~isnan(temp_elev)), interp_t(isnan(temp_elev))); % filling middle nan;

plot(temp_elev);hold on; plot(temp_elev_bf,'r')
line(1:length(temp_elev),repmat(50,length(temp_elev),1),'color','k')

merge_elev= temp_elev; 

for j=2001:2015
    if j == 2001
        pre_t_ax_mth(j-2000)=sum(eomday(j,6:12));
    else
        pre_t_ax_mth(j-2000)=sum(eomday(j,1:12));
    end
end

for i = 1:length(2001:2015)
t_ax_mth(i) = sum(pre_t_ax_mth(1:i));
end

% non - missing time stamp
zz =0 
for k = 2001:2015
for j = 1:12
for i = 1:eomday(k,j)
    zz=zz+1;
    t_ind(zz,:)=[num2str(k),'-', num2str(j,'%02d'), '-', num2str(i,'%02d')];
end
end
end
% t_ind(1:151,:)=[]; % 1987 start 06.01

%make it cell array
for i = 1:length(t_ind)
t_ind_c{i}= t_ind(i,:);
end

merge_date= temp_date; ;

% find matched date
j=0
k=0
for i = 1:length(t_ind_c)
   if  sum(strcmp(t_ind_c{i}, merge_date)) == 0
       j=j+1;   
        indx(j) = i;  %not matched index 
   end
   
   if  sum(strcmp(t_ind_c{i}, merge_date)) ~= 0
       k=k+1;   
        indx_match(k) = i; %matched index 
   end
   
end
% t_ind(indx,:) %not matched date

% make matrix including NaNs
merge_elev_fix=NaN(length(t_ind_c),1);
merge_elev_fix(indx_match) =  merge_elev;

% replacing NaN
interp_t = 1:length(merge_elev_fix);
merge_elev_fix(isnan(merge_elev_fix))=interp1(interp_t(~isnan(merge_elev_fix)),merge_elev_fix(~isnan(merge_elev_fix)), interp_t(isnan(merge_elev_fix))); % filling middle nan;

t_ax = [[1]; [t_ax_mth(1:end-1)+1]';]

return

figure;
plot(merge_elev_fix); hold on;
xlabel('year','fontsize',13)
ylabel('dam elev. (m)','fontsize',13)
title('Dam elevation.(Sacheon)')
grid on
line(1:length(merge_elev_fix),repmat(50,length(merge_elev_fix),1),'color','r')
set(gca,'xtick',t_ax(2:end));
set(gca,'xlim',[1 t_ax_mth(end)]);
set(gca,'xticklabel',2002:2015);
% gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
set(gca,'fontsize',13,'fontweight','bold')




return

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
