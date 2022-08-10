close all; clear; clc; 
cd C:\Users\user\Desktop\장기생태
[raw txt]=xlsread('가화천_수온_염분_DO.xls','sheet','');
[raw1 txt1]=xlsread('사천만_수위_2000to2002.xls','사천만_수위_2000to2002','');
[raw2 txt2]=xlsread('사천만_수위_2003to2005.xls','사천만_수위_2003to2005','');
[raw3 txt3]=xlsread('사천만_수위_2006to2008.xls','사천만_수위_2006to2008','');
[raw4 txt4]=xlsread('사천만_수위_2009to2011.xls','사천만_수위_2009to2011','');
[raw5 txt5]=xlsread('사천만_수위_2012to2014.xls','사천만_수위_2012to2014','');
[raw6 txt6]=xlsread('사천만_수위_2015to2017.xls','사천만_수위_2015to2017','');
[raw7 txt7]=xlsread('사천만_수위_2018.xls','사천만_수위_2018','');
[raw11 txt11]=xlsread('사천만_수위_1987to1989.xls','사천만_수위_1987to1989','');
[raw22 txt22]=xlsread('사천만_수위_1990to1992.xls','사천만_수위_1990to1992','');
[raw33 txt33]=xlsread('사천만_수위_1993to1995.xls','사천만_수위_1993to1995','');
[raw44 txt44]=xlsread('사천만_수위_1996to1998.xls','사천만_수위_1996to1998','');
[raw55 txt55]=xlsread('사천만_수위_1999.xls','사천만_수위_1999','');

% dash_c = '-';
% r_txt_ud = flipud(txt);
% r_date_txt=[char(r_txt_ud(1:end-2,4)) repmat(dash_c,length(r_txt_ud(1:end-2,4)),1) char(r_txt_ud(1:end-2,5)) ...
%     repmat(dash_c,length(r_txt_ud(1:end-2,4)),1) char(r_txt_ud(1:end-2,6))];

% r_temp_txt=[r_txt_ud(1:end-2,8)];
% r_temp=str2num(char(r_temp_txt));
temp_elev=[flipud(raw1(:,1)); flipud(raw2(:,1)); flipud(raw3(:,1)); flipud(raw4(:,1)); flipud(raw5(:,1)); flipud(raw6(:,1)); ...
    flipud(raw7(:,1));];  %temp_elev is present Namgang dem
temp_date=[flipud(txt1(3:end,1)); flipud(txt2(3:end,1)); flipud(txt3(3:end,1)); flipud(txt4(3:end,1)); flipud(txt5(3:end,1)); flipud(txt6(3:end,1)); ...
    flipud(txt7(3:end,1));];

temp_elev_bf=[flipud(raw11(:,1)); flipud(raw22(:,1)); flipud(raw33(:,1)); flipud(raw44(:,1)); flipud(raw55(:,1));];  %temp_elev is present Namgang dem
temp_date_bf=[flipud(txt11(3:end,1)); flipud(txt22(3:end,1)); flipud(txt33(3:end,1)); flipud(txt44(3:end,1)); flipud(txt55(3:end,1));];

temp_elev(find(temp_elev == -999)) =NaN; %-999 is NaN 
temp_elev(find(temp_elev >= 50)) =NaN;
temp_elev_bf(find(temp_elev_bf >= 50)) = NaN;

% interp_t = 1:length(temp_elev);
% temp_elev(isnan(temp_elev))=interp1(interp_t(~isnan(temp_elev)),temp_elev(~isnan(temp_elev)), interp_t(isnan(temp_elev))); % filling middle nan;

plot(temp_elev);hold on; plot(temp_elev_bf,'r')
line(1:length(temp_elev),repmat(50,length(temp_elev),1),'color','k')

merge_elev= [temp_elev_bf; temp_elev;]; 

for j=1987:2018
    if j == 1987
        pre_t_ax_mth(j-1986)=sum(eomday(j,6:12));
    else
        pre_t_ax_mth(j-1986)=sum(eomday(j,1:12));
    end
end

for i = 1:length(1987:2018)
t_ax_mth(i) = sum(pre_t_ax_mth(1:i));
end

% non - missing time stamp
zz =0 
for k = 1987:2018
for j = 1:12
for i = 1:eomday(k,j)
    zz=zz+1;
    t_ind(zz,:)=[num2str(k),'-', num2str(j,'%02d'), '-', num2str(i,'%02d')];
end
end
end
t_ind(1:151,:)=[]; % 1987 start 06.01

%make it cell array
for i = 1:length(t_ind)
t_ind_c{i}= t_ind(i,:);
end

merge_date=[temp_date; temp_date_bf;];

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
set(gca,'xticklabel',1988:2018);
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
