%===================================================================
% ����������п�(KODC) �������� �ڷ� - ���縸 �߽� �ڷ� ����ȭ (����)
% sorted_data �� �ڷ�(���ϴ� ����)���� ���ϴ� �⵵ syr = 19840000 ���� ������
% ���δ� �⵵, ���δ� ¦���޷� ���ϴ� ���ɿ��� ���� �̾Ƴ���
% 
% copywrite EunByeol Cho 2015.05.24
%===================================================================

clc; clear all;
close all

%--- ������ ���� ���� �������� ---------------------------------------------
[ST_info] = textread('NSOposition_silver.txt','%s %*[^\n]','headerlines',1);
ST_info = str2mat(ST_info); 
ST_info = ST_info(:,[1:3,5:6]);
ST_info = str2num(ST_info);
fname = ST_info;
%--- data load Ư�� ������ -------
% fpath = ['.\sorted_data\'];
% fname = ['20500'];
%================================
m = length(fname);
syr = 19840000; % �̾Ƴ��� ���� ���� �⵵
for ss = 1:1 % 60: 40012
    fpath = ['.\sorted_data\'];
    temp = load([fpath, num2str(fname(ss)),'_re.txt']);
    temp = temp(2:end,:);
    temp = sortrows(temp,2);
    [a b]=find(temp(:,2)>=syr);
    temp = temp(min(a):end,:);
    
%======================================================
%--- ���ϴ� ���� �ڷ� �̱� -----------------------------
p_dep = 5; % 5 = sst
[a b] = find(temp(:,5) <= p_dep);
data = temp(a,:);

%--- ��Һ� ������ ���� --------------------------------
st = data(:,1); % station number
date = data(:,2);  
lon = data(:,3); % location
lat = data(:,4);
dep = data(:,5);
tempe = data(:,6);
salt = data(:,7);
%======================================================

%--- ���� �ڷḦ �̱� ���� ��¥ �� ó�� -----------------
ye = fix(date/10000);
mmdd = mod(date,10000);
mm = fix(mmdd/100);
dd = mod(mmdd,100);
%======================================================

red(:,1) = st;
red(:,2) = ye;
red(:,3) = mm;
red(:,9) = dd;
red(:,4) = lon;
red(:,5) = lat;
red(:,6) = dep;
red(:,7) = tempe;
red(:,8) = salt;

%--- total year mean -------------------------------------------
    for k = 1:length(ye)
        S = num2str(date(k));
        D(k,1:10) = [S(5:6), '/', S(7:8), '/', S(1:4)];
        NumDays(k,1) = days365(['01/01/',S(1:4)], D(k,1:10));
    end

x = NumDays;
y = tempe;

%--- poly fit ����ϴ� ��� -------------
% p = polyfit(x,y,3);
% y1 = polyval(p,x); % year mean
%=======================================
%--- Spline fit ����ϴ� ��� -----------
SplineFit = fit(x, y, 'smoothingspline','SmoothingParam',0.00015);
y1 = feval(SplineFit,x);

ym_red = red;
% ym_red(:,7) = red(:,7)-y1;
red_10d = ym_red; % ����� �� sst�� �� �ִ� �ڷ� 

%--- resorted -------------------------------------------------------------
re_red = red_10d(:,1:8);
ly = max(ye)-min(ye)+1;
sor = nan(ly+1,7);
sor(1:ly,1) = [min(ye):max(ye)];
l = length(re_red);
    for i = 1:l
        sor((re_red(i,2)-(min(ye)-1)),re_red(i,3)) = re_red(i,7);
    end

re_sort = sor(:,[1,2:2:12]);
re_sort(re_sort == 0) = NaN;
re_sort(1:ly,1) = [min(ye):max(ye)];
%==========================================================================

outpath = '.\sorted_yearly\';
out=[outpath,num2str(fname(ss))];
save([out '_1984.dat'],'re_sort','-ascii');
% save([out '_ye_trend.dat'],'y','-ascii');
% clear red re_sort sor x y p NumDays
% close; 

%%

x = re_sort(:,1);
y1 = re_sort(:,2);
y2 = re_sort(:,3);
y3 = re_sort(:,4);
y4 = re_sort(:,5);
y5 = re_sort(:,6);
y6 = re_sort(:,7);
%--- linear fit --------------------------
[a b]=find(~isnan(y1)); x1 = x(a); ya1 = y1(a);
P = polyfit(x1,ya1,1); f1 = P(1)*x1+P(2); 
GD(ss,1) = P(1);
[a b]=find(~isnan(y2)); x2 = x(a); ya2 = y2(a);
P = polyfit(x2,ya2,1); f2 = P(1)*x2+P(2);
GD(ss,2) = P(1);
[a b]=find(~isnan(y3)); x3 = x(a); ya3 = y3(a);
P = polyfit(x3,ya3,1); f3 = P(1)*x3+P(2);
GD(ss,3) = P(1);
[a b]=find(~isnan(y4)); x4 = x(a); ya4 = y4(a);
P = polyfit(x4,ya4,1); f4 = P(1)*x4+P(2);
GD(ss,4) = P(1);
[a b]=find(~isnan(y5)); x5 = x(a); ya5 = y5(a);
P = polyfit(x5,ya5,1); f5 = P(1)*x5+P(2);
GD(ss,5) = P(1);
[a b]=find(~isnan(y6)); x6 = x(a); ya6 = y6(a);
P = polyfit(x6,ya6,1); f6 = P(1)*x6+P(2);
GD(ss,6) = P(1);
%============================================

figure('Position', [100, 100, 600, 400])
hold on;
sz = 25;
scatter(x,y1,sz,[0.2 0.0 0.6],'filled'); 
scatter(x,y2,sz,[65 105 225]/255,'filled');
scatter(x,y3,sz,[34 139 34]/255,'filled'); 
scatter(x,y4,sz,[255 215 0]/255,'filled'); 
scatter(x,y5,sz,[255 97 3]/255,'filled');
scatter(x,y6,sz,[64 224 208]/255,'filled');

plot(x1,f1,'color',[0.2 0.0 0.6]', 'LineWidth', 2 );
plot(x2,f2,'color',[65 105 225]/255', 'LineWidth', 2 );
plot(x3,f3,'color',[34 139 34]/255', 'LineWidth', 2 );
plot(x4,f4,'color',[255 215 0]/255', 'LineWidth', 2 );
plot(x5,f5,'color',[255 97 3]/255', 'LineWidth', 2 );
plot(x6,f6,'color',[64 224 208]/255', 'LineWidth', 2 );

xlabel('Time (year)','fontsize',14);
ylabel('Temperature (0m)','fontsize',14);
legend('2��','4��','6��','8��','10��','12��','location','south','Orientation','horizontal')
% axis([min(ye)-1 max(ye)+1 -5 5]);
t_name=[num2str(fname(ss)),'  ',num2str(p_dep), ' m depth'];
title(t_name,'fontweight','bold','fontsize',14,'color','b');

%--- save figure --------------------
outpath = '.\figure\';
out=[outpath,num2str(fname(ss)),'_lf_',num2str(p_dep)];
set(gcf,'renderer','painter');
set(gcf, 'PaperUnits', 'inches');
x_width=6;
y_width=4;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
saveas(gcf,out,'tif');
%------------------------------------

clear red re_sort sor x y p NumDays
close; 
%}
end

%% ���� ��¥ Ȯ��
%{
z = ye;
response = z;

% Make a color index for the ozone levels
nc =max(ye)-min(ye)+1;
offset = 1;
c = response - min(response);
% c = round((nc-1-2*offset)*c/max(c)+1+offset);
% c = [nc:-1:1];
% Create a 3D scatter plot using the scatter3 function
figure
scatter3(mm, dd, ye, 30, c, 'filled')
view(-34, 14)

hcb=colorbar;
set(hcb,'YTick',[1:5:48],'yticklabel',[1965:5:2015]);

xlabel('Time (month)','fontsize',14);
ylabel('Time (Day)','fontsize',14);
title(t_name,'fontweight','bold','fontsize',14,'color','b');
ylim([0 31]);
az = 0;
el = 90;
view(az, el);
%--- save figure --------------------
% outpath = '.\figure\';
% out=[outpath,num2str(fname(ss)),'day_',num2str(p_dep)];
% set(gcf,'renderer','painter');
% set(gcf, 'PaperUnits', 'inches');
% x_width=6;
% y_width=4;
% set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
% saveas(gcf,out,'tif');
%------------------------------------
%}