%===================================================================
% 국립수산과학원(KODC) 정선관측 자료 - 광양만 중심 자료 가시화 (남해)
% 
% 
% 
% copywrite EunByeol Cho 2015.05.24
%===================================================================

clc; clear all; close all

%--- 정선별 정점 정보 가져오기 ---------------------------------------------
[ST_info] = textread('NSOposition_silver.txt','%s %*[^\n]','headerlines',1);
ST_info = str2mat(ST_info); ST_info = ST_info(:,[1:3,5:6]);
ST_info = str2num(ST_info);
fname = ST_info;
%--- data load ------------------------------------
% fpath = ['.\sorted_data\'];
% fname = ['20500'];
m = length(fname);

for ss = 1:length(ST_info)
% for ss = 12:12
    fpath = ['.\sorted_data\'];
    temp = load([fpath, num2str(fname(ss)),'_re.txt']);
%======================================================

%--- 원하는 수심 자료 뽑기 -----
p_dep = 0;
[a b] = find(temp(:,5) <= p_dep+3);
data = temp(a,:);

%--- 요소별 데이터 정리 --------------------------------
st = data(:,1);
date = data(:,2);
lon = data(:,3);
lat = data(:,4);
dep = data(:,5);
temp = data(:,6);
salt = data(:,7);
%======================================================

%--- 월별 자료를 뽑기 위해 날짜 열 처리 -----------------
ye = fix(date/10000);
mmdd = mod(date,10000);
mm = fix(mmdd/100);
dd = mod(mmdd,100);
%======================================================

%--- red: resorted 될 행렬 ---------
red(:,1) = st;
red(:,2) = ye;
red(:,3) = mm;
red(:,4) = lon;
red(:,5) = lat;
red(:,6) = dep;
red(:,7) = temp;
red(:,8) = salt;

%--- resorted ----------------------
ly = max(ye)-min(ye)+1;
sor = nan(ly+1,7);
sor(1:ly,1) = [min(ye):max(ye)];
l = length(data);

    for i = 1:l
        sor((red(i,2)-(min(ye)-1)),red(i,3)) = red(i,7);
    end

re_sort = sor(:,[1,2:2:12]);
re_sort(re_sort == 0) = NaN;
re_sort(1:ly,1) = [min(ye):max(ye)];
%%
x = re_sort(:,1);
y1 = re_sort(:,2);
y2 = re_sort(:,3);
y3 = re_sort(:,4);
y4 = re_sort(:,5);
y5 = re_sort(:,6);
y6 = re_sort(:,7);

figure('Position', [100, 100, 600, 400])
hold on;
sz = 25;
scatter(x,y1,sz,[0.2 0.0 0.6],'filled'); 
scatter(x,y2,sz,[65 105 225]/255,'filled');
scatter(x,y3,sz,[34 139 34]/255,'filled'); 
scatter(x,y4,sz,[255 97 32]/255,'filled'); 
scatter(x,y5,sz,[255 215 0]/255,'filled');
scatter(x,y6,sz,[64 224 208]/255,'filled');
ls = lsline;
set(ls(6),'color',[0.2 0.0 0.6]', 'LineWidth', 2 );
set(ls(5),'color',[65 105 225]/255', 'LineWidth', 2 );
set(ls(4),'color',[34 139 34]/255', 'LineWidth', 2 );
set(ls(3),'color',[255 97 3]/255', 'LineWidth', 2 );
set(ls(2),'color',[255 215 0]/255', 'LineWidth', 2 );
set(ls(1),'color',[64 224 208]/225', 'LineWidth', 2 );

xlabel('Time (year)','fontsize',14);
ylabel('Temperature (0m)','fontsize',14);
legend('2월','4월','6월','8월','10월','12월','location','south','Orientation','horizontal')
axis([min(ye)-1 max(ye)+1 0 35]);
t_name=[num2str(fname(ss)),'  ',num2str(p_dep), ' m depth'];
title(t_name,'fontweight','bold','fontsize',14,'color','b');

%--- save figure --------------------
outpath = '.\figure\';
% out_vari = fname(i);
out=[outpath,num2str(fname(ss)),'_',num2str(p_dep)];
set(gcf,'renderer','painter');
set(gcf, 'PaperUnits', 'inches');
x_width=6;
y_width=4;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
saveas(gcf,out,'tif');
%------------------------------------

clear red re_sort sor
close
end

%% 관측 날짜 확인
%{
figure()
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
%}
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