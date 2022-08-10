clear; clc; close all

% station_list = {'거문도', '거제도', '고흥발포', '목포', '부산', '여수', '완도', '진도', '통영'};
station_list = {'여수'};
year_list = [2012:2015];
figpath = 'E:\11.사업\장기생태_2단계\2차년도_2015\03_kickoffmeeting\준비\tide_station\조위관측소 수온\';

for si = 1:length(station_list)
    
    
station = cell2mat(station_list(si))

fpath = 'E:\11.사업\장기생태_2단계\2차년도_2015\03_kickoffmeeting\준비_데이터처리\tide_station\';
%flist = dir(fullfile(fpath, '*.txt'));
fname = [station, '_2000-01_2015-09.txt'];
file = [fpath, fname];

[num yyyymmdd hh temp] = textread(file, '%d %s %s %f', 'headerlines', 6);
yyyymmdd = cell2mat(yyyymmdd);
yyyy = yyyymmdd(:, 1:4);
mm = yyyymmdd(:,6:7);
dd = yyyymmdd(:, 9:10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fmmdd = str2num([mm, dd]);
fdata = [fmmdd, temp];
fdatenum = datenum(2012, str2num(mm), str2num(dd));

fmmdd_u = unique(fmmdd);
fdata_mean = [];
for i = 1:length(fmmdd_u)
    ind = find(fmmdd == fmmdd_u(i));
    fdata_mean = [fdata_mean; mean(fdatenum(ind)) mean(fdata(ind, 2))];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for yi = 1:length(year_list)
    
yi2 = find(str2num(yyyy) == year_list(yi));
yyyy2 = yyyy(yi2, :);
mm2 = mm(yi2, :);
dd2 = dd(yi2, :);
temp2 = temp(yi2, :);

fmmdd2 = str2num([mm2, dd2]);
fdata2 = [fmmdd2, temp2];
fdatenum2 = datenum(2012, str2num(mm2), str2num(dd2));

fmmdd_u2 = unique(fmmdd2);
fdata_mean2 = [];
for ii = 1:length(fmmdd_u2)
    ind2 = find(fmmdd2 == fmmdd_u2(ii));
    fdata_mean2 = [fdata_mean2; mean(fdatenum2(ind2)) mean(fdata2(ind2, 2))];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold on;
plot(fdata_mean(:,1), fdata_mean(:,2), '.k', 'MarkerSize', 15)
plot(fdata_mean2(:,1), fdata_mean2(:,2), '.r', 'MarkerSize', 15)

%---- 은별 수정 부분 ------------
l = length(fdata_mean2);
ss(yi,1:l) = fdata_mean2(:,2);
% --- 은별 수정 부분 끝 ---------

datetick('x', 'mm')
xlabel('Month'); ylabel('Temperature')
set(gca, 'fontsize', 15)
ylim([0 30])
title(['Mean vs. ',yyyy2(1,:), ' ', station])
legend([yyyy(1,:), '-', yyyy(end,:)], yyyy2(1,:), 'location', 'NorthWest')

saveas(gcf, [figpath, station, '_', num2str(year_list(yi)), '.png'])

end
end