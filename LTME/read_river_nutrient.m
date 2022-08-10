% ������ ��� ���翰
clc; clear all; close all;

filename = '�ϵ�_��õ����������.xlsx';
sname = '�ϵ�'; % �ϵ�, ����
[num,txt,raw] = xlsread(filename,sname);

ph = num(:,1);
DO = num(:,2);
BOD = num(:,3);
COD = num(:,4);
ss = num(:,5);
total_n = num(:,6);
total_p = num(:,7);

DateStrings = txt(3:end,1);
t = datetime(DateStrings,'InputFormat','yyyy-MM-dd');
DateNumber = datenum(t);

yData = total_n;
%%
figure
plot(t,yData,'linewidth',1.5)
ylabel('������ (ml/l)','fontsize',16);
xlabel('�ð�(����)','fontsize',16);
startDate = datenum('01-01-1980');
endDate = datenum('12-31-2018');
xlim([startDate endDate]);
set(gca,'fontsize',14);
%% --- save figure --------------------------------------
outpath = '.\figure\';
out_vari = sname;
out=[outpath,out_vari];
set(gcf,'renderer','painter');
set(gcf, 'PaperUnits', 'inches');
x_width=10;
y_width=2.5;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
saveas(gcf,out,'tif');
%------------------------------------------------------


% startDate = datenum('01-01-1980');
% endDate = datenum('12-31-2018');
% xData = linspace(startDate,endDate,50);
% 
% datetick('x','yyyy','keeplimits')
% dateFormat = 12;
% datetick('x',dateFormat)
