% 수위자료 2000- 2005
% 데이터 이름: 2000.dat
% 첫번째 열 날짜, 1시간 간격 자료
% 화개
% clc; clear all; close all;

% csv read test
% 송정
clc; clear all; close all;



%% read multiple excel file
folder='D:\mepl\data\03_Sumjin_river\songjung\sj_00_16_year\';
filetype='*.xls';  % or xlsx
f=fullfile(folder,filetype);
d=dir(f);
for k=1:numel(d);
  data{k}=xlsread(fullfile(folder,d(k).name));
  temp = data{k};
  temp = temp(:,2:end);
  [a b] = size(temp);
  temp2 = reshape(temp, a*b,1);
  eval(['RD' num2str(1999+k) '= temp2;']) 
  eval(['save RD' num2str(1999+k) ';'])
end
%}

%% River discharge
figure();
% for i = 1:k
%      eval(['h=plot(1:a*b,RD' num2str(1999+i) ');']) 
%      set(h,'linewidth',2,'color',[0.6 0.6 0.6]+i/70);
% hold on;
% end
% ss = [RD2000 RD2001 RD2002 RD2003 RD2004 RD2005 RD2006 RD2007 RD2008 RD2009 RD2010 RD2011 RD2012 RD2013];
ss = [RD2000 RD2001 RD2002 RD2003 RD2004 RD2005 RD2006 RD2007 RD2008 RD2009 RD2010 RD2011 RD2012 RD2013 ...
    RD2014 RD2015 RD2016 ];
plot( nanmean(ss,2),'k-','linewidth',2);
hold on;
for i = 1:k
     eval(['h=plot(1:a*b,RD' num2str(1999+i) ');']) 
     set(h,'linewidth',2,'color',[0.6 0.6 0.6]+i/70);
hold on;
end
plot( nanmean(ss,2),'k-','linewidth',2);
legend('Mean','2000-2013');
ylabel('River discharge (m^3s^-^1)','fontsize',14);
xlabel('Time (Month)','fontsize',14);
xlim([0 365]);
set(gca,'xtick',[1:31:371],'xticklabel',[1:12]);
%% mean 값과 44cms
figure()
plot(nanmean(ss,2),'k-','linewidth',2);
hold on;
% plot(1:5:372,44,'.','color',[0.7 0.7 0.7],'linewidth',1);
xlabel('Time (Month)','fontsize',16);
ylabel('[m^3s^-^1]','fontsize',16);
title('River discharge','fontsize',16);
xlim([1 365]);
set(gca,'xtick',[1:31:371],'xticklabel',[1:12]);
% legend('Mean (2000-2013)','44m^3s^-^1');