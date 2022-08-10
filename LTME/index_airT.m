clc; clear all; close all;

%--- PDO ----------------
temp = load('PDO.latest.txt'); 
%--- 1995년부터  2015년 까지
data = temp(96:end,2:end);
[a b] = size(data);
ss =reshape(data,[],1);
l = length(ss);
%
FigHandle = figure;
set(FigHandle, 'Position', [100, 200, 750, 200]);
plot(ss,'linewidth',2);
set(gca, 'xtick',[1:12*5:l],'xticklabel',[1995:5:2015]);
xlabel('Time (year)','fontsize',14);
ylabel('PDO','fontsize',14);
xlim([1 l]);

%% --- AO ----------
temp = load('monthly.ao.index.b50.current.ascii.txt'); 
%--- 1995년부터  2015년 까지
data = temp(541:792,3);
ss = data;
l = length(ss);
FigHandle = figure;
set(FigHandle, 'Position', [100, 200, 750, 200]);
plot(ss,'linewidth',2);
set(gca, 'xtick',[1:12*5:l],'xticklabel',[1995:5:2015]);
xlabel('Time (year)','fontsize',14);
ylabel('AO','fontsize',14);
xlim([1 l]);
ylim([-5 5]);
%% --- MEI ------------
temp = load('mei_index.txt'); 
%--- 1995년부터  2015년 까지
data = temp(46:end,2:end);
[a b] = size(data);
ss =reshape(data,[],1);
l = length(ss);
%
FigHandle = figure;
set(FigHandle, 'Position', [100, 200, 750, 200]);
plot(ss,'linewidth',2);
set(gca, 'xtick',[1:12*5:l],'xticklabel',[1995:5:2015]);
xlabel('Time (year)','fontsize',14);
ylabel('MEI','fontsize',14);
xlim([1 l]);


