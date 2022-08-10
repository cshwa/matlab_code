clc;clear all;close all;

%     [f, p]=uigetfile('adcp.xlsx');
%         filedir=[p,f];
%         disp(filedir);
    data=xlsread('adcp.xlsx');data=data(8:end,11:end);
    data=data(:,2:end);

% ------------------- 1시간간격 moving average 한 것 그려보았으나, 타원모양은 나타나지 않음.
%     mean_data=data;
%     for j=1:1:length(data(1,:))
%     for i=4:1:length(data(:,1))-2
%         mean_data(i,j)=mean(data(i-3:i+2,j));
%     end
%     end
up_data=[];
for i=1:2:21
  j=i+1 ;
%  plot(data(:,i),data(:,j),'r.')
%  feather(data(:,2),data(:,3))
  [cc]=princomp([data(:,i),data(:,j)]);
  if cc(1,1) > 0
    cc=[cc(1,1),cc(1,2)*(-1);cc(2,1)*(-1),cc(2,2)];
  else
    cc=[cc(1,1),cc(1,2)*(-1);cc(2,1)*(-1),cc(2,2)]*(-1);
  end
  p_data=cc*[data(:,i),data(:,j)]';
  up_data=[up_data,p_data(1,:)'];
end

for i=1:1:length(up_data(1,:))
  up_mean_data(i,1)=mean(up_data(1:144,i));
end

plot(up_mean_data,[1:1:11])