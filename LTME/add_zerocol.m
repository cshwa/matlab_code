% 1 coloum data add method
% m-file: add_zerocol.m
% inputdata: busan0511.dat
% output data: busan0511_out.dat

% load raw data
data = textread('busan0511.dat','','headerlines',4);   % ��ü���� load
% data name and data
year = data(:,1);
month = data(:,2);
day = data(:,3);
hour = data(:,4);
t_hight = data(:,5);   % tide height

n = size(data,1);          % count a line
zerocol = zeros(n,1);      % add this column (n by 1)

% ���ο� �������� �ٲٰ� ����
Ndata = [year,month,day,hour,zerocol,t_hight];  %Ndata ������� ����
fid = fopen('busan0511_new.dat','w');
fprintf(fid,'%4d%3d%3d%3d%3d%5d\n',Ndata');
fclose(fid)