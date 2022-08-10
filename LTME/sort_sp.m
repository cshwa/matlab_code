

close all
clear all
clc


FILE = 'FlowQuest600_PC_1.xlsx';
[NUMERIC,TXT,RAW]=xlsread(FILE);

data = NUMERIC;
[n id] = size(data);

sp = zeros(n, id/2);
dir = zeros(n, id/2);

n=0;
m=0;
for i = 1:id
    if mod(i,2)==1
        n = n+1;
        sp(:,n) = data(:,(n-1)*2+1);
    else
        m=m+1;
        dir(:,m) = data(:,m*2);
    end
end

dir = 450 - dir;
dir = mod(dir,360);

u = sp.*cos(dir*pi/180);
v = sp.*sin(dir*pi/180);

FILE = 'tide_PC1.xlsx';
[NUMERIC,TXT,RAW]=xlsread(FILE);

ele = NUMERIC;

ind = round(ele-6.0);
s_lay = 20;
uus2=[]; vvs2=[];
for i = 1:length(u)
    us(i) = u(i,ind(i));
    vs(i) = v(i,ind(i));
    
    uu = squeeze(u(i,1:ind(i)));
    vv = squeeze(v(i,1:ind(i)));
    dd = (ind(i)-1)/(s_lay-1);
    uus = interp1(1:ind(i),uu,1:dd:ind(i));
    vvs = interp1(1:ind(i),vv,1:dd:ind(i));
    uus2 = [uus2; uus]; vvs2 = [vvs2; vvs];
end



% ub = u(:,1);
% vb = v(:,1);
% plot(us)
% hold on
% plot(vs,'r')

% figure
% plot(ub)
% hold on
% plot(vb,'r')
% data = [u(:,20) v(:,20) ele];
% data = data';
% fid = fopen('sur_uv_ele.dat','w');
% fprintf(fid,'%10.2f%10.2f%10.3f\n',data);
% fclose(fid);
% 
% data = [u(:,1) v(:,1) ele];
% data = data';
% fid = fopen('bot_uv_ele.dat','w');
% fprintf(fid,'%10.2f%10.2f%10.3f\n',data);
% fclose(fid);

% plot(u(:,1),'.r')
% hold on
% plot(u(:,20),'-r')
% hold on
% 
% plot(v(:,1),'.b')
% hold on
% plot(v(:,20),'-b')
% hold on
% plot(ele,'k')








