% ����data�� ������ ǥ���ϴ� ������� Scatter plot, histogram, rose���� Ȱ��
% m-file: tidalcurrent_exam1.m
% ����: ���� rose �Լ����� x���ذ����� ���ϱ��ذ����� �ٲ��� ����

clc;clear all;close all
% data load
[yy,mm,dd,hr,min,temp,speed,n_dir,u_comp,v_comp]=textread('suyeong9805.dat','');
% scatterplot ǥ��
subplot(2,2,1)
plot(u_comp,v_comp,'b.','MarkerSize',5)   % scatter plot 
hold on
plot([-50 50],[0 0],'k');plot([0 0],[-50 50],'k');% �߽�ǥ��
xlim([-100 100]);ylim([-100 100]);
axis equal
title('Scatter plot');xlabel('U comp (cm/sec)');ylabel('V comp (cm/sec)')

n=length(yy); %data�� ���
[rad_deg,spd]=cart2pol(u_comp,v_comp); %������ǥ�谪�� �� ��ǥ�� ������ (���� radian)

% Rose �Լ��� �̿��� histogram
subplot(2,2,2)
[tout, rout] = rose(rad_deg,12); % 30�� �������� ǥ��
% ä����� �����ϱ� ���Ͽ� �ʿ��� �۾�
polar(tout,rout);
set(gca,'nextplot','add');
[xout, yout] = pol2cart(tout,rout); % ����ǥ�踦 ������ǥ��� ��ȯ
set(gca,'nextplot','add');
fill(xout, yout, 'g');               % ä����� green
hline = findobj(gca,'Type','line');  % line ����
set(hline,'LineWidth',2.0,'Color','k')  % line�� �β��� 2.0, ���� blue��
title('Degrees Rose histogram')

% Histogram �׸���
% radian �� ������ degree
degree=rad_deg*180/pi;     % degree�� ǥ��(�׷��� ���� x�࿡�� �ݽð���� ������)
degree=450-degree;         % ���ϱ��ذ����� ����
for i=1:n,
    if degree(i)<0
        degree(i)=degree(i)+360;
    elseif degree(i)>360
        degree(i)=degree(i)-360;
    end
end
% current speed�� histogram
subplot(2,2,3)
hist(spd,10);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','c','EdgeColor','w')
xlabel('current speed(cm/sec)');  
ylabel('number');         
title('Speed histogram');

% ���ϱ��� ������ histogram
subplot(2,2,4)
hist(degree,12);   % 30�� �������� hisgogram �� 360�� 12�� ���� ����
h = findobj(gca,'Type','patch');
set(h,'FaceColor','c','EdgeColor','w')
xlabel('degree(\circ)'); ylabel('counts');
xlim([0 360])
title('Degrees histogram');