% 조류data를 간단히 표현하는 방법으로 Scatter plot, histogram, rose등을 활용
% m-file: tidalcurrent_exam1.m
% 주의: 아직 rose 함수에서 x기준각도를 정북기준각도로 바꾸지 못함

clc;clear all;close all
% data load
[yy,mm,dd,hr,min,temp,speed,n_dir,u_comp,v_comp]=textread('suyeong9805.dat','');
% scatterplot 표현
subplot(2,2,1)
plot(u_comp,v_comp,'b.','MarkerSize',5)   % scatter plot 
hold on
plot([-50 50],[0 0],'k');plot([0 0],[-50 50],'k');% 중심표현
xlim([-100 100]);ylim([-100 100]);
axis equal
title('Scatter plot');xlabel('U comp (cm/sec)');ylabel('V comp (cm/sec)')

n=length(yy); %data의 행수
[rad_deg,spd]=cart2pol(u_comp,v_comp); %직각좌표계값을 극 좌표계 값으로 (각은 radian)

% Rose 함수를 이용한 histogram
subplot(2,2,2)
[tout, rout] = rose(rad_deg,12); % 30도 간격으로 표현
% 채움색을 설정하기 위하여 필요한 작업
polar(tout,rout);
set(gca,'nextplot','add');
[xout, yout] = pol2cart(tout,rout); % 극좌표계를 직각좌표계로 변환
set(gca,'nextplot','add');
fill(xout, yout, 'g');               % 채움색을 green
hline = findobj(gca,'Type','line');  % line 선택
set(hline,'LineWidth',2.0,'Color','k')  % line의 두께를 2.0, 색을 blue로
title('Degrees Rose histogram')

% Histogram 그리기
% radian 을 각도로 degree
degree=rad_deg*180/pi;     % degree로 표현(그러나 아직 x축에서 반시계방향 각도임)
degree=450-degree;         % 정북기준각도로 변경
for i=1:n,
    if degree(i)<0
        degree(i)=degree(i)+360;
    elseif degree(i)>360
        degree(i)=degree(i)-360;
    end
end
% current speed의 histogram
subplot(2,2,3)
hist(spd,10);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','c','EdgeColor','w')
xlabel('current speed(cm/sec)');  
ylabel('number');         
title('Speed histogram');

% 정북기준 유향의 histogram
subplot(2,2,4)
hist(degree,12);   % 30도 간격으로 hisgogram 즉 360을 12로 나눈 각도
h = findobj(gca,'Type','patch');
set(h,'FaceColor','c','EdgeColor','w')
xlabel('degree(\circ)'); ylabel('counts');
xlim([0 360])
title('Degrees histogram');