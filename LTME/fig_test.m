

close all
clear all
clc

foot = '_u.txt';

for i=1:6:2201
    num = num2char(i,4);
    file = ([num, foot]);
    data = load(file);
    dep = data(:,2);
    u = data(:,3);
    hold on
    if(i>=1 & i<=300)
        plot(u, dep,'r');
    elseif(i>=301 & i<=600)
        plot(u, dep,'m');
    elseif(i>=601 & i<=840)
        plot(u, dep,'b');
    elseif(i>=1021 & i<=1200)
        plot(u, dep,'g');
    elseif(i>=1201 & i<=1500)
        plot(u, dep,'y');
    elseif(i>=1501 & i<=1800)
        plot(u, dep,'k');
    elseif(i>=1801 & i<=2200)
        plot(u, dep,'c');
    end
end
xlim = [-2,2];
set(gca,'XLim',xlim(:))
set(gca,'fontsize',14)
xlabel('Velocity (m/s)','fontsize',18)
ylabel('Depth from bottom (m)','fontsize',18)