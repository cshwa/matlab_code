clear all
c=20.0
L=2000000.0;
dx=200000.0;
dt=3600.;
a=dt*c/dx
for i=1:11
    x(i)=(i-1)*dx;
    u1(i,1)=c + 10.0*sin(2*pi/L*x(i));
    u2(i,1)=c + 10.0*sin(2*pi/L*x(i));
end

for t=1:10000
    for i=1:11
        if (i==1)
            u1(1,t+1)=u1(1,t)-a*(u1(2,t)-u1(10,t))/2.0;
        elseif(i==11)
            u1(11,t+1)=u1(11,t)-a*(u1(2,t)-u1(10,t))/2.0;
        else
            u1(i,t+1)=u1(i,t)-a*(u1(i+1,t)-u1(i-1,t))/2.0;
        end
    end
end
t=1
for i=1:11
    if (i==1)
        u2(1,t+1)=u2(1,t)-a*(u2(2,t)-2*u2(1,t)+u2(10,t))/2.0;
    elseif(i==11)
        u2(11,t+1)=u2(11,t)-a*(u2(2,t)-2*u2(11,t)+u2(10,t))/2.0;
    else
        u2(i,t+1)=u2(i,t)-a*(u2(i+1,t)-2*u2(i,t)+u2(i-1,t))/2.0;
    end
end

for t=2:10000
    for i=1:11
        if (i==1)
            u2(1,t+1)=u2(1,t-1)-a*(u2(2,t)-u2(10,t));
        elseif(i==11)
            u2(11,t+1)=u2(1,t-1)-a*(u2(2,t)-u2(10,t));
        else
            u2(i,t+1)=u2(i,t-1)-a*(u2(i+1,t)-u2(i-1,t));
        end
    end
end
% 
% plot (u1(:,10001))
% hold on
% plot (u2(:,10001))
% for i=9000:10000
%     plot(u1(:,i))
%     pause(0.1)
% end
plot(u1(:,1))
hold on
plot(u1(:,5))
figure;
plot(u1(:,5001))