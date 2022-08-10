clear all
k=10000
L=10000000.0;
dx=100000.0;
dt=500000;
for i=1:101
    x(i)=(i-1)*dx;
    if (i<51)
        u(i,1)=x(i)/L;
    elseif(i>=51)
        u(i,1)=(L-x(i))/L;
    end
end

for t=1:1000
    u(1,t)=0;
    u(101,t)=0;
    for i=2:100
        u(i,t+1)=u(i,t)+k*(u(i+1,t)-2*u(i,t)+u(i-1,t))/dx^2*dt;
    end
end

plot(u(:,1))
hold on
plot(u(:,1001))