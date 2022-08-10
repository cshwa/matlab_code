clear all
pi=3.141592
c=0.02
A=0.01
dx=200
L=10*dx
k=2*pi/L
x=0:dx:10*dx
u0=c+A*sin(k*x)

laxu=u0;
lelu=u0;
leapu=u0;
dt=100;
nn=1:11
templelu=lelu
templaxu=laxu
for i=1:1000
    for n=1:11
        if (n==1)
            laxu(n)=(1/2-c*dt/(2*dx))*templaxu(n+1) + (1/2+c*dt/(2*dx))*templaxu(n+10);
            lelu(n)=templelu(n)-(c*dt/(dx))*(templelu(n)-templelu(n+9));
        elseif(n==11)
            laxu(n)=(1/2-c*dt/(2*dx))*templaxu(n-10)+(1/2+c*dt/(2*dx))*templaxu(n-1);
            lelu(n)=templelu(n)-(c*dt/(dx))*(templelu(n)-templelu(n-1));
        else
            laxu(n)=(1/2-c*dt/(2*dx))*templaxu(n+1)+(1/2+c*dt/(2*dx))*templaxu(n-1);
            lelu(n)=templelu(n)-(c*dt/(dx))*(templelu(n)-templelu(n-1));
        end
    end
    templelu=lelu;
    templaxu=laxu;

end

bleapu=leapu;
for n=1:11       
        if (n==1)
            leapu(n)=(1/2-c*dt/(2*dx))*leapu(n+1) + (1/2+c*dt/(2*dx))*leapu(n+10);
        elseif(n==11)
            leapu(n)=(1/2-c*dt/(2*dx))*leapu(n-10)+(1/2+c*dt/(2*dx))*leapu(n-1);
        else
            leapu(n)=(1/2-c*dt/(2*dx))*leapu(n+1)+(1/2+c*dt/(2*dx))*leapu(n-1);
        end
end
templeapu=leapu;

for i=1:1000
    for n=1:11
        if (n==1)
            leapu(n)=bleapu(n)-dt/dx*c*(templeapu(n+1)-templeapu(n+10));
        elseif(n==11)
            leapu(n)=bleapu(n)-dt/dx*c*(templeapu(n-10)-templeapu(n-1));
        else
            leapu(n)=bleapu(n)-dt/dx*c*(templeapu(n+1)-templeapu(n-1));
        end
    end
    bleapu=templeapu;
    templeapu=leapu;
%     plot(nn,leapu)
%     pause(0.01)
end

plot(nn,laxu)
hold on
plot(nn,lelu)
plot(nn,leapu)