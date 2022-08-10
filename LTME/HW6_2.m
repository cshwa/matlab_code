clear all
N=8:100;
pi=3.14;
o=2.*pi./N;
a=sin(o)./o;
b=(4./3.*sin(o)-1./6.*sin(o.*2))./o;
c=(2./3.*sin(o.*2)-1./12.*sin(4.*o))./o;
d=(1./3.*sin(o.*4)-1./24.*sin(8.*o))./o;

plot(N,a)
hold on
plot(N,b)
plot(N,c)
plot(N,d)