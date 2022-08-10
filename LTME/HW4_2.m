clear all
x(1:1000,1)=0;
x(1:1000,2)=0;
x(1:1000,3)=0;


for i=1:1000
    x(i+1,1)=1-2*x(i,2)+2*x(i,3);
    x(i+1,2)=3-x(i+1,1)-x(i,3);
    x(i+1,3)=5-2*x(i+1,1)-2*x(i+1,2);
end

gsx=x;

x(1:1000,1)=0;
x(1:1000,2)=0;
x(1:1000,3)=0;

for i=1:1000
    tempx=x;
    x(i+1,1)=1-2*tempx(i,2)+2*tempx(i,3);
    x(i+1,2)=3-tempx(i+1,1)-tempx(i,3);
    x(i+1,3)=5-2*tempx(i+1,1)-2*tempx(i+1,2);
end

jacobix=x;