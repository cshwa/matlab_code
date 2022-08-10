clear all
maxnum=1000;
u1(1:maxnum)=0;
u2(1:maxnum)=0;
u3(1:maxnum)=0;
u4(1:maxnum)=0;
u5(1:maxnum)=0;
u6(1:maxnum)=0;
u7(1:maxnum)=0;
u8(1:maxnum)=0;
u9(1:maxnum)=0;

%%% direct
A(1:9,1:9)=0;
B(1:9)=0;
B(1)=-1;
B(2)=-8;
B(3)=-75;
B(6)=-32;
B(9)=-16

for i=1:9
    A(i,i)=-4;
end
for i=1:8
    A(i,i+1)=1;
    A(i+1,i)=1;
end
for i=1:6
    A(i+3,i)=1;
    A(i,i+3)=1;
end
A(3,4)=0;
A(4,3)=0;
A(7,6)=0;
A(6,7)=0
invA=inv(A)
x=invA*transpose(B)

directx=x;

%%% Jacobi
n=1
x(1:9)=0;
tempx(1:9)=x;
x(1)=tempx(2)/4 + tempx(4)/4 + 1/4;
x(2)=tempx(1)/4 + tempx(3)/4 + tempx(5)/4 + 2*2*2/4;
x(3)=tempx(2)/4 + 16*3/4 + tempx(6)/4 + 3*3*3/4;
x(4)=0 + tempx(5)/4 + tempx(7)/4 + tempx(1)/4;
x(5)=tempx(4)/4 + tempx(6)/4 + tempx(8)/4 + tempx(2)/4;
x(6)=tempx(5)/4 + 16*2/4 + tempx(9)/4 + tempx(3)/4;
x(7)=0/4 + tempx(8)/4 + 0/4 + tempx(4)/4;
x(8)=tempx(7)/4 + tempx(9)/4 + 0/4 + tempx(5)/4;
x(9)=tempx(8)/4 + 16*1/4 + 0/4 + tempx(6)/4;
while sum(x)-sum(tempx)>0.00001
    n=n+1
    tempx=x;
    x(1)=tempx(2)/4 + tempx(4)/4 + 1/4;
    x(2)=tempx(1)/4 + tempx(3)/4 + tempx(5)/4 + 2*2*2/4;
    x(3)=tempx(2)/4 + 16*3/4 + tempx(6)/4 + 3*3*3/4;
    x(4)=0 + tempx(5)/4 + tempx(7)/4 + tempx(1)/4;
    x(5)=tempx(4)/4 + tempx(6)/4 + tempx(8)/4 + tempx(2)/4;
    x(6)=tempx(5)/4 + 16*2/4 + tempx(9)/4 + tempx(3)/4;
    x(7)=0/4 + tempx(8)/4 + 0/4 + tempx(4)/4;
    x(8)=tempx(7)/4 + tempx(9)/4 + 0/4 + tempx(5)/4;
    x(9)=tempx(8)/4 + 16*1/4 + 0/4 + tempx(6)/4;
end
jacobin=n;
jacobix=x;

%%% gauss-seidel
n=1
x(1:9)=0;
tempx(1:9)=x;
x(1)=x(2)/4 + x(4)/4 + 1/4;
x(2)=x(1)/4 + x(3)/4 + x(5)/4 + 2*2*2/4;
x(3)=x(2)/4 + 16*3/4 + x(6)/4 + 3*3*3/4;
x(4)=0 + x(5)/4 + x(7)/4 + x(1)/4;
x(5)=x(4)/4 + x(6)/4 + x(8)/4 + x(2)/4;
x(6)=x(5)/4 + 16*2/4 + x(9)/4 + x(3)/4;
x(7)=0/4 + x(8)/4 + 0/4 + x(4)/4;
x(8)=x(7)/4 + x(9)/4 + 0/4 + x(5)/4;
x(9)=x(8)/4 + 16*1/4 + 0/4 + x(6)/4;
while sum(x)-sum(tempx)>0.00001
    n=n+1
    tempx=x;
    x(1)=x(2)/4 + x(4)/4 + 1/4;
    x(2)=x(1)/4 + x(3)/4 + x(5)/4 + 2*2*2/4;
    x(3)=x(2)/4 + 16*3/4 + x(6)/4 + 3*3*3/4;
    x(4)=0 + x(5)/4 + x(7)/4 + x(1)/4;
    x(5)=x(4)/4 + x(6)/4 + x(8)/4 + x(2)/4;
    x(6)=x(5)/4 + 16*2/4 + x(9)/4 + x(3)/4;
    x(7)=0/4 + x(8)/4 + 0/4 + x(4)/4;
    x(8)=x(7)/4 + x(9)/4 + 0/4 + x(5)/4;
    x(9)=x(8)/4 + 16*1/4 + 0/4 + x(6)/4;
end
gsn=n;
gsx=x;

%%% S.O.R
n=1;
omega=1.18
x(1:9)=0;
tempx(1:9)=x;
x(1)=omega*( x(2)/4 + x(4)/4 + 1/4 ) + (1-omega)* tempx(1);
x(2)=omega*( x(1)/4 + x(3)/4 + x(5)/4 + 2*2*2/4) + (1-omega)* tempx(2);
x(3)=omega*( x(2)/4 + 16*3/4 + x(6)/4 + 3*3*3/4) + (1-omega)* tempx(3);
x(4)=omega*( 0 + x(5)/4 + x(7)/4 + x(1)/4) + (1-omega)* tempx(4);
x(5)=omega*( x(4)/4 + x(6)/4 + x(8)/4 + x(2)/4) + (1-omega)* tempx(5);
x(6)=omega*( x(5)/4 + 16*2/4 + x(9)/4 + x(3)/4) + (1-omega)* tempx(6);
x(7)=omega*( 0/4 + x(8)/4 + 0/4 + x(4)/4) + (1-omega)* tempx(7);
x(8)=omega*( x(7)/4 + x(9)/4 + 0/4 + x(5)/4) + (1-omega)* tempx(8);
x(9)=omega*( x(8)/4 + 16*1/4 + 0/4 + x(6)/4) + (1-omega)* tempx(9);
while sum(x)-sum(tempx)>0.00001
    n=n+1
    tempx=x;
    x(1)=omega*( x(2)/4 + x(4)/4 + 1/4 ) + (1-omega)* tempx(1);
    x(2)=omega*( x(1)/4 + x(3)/4 + x(5)/4 + 2*2*2/4) + (1-omega)* tempx(2);
    x(3)=omega*( x(2)/4 + 16*3/4 + x(6)/4 + 3*3*3/4) + (1-omega)* tempx(3);
    x(4)=omega*( 0 + x(5)/4 + x(7)/4 + x(1)/4) + (1-omega)* tempx(4);
    x(5)=omega*( x(4)/4 + x(6)/4 + x(8)/4 + x(2)/4) + (1-omega)* tempx(5);
    x(6)=omega*( x(5)/4 + 16*2/4 + x(9)/4 + x(3)/4) + (1-omega)* tempx(6);
    x(7)=omega*( 0/4 + x(8)/4 + 0/4 + x(4)/4) + (1-omega)* tempx(7);
    x(8)=omega*( x(7)/4 + x(9)/4 + 0/4 + x(5)/4) + (1-omega)* tempx(8);
    x(9)=omega*( x(8)/4 + 16*1/4 + 0/4 + x(6)/4) + (1-omega)* tempx(9);
end
sorn=n;
sorx=x;
x
