a0=31.4
a1=0.1
a2=0.1
A(1:3, 1:3)=0
for i=1:n
    A(1,1)=A(1,1)+1;
    A(1,2)=A(1,2)-(exp(a2*x(i)))
    A(1,3)=A(1,3)-(y(i)/a2)
    
    A(2,1)=A(2,1)+(