x=1:100
sum(1)=1
i=2
while i<101;
    sum(i)=sum(i-1)+i
    i=i+1
end
plot(x,sum)

!!!!

figure;
x=1:100
sum2(1)=1
for i=2:100
    sum2(i)=sum2(i-1)+i
end
plot(x,sum2)