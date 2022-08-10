x=1:100
y(1)=1
for i=1:99
    y(i+1)=y(i)+i+1
end
plot(x,y)