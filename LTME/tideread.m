clear all
fid=fopen('¿Œ√µ_DT_201310.txt')
data=textscan(fid,'%*8s%2f%3f%*1s%2f%*3s%5f%*2s%*2s%*2s%*2s%*2s%*2s%*4s%*6s%*3s%*5s%*7s%*2s','treatAsEmpty',{'-'},'headerlines',4)
day=data{1}
hour=data{2}
minute=data{3}
tide=data{4}

for i=1:25920
    if (tide(i)==0)
        tide(i)=tide(i-1)
    end
end
j=1
for i=1:1440
    if(mod(i,60)==0)
        exday(j)=day(i)
        exhour(j)=hour(i)
        exminute(j)=hour(i)
        extide(j)=tide(i)
        j=j+1
    end
end
figure;
plot(datenum(2013,10,exday,exhour,exminute,0),extide)
datetick('x','HH','keepticks')
title('2013.10.1.incheon')
xlabel('time(hour)')
ylabel('tide(m)')
fclose(fid)

ysum=0
xsum=0
Sxx=0
Sxy=0
n=1440/60
x=1:n
for i=1:n
    ysum=ysum+extide(i)
    xsum=xsum+i
end
yavr=ysum/n
xavr=xsum/n
for i=1:n
    Sxx=Sxx+(x(i)-xavr).^2
    Sxy=Sxy+(x(i)-xavr).*(extide(i)-yavr)
end
a_1=Sxy/Sxx
a_0=yavr-a_1*xavr
y=a_0+a_1*x
hold on 
%plot(datenum(2013,10,exday,exhour,exminute,0),y)