clear all
fid=fopen('¸ñÆ÷_DT_201310.txt')
data=textscan(fid,'%*8s%2f%3f%*1s%2f%*3s%5f%6f%6f%*2s%*2s%*2s%*2s%*4s%*6s%*3s%*5s%*7s%*2s','treatAsEmpty',{'-'},'headerlines',4)
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
title('2013.10.1.mokpo')
xlabel('time(hour)')
ylabel('tide(m)')
fclose(fid)