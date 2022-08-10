sta = load('kodc_sta.dat');
data  = load('kodc_hydro.dat');
[len_t,error]=size(data);
fid = fopen('kodc_corrected.dat','w');
for t=1:len_t
    if(data(t,6)==0)
        index = find(sta(:,1)==data(t,4) & sta(:,2)==data(t,5));
        fprintf(fid,'%5d%5d%5d%5d%5d%9.4f%9.4f%10.4f%6.1f%7.2f%7.2f\n',data(t,1:5),sta(index(1),3),sta(index(1),4),data(t,8:11));
     else
        fprintf(fid,'%5d%5d%5d%5d%5d%9.4f%9.4f%10.4f%6.1f%7.2f%7.2f\n',data(t,:));
    end
end

fclose(fid);