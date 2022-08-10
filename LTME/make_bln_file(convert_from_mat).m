close all; clear; clc; 
load coastline_80s.mat
temp=[lon lat];

index_lon=find(isnan(lon)==1);
for i = 1:length(index_lon)
    if i == 1
        number_point(i)=size(1:index_lon(i),2)-1;
    else
        number_point(i)=size(1+index_lon(i-1):index_lon(i),2)-1;
    end
end

number_point(end+1)=size(temp(index_lon(end)+1:end,1),1);

fid =fopen('1980s_coast_fix.bln','w');
for i = 1:length(number_point)
    if i==85
    temp1 = temp(1+index_lon(i-1):end,:);
    fprintf(fid,'%d\n',number_point(i));
    fprintf(fid,'%12.7f %12.7f\n',temp1');
%     clearvars temp1;
    else  
        if i == 1
            temp1 = temp(1:index_lon(i)-1,:);
            fprintf(fid,'%d\n',number_point(i));
            fprintf(fid,'%12.7f %12.7f\n',temp1');
            clearvars temp1;
        else
            temp1 = temp(1+index_lon(i-1):index_lon(i)-1,:);
            fprintf(fid,'%d\n',number_point(i));
            fprintf(fid,'%12.7f %12.7f\n',temp1');
            clearvars temp1;
        end
    end
end
fclose(fid);

