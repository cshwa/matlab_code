
close all
clear all
clc

days = [31 28 31 30 31 30 31 31 30 31 30 31];
days2 = [31 29 31 30 31 30 31 31 30 31 30 31];

 fid = fopen('time_ind.dat','w');
dd = 1;
for i = 2:length(days)+1
   
    fprintf(fid,'%10d\n',dd);
    dd = days(i-1)+dd;
    
    
end
fclose(fid);