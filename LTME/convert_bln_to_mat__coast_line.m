% 1980s
close all; clear all; clc;
filename='1980s_coast_fix_force2close_final.txt';
fileID = fopen(filename);
C = textscan(fileID,'%f %f');
fclose(fileID);

lon=C{1};lat=C{2};
cut_index=find(isnan(lat)==1);

% lon(cut_index)=[];lat(cut_index)=[];
lon(cut_index)=NaN;
save('coast_line_1980s_final.mat','lon','lat');

for i=1:length(cut_index)-1
if i ==1
    temlon=lon(cut_index(i):cut_index(i+1)); 
    temlon(isnan(temlon))=[];
    island_lon{i}=temlon;
    temlat=lat(cut_index(i):cut_index(i+1));
    temlat(isnan(temlat))=[];
    island_lat{i}=temlat;
else
    temlon=lon(cut_index(i)+1:cut_index(i+1));
    temlon(isnan(temlon))=[];
    island_lon{i}=temlon;
    temlat=lat(cut_index(i)+1:cut_index(i+1));
    temlat(isnan(temlat))=[];
    island_lat{i}=temlat;  
end
temlon=[];temlat=[];
end
save('coast_line_1980s_island.mat','island_lon','island_lat');

%% present
close all; clear all; clc;
filename='fine_coast_KOR.txt';
fileID = fopen(filename);
C = textscan(fileID,'%f %f');
fclose(fileID);

p_lon=C{1};p_lat=C{2};
cut_index=find(isnan(p_lat)==1);

% lon(cut_index)=[];lat(cut_index)=[];
p_lon(cut_index)=NaN;
save('coast_line_present.mat','p_lon','p_lat');

for i=1:length(cut_index)-1
if i ==1
    temp_lon=p_lon(cut_index(i):cut_index(i+1)); 
    temp_lon(isnan(temp_lon))=[];
    p_island_lon{i}=temp_lon;
    temp_lat=p_lat(cut_index(i):cut_index(i+1));
    temp_lat(isnan(temp_lat))=[];
    p_island_lat{i}=temp_lat;
else
    temp_lon=p_lon(cut_index(i)+1:cut_index(i+1));
    temp_lon(isnan(temp_lon))=[];
    p_island_lon{i}=temp_lon;
    temp_lat=p_lat(cut_index(i)+1:cut_index(i+1));
    temp_lat(isnan(temp_lat))=[];
    p_island_lat{i}=temp_lat;  
end
temp_lon=[];temp_lat=[];
end
save('coast_line_present_island.mat','p_island_lon','p_island_lat');

