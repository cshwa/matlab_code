close all; clear; clc; 

load KOEM_st_info_(decimal_deg).mat
lat_koem = [lat_1; lat_2; lat_3;];
lon_koem = [lon_1; lon_2; lon_3;];

% station name tag
for i = 1:3
name_tag_1{i} = ['여수신항 H' num2str(i)] 
end

name_tag_2{1} = ['광양항 H' num2str(1)] 

name_tag_3{1} = ['삼천포항 H' num2str(1)] 

for i = 1:5
name_tag_4{i} = ['가막만 ' num2str(i,'%02d')] 
end

for i = 1:25
name_tag_5{i} = ['섬진강하구 ' num2str(i,'%02d')] 
end

for i = 1:2
name_tag_6{i} = ['진주만 ' num2str(i,'%02d')] 
end

for i = 1:28
name_tag_7{i} = ['대한해협연안 ' num2str(i,'%02d')] 
end

name_tag = name_tag_1'; 
name_tag{end+1:end+length(name_tag_2)} = name_tag_2; 
name_tag{end+1:end+length(name_tag_3)} = name_tag_3; 
size_tag = length(name_tag);

for i = 1:length(name_tag_4)
name_tag{size_tag+i} = name_tag_4{i}; 
end
size_tag = length(name_tag);

for i = 1:length(name_tag_5)
name_tag{size_tag+i} = name_tag_5{i}; 
end
size_tag = length(name_tag);

for i = 1:length(name_tag_6)
name_tag{size_tag+i} = name_tag_6{i}; 
end
size_tag = length(name_tag);

for i = 1:length(name_tag_7)
name_tag{size_tag+i} = name_tag_7{i};
end

% [raw_h1 txt_h1]=xlsread('해양환경측정망(항만측정망).xls','sheet1');  
% [raw_h1 txt_h1]=xlsread('해양환경측정망(항만측정망).xls','sheet1');


% name_tag_97 = zeros(65,1); %mask for 1997 starting

% plot masking
figure; hold on; pcolor(lon,lat,(mask./mask)); shading flat;
set(gca,'fontsize',15,'fontweight','bold'); 
grid on; plot(lon_koem,lat_koem,'.','color','k')
for i=1:length(name_tag)
tx=text(lon_koem(i),lat_koem(i), char(name_tag{i}));
% set(tx,'Rotation',25)
% set(tx,'FontSize',9)
end
% 연안 3, 6, 13:28 has to be remove (18)  [it's located at out of domain] 

name_tag_n=load('KOEM_name_tag.mat','name_tag');
name_tag_n = name_tag_n.name_tag;
for i = 1:length(name_tag); tag_mask{i,1}=name_tag_n{i,2}; end

non_97=find(cellfun('isempty',tag_mask))  % not 97'
%yeah_97=find(~cellfun('isempty',tag_mask))  % 97' is 40 st.

for i = 1:length(non_97); tag_mask{non_97(i)} = 0; end  %miss cells replace to be 0
