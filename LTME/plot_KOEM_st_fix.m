close all; clear; clc;   % -v3
lon=ncread('grid_sumjin_v1970_fix_3m.nc','lon_rho');
lat=ncread('grid_sumjin_v1970_fix_3m.nc','lat_rho');
mask=ncread('grid_sumjin_v1970_fix_3m.nc','mask_rho');
h=ncread('grid_sumjin_v1970_fix_3m.nc','h');
load KOEM_st_info_(decimal_deg).mat
lat_koem = [lat_1; lat_2; lat_3;];
lon_koem = [lon_1; lon_2; lon_3;];

% station name tag
for i = 1:3
name_tag_1{i} = ['여수항' num2str(i,'%02d')] 
end

name_tag_2{1} = ['광양항' num2str(1,'%02d')] 

name_tag_3{1} = ['삼천포' num2str(1,'%02d')] 

for i = 1:5
name_tag_4{i} = ['가막만' num2str(i,'%02d')] 
end

for i = 1:25
name_tag_5{i} = ['섬진강' num2str(i,'%02d')] 
end

for i = 1:2
name_tag_6{i} = ['진주만' num2str(i,'%02d')] 
end

for i = 1:28
name_tag_7{i} = ['연안' num2str(i,'%02d')] 
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

name_tag_97 = zeros(65,1); %mask for 1997 starting

% plot masking
figure; hold on; pcolor(lon,lat,(mask./mask)); shading flat;
set(gca,'fontsize',15,'fontweight','bold'); 
grid on; plot(lon_koem,lat_koem,'.','color','k')
for i=1:length(name_tag)
tx=text(lon_koem(i),lat_koem(i), char(name_tag{i}),'color',[1 1 1]);
% set(tx,'Rotation',25)
% set(tx,'FontSize',9{
end
% 연안 3, 6, 13:28 has to be remove (18)  [it's located at out of domain] 

name_tag_n=load('KOEM_name_tag.mat','name_tag');
name_tag_n = name_tag_n.name_tag;
for i = 1:length(name_tag); tag_mask{i,1}=name_tag_n{i,2}; end

non_97=find(cellfun('isempty',tag_mask))  % not 97'
yeah_97=find(~cellfun('isempty',tag_mask))  % 97' is 40 st.
return
for i = 1:length(non_97); tag_mask{non_97(i)} = 0; end  %miss cells replace to be 0


%plot only in the domain
figure; hold on; %pcolor(lon,lat,(mask./mask)); shading flat;
set(gca,'fontsize',15,'fontweight','bold'); 
grid on; 
for i=1:length(name_tag)
    if tag_mask{i} == 1
        text(lon_koem(i),lat_koem(i), char(name_tag{i}),'color',[1 1 1]); %plot only 97
        plot(lon_koem(i),lat_koem(i),'.','color','r'); %plot only 97
    end
end
ylim([34.723 35])
xlim([min(min(lon)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','satellite')  % overlay google map
ylim([34.723 35])
xlim([min(min(lon)) 127.9])

% full plot
figure; hold on; %pcolor(lon,lat,(mask./mask)); shading flat;
set(gca,'fontsize',15,'fontweight','bold'); 
grid on; 
for i=1:length(name_tag)
    if tag_mask{i} == 1
        text(lon_koem(i),lat_koem(i), char(name_tag{i}),'color',[1 1 1]); %plot only 97
        plot(lon_koem(i),lat_koem(i),'.','color','r'); %plot only 97
    end
end
ylim([34.5 35])
xlim([min(min(lon)) 128.7])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','satellite')  % overlay google map
% ylim([34.5 35])
% xlim([min(min(lon)) 128.7])
xlabel('LON'); ylabel('LAT')
% 
% plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map


%%
%% plot point trend
load('kodc_just_mean_yearly_linear_trend.mat','coeff_*');

%temp : 0 ~ 0.12 
max(coeff_temp(:,1))
min(coeff_temp(:,1))

%salt : 0 ~ 0.9
max(coeff_salt(:,1))
min(coeff_salt(:,1))

%no3 : -4.5~0
max(coeff_no3(:,1))
min(coeff_no3(:,1))


%--- color map ------------------------------------------------------------
% cmap = [255 85 0; 255 170 0;255 255 0;170 255 85;85 255 170; 0 255 255; ...
%         0 170 255; 0 85 255;0 0 255];
% cmap = flipud(cmap)/255;

%% no3
cmap = [255 85 0; 255 170 0;255 255 0; ...
        0 170 255; 0 85 255;0 0 255];
cmap = flipud(cmap)/255;
figure; hold on; %pcolor(lon,lat,(mask./mask)); shading flat;
set(gca,'fontsize',15,'fontweight','bold'); 
grid on; 
ms = 100; % marker size
clearvars data data_pre lon_plt lat_plt
data_pre = coeff_no3(:,1).*(length(1997:2018)-1); % rate to absolute value 
data = coeff_no3(:,1).*(length(1997:2018)-1); 
lon_plt=lon_koem; lat_plt=lat_koem;
lon_plt(find( data_pre == 0)) = []; % remove lack of timeseries point
lat_plt(find( data_pre == 0)) = [];
data(find( data_pre == 0)) = [];
max(data)
min(data)
c_lim = -100;
for i = 1:length(cmap)-1
%     c_int = [linspace(-100,0,9)];
    c_int = -100:20:0;
    [a b] = find(data > 0);
%     m_plot(lon_koem(a),lat_koem(a),'MarkerEdgeColor',cmap(9,:), ...
%             'MarkerFaceColor',cmap(9,:),'MarkerSize',ms,'Marker','o');
    scatter(lon_plt(a),lat_plt(a),ms,'MarkerEdgeColor',cmap(end,:), ...
            'MarkerFaceColor',cmap(end,:),'Marker','o');
    [a b] = find(data > c_int(i) & data <= c_int(i+1));
%     m_plot(lon_koem(a),lat_koem(a),'Marker', 'o','MarkerEdgeColor',cmap(i,:), ...
%             'MarkerFaceColor',cmap(i,:),'MarkerSize',ms);
    scatter(lon_plt(a),lat_plt(a),ms,'MarkerEdgeColor',cmap(i,:), ...
            'MarkerFaceColor',cmap(i,:),'Marker', 'o');
hold on; 
end
cbh = colorbar;
set(cbh,'ytick',[0:1/5:1],'YTicklabel',c_int);
set(get(cbh,'Title'),'String','ug/L','fontsize',12,'fontweight','bold')
%--- 제목 수정하기 -------------------------------------
title('KOEM 표층NO3 감소량 (1997-2018)','fontsize',18);

ylim([34.5 35])
xlim([min(min(lon)) 128.7])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)  % overlay google map
ylim([34.5 35])
xlim([min(min(lon)) 128.7])
xlabel('경도'); ylabel('위도')
% for i = 1:length(c_int)
%     c_s_int{i}=num2str(c_int(i),'%.2f'); 
% end;
set(cbh,'ytick',[0:1/5:1],'YTicklabel',c_int,'fontweight','bold','fontsize',14);

%% salt
cmap = [255 85 0; 255 170 0;255 255 0;170 255 85; ...
        0 170 255; 0 85 255;0 0 255]
cmap = flipud(cmap)/255;

figure; hold on; %pcolor(lon,lat,(mask./mask)); shading flat;
set(gca,'fontsize',15,'fontweight','bold'); 
grid on; 
ms = 100; % marker size
clearvars data data_pre lon_plt lat_plt
data_pre = coeff_salt(:,1).*(length(1997:2018)-1); % rate to absolute value
data = coeff_salt(:,1).*(length(1997:2018)-1); % rate to absolute value 
lon_plt=lon_koem; lat_plt=lat_koem;
lon_plt(find( data_pre == 0)) = []; % remove lack of timeseries point
lat_plt(find( data_pre == 0)) = [];
data(find( data_pre == 0)) = [];
max(data)
min(data)
c_lim = 1.75;
for i = 1:length(cmap)
%     c_int = [linspace(0,1.75,9)];
    c_int = 0:0.25:1.75;
    [a b] = find(data > c_lim);
    scatter(lon_plt(a),lat_plt(a),ms,'MarkerEdgeColor',cmap(end,:), ...
            'MarkerFaceColor',cmap(end,:),'Marker','o');
    [a b] = find(data > c_int(i) & data <= c_int(i+1));
    scatter(lon_plt(a),lat_plt(a),ms,'MarkerEdgeColor',cmap(i,:), ...
            'MarkerFaceColor',cmap(i,:),'Marker', 'o');
hold on; 
end
cbh = colorbar;
set(cbh,'ytick',0:1/7:1,'YTicklabel',c_int);
set(get(cbh,'Title'),'String','o/oo','fontsize',12,'fontweight','bold')
%--- 제목 수정하기 -------------------------------------
title('KOEM 표층염분 증가량 (1997-2018)','fontsize',18);

ylim([34.5 35])
xlim([min(min(lon)) 128.7])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.5 35])
xlim([min(min(lon)) 128.7])
xlabel('경도'); ylabel('위도')
clearvars c_s_int
for i = 1:length(c_int)
    if i==1
        c_s_int{i}=num2str(c_int(i)); 
    else
    c_s_int{i}=num2str(c_int(i),'%.2f'); 
    end
end;
set(cbh,'ytick',0:1/7:1,'YTicklabel',c_s_int,'fontweight','bold','fontsize',14);


%% temp
cmap = [255 85 0; 255 170 0;255 255 0; ...
        0 85 255;0 0 255];
cmap = flipud(cmap)/255;

figure; hold on; %pcolor(lon,lat,(mask./mask)); shading flat;
set(gca,'fontsize',15,'fontweight','bold'); 
grid on; 
ms = 100; % marker size
clearvars data data_pre lon_plt lat_plt
data_pre = coeff_temp(:,1).*(length(1997:2018)-1); % rate to absolute value
data = coeff_temp(:,1).*(length(1997:2018)-1); % rate to absolute value
lon_plt=lon_koem; lat_plt=lat_koem;
lon_plt(find( data_pre == 0)) = []; % remove lack of timeseries point
lat_plt(find( data_pre == 0)) = [];
data(find( data_pre == 0)) = [];
max(data)
min(data)
c_lim = 2.5;
for i = 1:length(cmap)
%     c_int = [linspace(0,2.5,9)];
    c_int = 0:0.5:2.5;
    [a b] = find(data > c_lim);
    scatter(lon_plt(a),lat_plt(a),ms,'MarkerEdgeColor',cmap(end,:), ...
            'MarkerFaceColor',cmap(end,:),'Marker','o');
    [a b] = find(data > c_int(i) & data <= c_int(i+1));
    scatter(lon_plt(a),lat_plt(a),ms,'MarkerEdgeColor',cmap(i,:), ...
            'MarkerFaceColor',cmap(i,:),'Marker', 'o');
hold on; 
end
cbh = colorbar;
% set(cbh,'ytick',[0:1/8:1],'YTicklabel',c_int);
set(cbh,'ytick',0:1/5:1,'YTicklabel',c_int);
set(get(cbh,'Title'),'String','^oC','fontsize',12,'fontweight','bold')
%--- 제목 수정하기 -------------------------------------
title('KOEM 표층수온 증가량 (1997-2018)','fontsize',18);

ylim([34.5 35])
xlim([min(min(lon)) 128.7])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.5 35])
xlim([min(min(lon)) 128.7])
xlabel('경도'); ylabel('위도')
% set(cbh,'ytick',[0:1/8:1],'YTicklabel',c_int);
for i = 1:length(c_int)
    if i==1
        c_s_int{i}=num2str(c_int(i)); 
    else
    c_s_int{i}=num2str(c_int(i),'%.1f'); 
    end
end;
set(cbh,'ytick',0:1/5:1,'YTicklabel',c_s_int,'fontweight','bold','fontsize',14);

% save('KOEM_name_tag.mat','name_tag');