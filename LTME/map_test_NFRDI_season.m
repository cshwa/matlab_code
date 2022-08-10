clc;clear all;
close all;

%--- load 정선 관측 자료 짝수달 값  -------------------------------
load E:\11.사업\장기생태_3단계\1차년도\Data\국립수산과학원(KODC)_정선관측\GD.mat

%--- load 정선 자료 위치 정보 -------------------------------------
filename = fullfile('E:\11.사업\장기생태_3단계\1차년도\Data\국립수산과학원(KODC)_정선관측\NSOposition_silver.txt');
fileID = fopen(filename);
formatSpec = '%s';
N = 1;
C_text = textscan(fileID,formatSpec,N,'Delimiter','|');
C_data = textscan(fileID,'%s %f %f');
lat = C_data{2};
lon = C_data{3};

%--- 연안정지 관측 자료도 같이 나타내기 -----------------------------------
load E:\11.사업\장기생태_3단계\1차년도\draw_figure\data_from_yk_season\all_trend_season.txt
lon2 = all_trend_season(:,2);
lat2 = all_trend_season(:,1);

%--- jet color ---------------------------
% cmap = [0 0 255; 0 170 255; 255 255 0; 255 170 0; 255 85 0];
cmap = [122 16 228; 73 9 239; 10 92 253; 0 175 255; 139 249 255; 255 255 255; ...
        255 255 0; 255 192 0; 255 128 0; 255 65 0; 192 32 0; 128 0 0];
cmap = (cmap)/255;
%%
%--- m_map 으로 coastal line load -----------------------------------------
figure1 = figure('Position', [0,0, 1000, 700]);
ss = 1;
for j = [1 2 4 5]
subplot(2,2,ss)
data = GD(:,j);
%     set(figure1,'PaperUnits','inches','PaperPosition',[0 0 1 1.5]*6);
    m_proj('mercator','lon',[126.125 129.125],'lat',[33.125 35.125]);
    m_gshhs_i('color','k');
    m_gshhs_i('patch',[.9,.9,.9]);
    m_grid('box','fancy','tickdir','in','linewidth',1);
    hold on;
    % m_plot(lon, lat,'k.','markersize',15);
    ms = 18; % marker size
    c_lim = 0.06;
    c_min = -0.06;
        for i = 1:length(cmap)-1    
            c_int = [c_min:(c_lim-(c_min))/12:c_lim];
            [a b] = find(data > 0.06);
            m_plot(lon(a),lat(a),'MarkerEdgeColor',cmap(end,:), ...
                    'MarkerFaceColor',cmap(end,:),'MarkerSize',ms,'Marker','o',...
                    'MarkerEdgeColor','k','LineStyle','none');
            [a b] = find(data < c_min);
            m_plot(lon(a),lat(a),'MarkerEdgeColor',cmap(1,:), ...
                    'MarkerFaceColor',cmap(1,:),'MarkerSize',ms,'Marker','o',...
                    'MarkerEdgeColor','k','LineStyle','none');
            [a b] = find(data > c_int(i) & data <= c_int(i+1));
            m_plot(lon(a),lat(a),'Marker', 'o','LineStyle','none','MarkerEdgeColor','k', ...
                    'MarkerFaceColor',cmap(i,:),'MarkerSize',ms);
        hold on; 
        end
   

%--- 연안정지자료----------------------------------------------------------
    data2=all_trend_season(:,j+2);
    for i = 1:length(cmap)-1    
%         c_int = [c_min:(c_lim-(c_min))/5:c_lim];
        [a b] = find(data2 > 0.06);
        m_plot(lon2(a),lat2(a),'MarkerEdgeColor',cmap(end,:), ...
                'MarkerFaceColor',cmap(end,:),'MarkerSize',ms,'Marker','s','MarkerEdgeColor','k','LineStyle','none');
        [a b] = find(data2 < c_min);
        m_plot(lon2(a),lat2(a),'MarkerEdgeColor',cmap(1,:), ...
                'MarkerFaceColor',cmap(1,:),'MarkerSize',ms,'Marker','s','MarkerEdgeColor','k','LineStyle','none');
        [a b] = find(data2 > c_int(i) & data2 <= c_int(i+1));
        m_plot(lon2(a),lat2(a),'Marker', 's','MarkerEdgeColor',cmap(i,:), ...
                'MarkerFaceColor',cmap(i,:),'MarkerSize',ms,'MarkerEdgeColor','k','LineStyle','none');
    hold on; 
    end

     colormap(cmap);
    caxis([-0.06,0.06]);
    cbh = colorbar;
    c_int2 = [c_min:(c_lim-(c_min))/6:c_lim];
    set(cbh,'ytick',c_int2,'YTicklabel',c_int2); 
    ylabel(cbh,'C^{\circ}/ yr','fontsize',14);
    if ss == 1
        m_text(126.5,34.95,['겨울 (2월)'],'color','k','FontSize',17,'fontweight','bold')
    elseif ss == 2
         m_text(126.5,34.95,['봄 (4월)'],'color','k','FontSize',17,'fontweight','bold')
    elseif ss == 3
         m_text(126.5,34.95,['여름 (8월)'],'color','k','FontSize',17,'fontweight','bold')
    else
         m_text(126.5,34.95,['가을 (10월)'],'color','k','FontSize',17,'fontweight','bold')
    end
    ss=ss+1;
end
%==========================================================================




%--- save figure --------------------------------------
%     outpath = '.\figure\';
%     out=[outpath,'NFRDI_KOHA_season_trend3'];
%     set(gcf,'renderer','painter');
%     set(gcf, 'PaperUnits', 'inches');
%     x_width=10;
%     y_width=7;
%     set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
%     saveas(gcf,out,'tif');
%------------------------------------------------------


