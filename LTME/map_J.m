function map_J(casename)

[lon_lim, lat_lim] = domain_J(casename);
lim = [lon_lim lat_lim];
fc = [.95 .95 .95 ];

%figure; hold on;
hold on;
set(gca,'ydir','nor');
m_proj('miller','lon',[lim(1) lim(2)],'lat',[lim(3) lim(4)]);
m_gshhs_i('line', 'Color','k') % plot coastline at intermediate resolution, as black line
m_gshhs_i('patch',fc )
m_grid('XTick', 6, 'YTick', 6, ...
    'LineWidth', 2, 'LineStyle', 'none', 'TickStyle', 'dd', 'TickDir', 'out', 'FontSize', 15, ...
    'FontWeight', 'bold','FontName', 'Times');

end