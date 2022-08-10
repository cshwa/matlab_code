% m-file: colorbar_exam.m
% Create a figure with a jet colormap
%
clc;clear all;close all
colormap jet
% Create a standard colorbar
set(axes,'Visible','off')
bar = colorbar;
set(get(bar,'Title'),'String','Standard Colorbar')
pos = get(bar,'Position');
  
% Create a terrain height scale using part of the colormap
colorscale([31 54], [0 9000], 1000, 'vert',...
            'Position',pos - [0.7 0 0 0])
title('Terrain Height')
ylabel('elevation above sea level, meters')

% Create a bathymetry scale that uses a different part of the
% colormap and runs in reverse (i.e., CMAPLIM(1) > CMAPLIM(2)).
colorscale([24 1], [0 12000], 1000, 'vert',...
          'Position',pos - [0.35 0.1 0.02 0],'YDir','reverse')
title('Bathymetry')
ylabel('depth below sea level, meters')