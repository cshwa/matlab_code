% �� ���� ������ �ֽ� ����(r2016a)������ ���� - �ӵ� ����
clear all
clc

% fname = 'F:\ROMS\Sumjin\01_2dm(basemap)\grid\Sumjin_roms_v11.2dm';
fname = '..\01_2dm(basemap)\grid\Sumjin_roms_v11.2dm';
[C a1 a2 a3 a4  ] = textread(fname,'%s  %f %f %f %f %*[^\n]', 'headerlines',2); % ���� *.2dm ������ headerlines �� 2����

data = [a1 a2 a3 a4];

cc = char(C);
for i = 1:length(C)
    if(cc(i,:) == 'ND ')
    idcc(i,:) = 1;
    end
end
ia = find(idcc ==1);

sub_grid = data(min(ia):max(ia),:);
% lm = a2(1) -1 ;
lm = a2(1) -2 ;
ln =  length(sub_grid) / lm ;

grid = reshape(sub_grid(:,1),lm,ln)';
lon_rho = reshape(sub_grid(:,2),lm,ln)';
lat_rho = reshape(sub_grid(:,3),lm,ln)';

% 
% save ideal_sumjin_v02.mat lon_rho lat_rho
% save e_v02_grid.mat lon_rho lat_rho
