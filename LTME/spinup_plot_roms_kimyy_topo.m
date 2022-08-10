
bathyfile=[figdir, 'bathy_nwp']; 
status=plot_bathy(bathyfile,workdir, lonlat, year, [0 5000], [0]);
status=plot_bathy([figdir, 'bathy_taiwan'],workdir, [117.5 123 22 27], year, [0 130], [30 50 80]);
status=plot_bathy([figdir, 'bathy_korea'],workdir, [128 132 32.5 36], year, [0 200],0:20:200);
status=plot_bathy([figdir, 'bathy_tsugaru'],workdir, [139 142 40 43], year, [0 200], 0:20:200);
status=plot_bathy([figdir, 'bathy_soya'],workdir, [141 143 45 47], year, [0 200],0:20:200);
status=plot_bathy([figdir, 'bathy_kuro_5000'],workdir, [132 143 27 35], year, [0 5000],0:500:5000);
status=plot_bathy([figdir, 'bathy_kuro_2000'],workdir, [132 143 27 35], year, [0 2000],0:200:2000);
status=plot_bathy([figdir, 'bathy_kuro_500'],workdir, [132 143 27 35], year, [0 500],0:50:500);
status=plot_bathy([figdir, 'bathy_kuro_200'],workdir, [132 143 27 35], year, [0 200],0:20:200);
status=plot_bathy([figdir, 'bathy_luzon_5000'],workdir, [117 124 16 23], year, [0 5000],0:500:5000);
status=plot_bathy([figdir, 'bathy_luzon_2000'],workdir, [117 124 16 23], year, [0 2000],0:200:2000);
status=plot_bathy([figdir, 'bathy_luzon_500'],workdir, [117 124 16 23], year, [0 500],0:50:500);
status=plot_bathy([figdir, 'bathy_luzon_200'],workdir, [117 124 16 23], year, [0 200],0:50:200);
status=plot_bathy([figdir, 'bathy_kuro_5000_nc'],workdir, [132 143 27 35], year, [0 5000],0);
status=plot_bathy([figdir, 'bathy_kuro_2000_nc'],workdir, [132 143 27 35], year, [0 2000],0);
status=plot_bathy([figdir, 'bathy_kuro_500_nc'],workdir, [132 143 27 35], year, [0 500],0);
status=plot_bathy([figdir, 'bathy_kuro_200_nc'],workdir, [132 143 27 35], year, [0 200],0);
status=plot_bathy([figdir, 'bathy_luzon_5000_nc'],workdir, [117 124 16 23], year, [0 5000],0);
status=plot_bathy([figdir, 'bathy_luzon_2000_nc'],workdir, [117 124 16 23], year, [0 2000],0);
status=plot_bathy([figdir, 'bathy_luzon_500_nc'],workdir, [117 124 16 23], year, [0 500],0);
status=plot_bathy([figdir, 'bathy_luzon_200_nc'],workdir, [117 124 16 23], year, [0 200],0);
% status=plot_etopo1([figdir, 'bathy_etopo1_kuro_5000'],workdir, [132 143 27 35], year, [0 5000],0:500:5000);
% status=plot_etopo1([figdir, 'bathy_etopo1_kuro_2000'],workdir, [132 143 27 35], year, [0 2000],0:200:2000);
% status=plot_etopo1([figdir, 'bathy_etopo1_kuro_500'],workdir, [132 143 27 35], year, [0 500],0:50:500);
% status=plot_etopo1([figdir, 'bathy_etopo1_kuro_200'],workdir, [132 143 27 35], year, [0 200],0:50:200);
% status=plot_etopo1([figdir, 'bathy_etopo1_luzon_5000'],workdir, [117 124 16 23], year, [0 5000],0:500:5000);
% status=plot_etopo1([figdir, 'bathy_etopo1_luzon_2000'],workdir, [117 124 16 23], year, [0 2000],0:200:2000);
% status=plot_etopo1([figdir, 'bathy_etopo1_luzon_500'],workdir, [117 124 16 23], year, [0 500],0:50:500);
% status=plot_etopo1([figdir, 'bathy_etopo1_luzon_200'],workdir, [117 124 16 23], year, [0 200],0:50:200);
etopodir ='D:\MEPL\project\SSH\1st_year\figure\etopo1\';
status=plot_etopo1([etopodir, 'bathy_etopo1_nwp_5000'],workdir, [115 164 15 52], year, [0 5000],0);
status=plot_etopo1([etopodir, 'bathy_etopo1_kuro_5000'],workdir, [132 143 27 35], year, [0 5000],0);
status=plot_etopo1([etopodir, 'bathy_etopo1_kuro_2000'],workdir, [132 143 27 35], year, [0 2000],0);
status=plot_etopo1([etopodir, 'bathy_etopo1_kuro_500'],workdir, [132 143 27 35], year, [0 500],0);
status=plot_etopo1([etopodir, 'bathy_etopo1_kuro_200'],workdir, [132 143 27 35], year, [0 200],0);
status=plot_etopo1([etopodir, 'bathy_etopo1_luzon_5000'],workdir, [117 124 16 23], year, [0 5000],0);
status=plot_etopo1([etopodir, 'bathy_etopo1_luzon_2000'],workdir, [117 124 16 23], year, [0 2000],0);
status=plot_etopo1([etopodir, 'bathy_etopo1_luzon_500'],workdir, [117 124 16 23], year, [0 500],0);
status=plot_etopo1([etopodir, 'bathy_etopo1_luzon_200'],workdir, [117 124 16 23], year, [0 200],0);
close all;
