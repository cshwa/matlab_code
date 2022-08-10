clf
% [t,eu1800,eu2120,ku1800,ku2283,ea0212,ea0363,ea1037,ea1378,ea1745,ea2314,kt1800,kt2158,et2085] = ...
% textread('flt.12h.0.dat','%f%f%f%f%f%f%f%f%f%f%f%f%f%f','headerlines',1);
% [t,infw,oufw,tran] = ...
% textread('UIG_inoutflow.dat','%f%f%f%f','headerlines',1);
% dat = load('wt9395.flt.dat');
% % u = infw;
% u = dat(:,15);  ouf = 'inflow';
% y = dat(:,2); m = dat(:,3); d = dat(:,4); h = dat(:,5);
% y = y + 1900;
% jo = date2jd(1990,1,1);
% t = (date2jd(y,m,d,h,0,0)-jo)*24;
dt = 1; % hourly

% clear eu1800 eu2120 ku1800 ku2283 ea0212 ea0363 ea1037 ea1378 ea1745 ea2314 kt1800 kt2158 et2085
load wavelet_exam.mat

iu = find(u > 999);
u(iu) = [];
t(iu) = [];
ti = [min(t):dt:max(t)];
ui = interp1(t,u,ti,'cubic');

variance = std(ui)^2;
uit = (ui - mean(ui))/sqrt(variance) ;
cor1 = corrcoef(uit(1:end-1), uit(2:end));
cor2 = corrcoef(uit(1:end-2), uit(3:end));
cor = (cor1(1,2)+sqrt(cor2(1,2)))/2;

n = length(uit);
dt = dt ;
time = ti;  % construct time array
xlim = [min(t),max(t)];  % plotting range
pad = 1;      % pad the time series with zeroes (recommended)
dj = 0.05;    % this will do 4 sub-octaves per octave (default = 0.25)
s0 = 6*dt;    % this says start at a scale of 1 hour
j1 = 7/dj;    % this says do 7 powers-of-two with dj sub-octaves each
lag1 = cor;  % lag-1 autocorrelation for red noise background
mother = 'Morlet';

% Wavelet transform:
[wave,period,scale,coi] = wavelet(uit,dt,pad,dj,s0,j1,mother);
power = (abs(wave)).^2 ;        % compute wavelet power spectrum

% Significance levels: (variance=1 for the normalized SST)
[signif,fft_theor] = wave_signif(1.0,dt,scale,0,lag1,0.95,-1,mother);
sig95 = (signif')*(ones(1,n));  % expand signif --> (J+1)x(N) array
sig95 = power ./ sig95;         % where ratio > 1, power is significant

levels = [1 2 4 8 12 16 24 48 96] ;
Yticks = [1 2 4 8 12 16 24 48 96]*24;
% Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
% contourf(time,log2(period),log2(power),log2(levels));  %*** or use 'contourfill'
% jo = date2jd(1996,1,1);
% ddy = dyear(t+date2jd(1900,1,1));
% mon = [1:12, 1 2 3];
% year = ones(1,12)*2003;
% year(end+1:end+3) = [1993 1994 1995];
% day = ones(1,15)*1;
% jjj = date2jd(year,mon,day)-jo;
% Xticks = jjj;
% Xticklab = ['Jan03'; 'Feb03'; 'Mar03'; 'Apr03'; 'May03'; 'Jun03'; 'Jul03'; 'Aug03'; 'Sep03'; 'Oct03'; 'Nov03'; 'Dec03'; 'Jan04'; 'Feb04'; 'Mar04'];
 
imagesc(time,log2(period),log2(power));  %*** uncomment for 'image' plot
caxis([-20 10])
xlabel('Time ')
ylabel('Period ')
% title([' Wavelet Power Spectrum ' ouf ' Velocity component of North-Eastern direction'])
% title([' Wavelet Power Spectrum of ' ouf ' through UIG'])
% title('Meridional Velocity Wavelet Power Spectrum')
set(gca,'XLim',xlim(:))
set(gca,'YLim',log2([min(period),max(period)]), ...
	'YDir','normal',...
	'YTick',log2(Yticks(:)), ...
	'YTickLabel',levels(:));
% set(gca,'YLim',log2([min(period),max(period)]), ...
% 	'YDir','normal', ...
% 	'YTick',log2(Yticks(:)), ...
% 	'YTickLabel',Yticks, ...
% 	'XTick',jjj, ...
% 	'XTickLabel',Xticklab)
% );
% significance contour, levels at -99 (fake) and 1 (95% signif)
hold on
contour(time,log2(period),sig95,[-99,1],'k');
hold on
% cone-of-influence, anything "below" is dubious
plot(time,log2(coi),'k')

% print('-depsc',[ouf '.eps'])

