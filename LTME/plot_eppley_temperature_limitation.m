close all; clear; clc;

figPos = [0 0 5 4];
fontSizeTick = 12;
fontSizeLab  = 12;
fontSizeCLab = 12;    

fennel_mu = 1.066.^[0:35];
tmax = 1.0 .* fennel_mu;
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos); hold on;
plot(0:35,tmax,'k--')
grid on
hold on
tmax_epp= 0.59 .* fennel_mu;
tmax0= 0.62 .* fennel_mu;
tmax2=2.0 .* fennel_mu;
tmax3=3.0 .* fennel_mu;
plot(0:35,tmax_epp,'r')
plot(0:35,tmax2,'k--')
plot(0:35,tmax3,'k--')
plot(0:35,tmax0,'k--')
xlim([0 35])

tm = 0.59 .* fennel_mu;
tm(1:24) = tm(25);

p=polyfit(0:35,tm,6);
y=polyval(p,0:35)

% plot(0:35,y,'m')

plot(0:35,tm,'b')

% https://aslopubs.onlinelibrary.wiley.com/doi/epdf/10.1002/lno.10523?src=getftr
% ni-tial  estimates  of  these  parameters  (a= 0.59  and b = 0.0633)have  been  updated  
% using  more  data  and  rigorous  quantileregression  methods (a = 0.81  and b=0.0613),  
% the values  weuse throughout the rest of this paper (Table 1; Bissinger et al.2008; 
% but see also Brush et al. 2002). 

umax = 0.81.*exp(0.0631.*[0:35]);
eppley = 0.59.*exp(0.0633.*[0:35]);

plot(0:35,umax,'g');
plot(0:35,eppley,'k--');


% legend({'1','2','3','0.62','6차 다항식'},'Location','northwest');
% legend({'1','2','3','0.62','new'},'Location','northwest');
ylabel('성장률(d^-1)','fontsize',12,'fontweight','bold');
xlabel('수온(^oC)','fontsize',12,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');
print(gcf,['growth_rate_compare.png'],'-dpng','-r500');
close;

xt = [0:35];
y2=round(p(1),12).* xt.^6 + round(p(2),12).* xt.^5 + round(p(3),12).* xt.^4 + ...
    + round(p(4),12).* xt.^3 + round(p(5),12).* xt.^2 + round(p(6),12).* xt + round(p(7),12);


plot(0:35,y,'k')
text(5,15,


1./(1+exp(-[0:35]))

1.066.^[0:35]

1.066.^[0:35]




