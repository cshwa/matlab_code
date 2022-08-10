close all; clear all;
% Example for Coherency spectrum
[seq,year,month,day,hour,min, temp, flt_temp, u1, v1, flt_u, flt_v, flt_spd, flt_drc, rot_u, rot_v] ...
    = textread('wt9395.flt_Korea_strait.dat','%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','headerlines',1);

[freq, coh, ph, ci, phi] = cmtm(u1,v1,1);
semilogx(freq,coh)
axis([0.0001 0.1 0 1])
hold on
semilogx(freq,ci,'r-') % 95% confidence level
xlabel('Frequency (CPH)','fontname','times','fontsize',14,'fontweight','bold');
ylabel('Coherency','fontname','times','fontsize',14,'fontweight','bold');

figure
semilogx(freq, ph)
axis([0.0001 0.1 -30 30])
xlabel('Frequency (CPH)','fontname','times','fontsize',14,'fontweight','bold');
ylabel('Phase (deg.)','fontname','times','fontsize',14,'fontweight','bold');
