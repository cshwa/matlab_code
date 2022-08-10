
function [f_final,s_final,phase]=fftspectrum(u,t)

%FFT 
% fft(1:10)�� �����ϸ�, 10���� ���� �����µ�, 1���� 0 frequency, 2~5���� + frequency, 6~10����
% - frequency�� �ǹ��Ѵ�. ���� ���ļ��� - , 0 , + ������ ������ �Ϸ���,
% fftshift(fft(1:10))���� ���ָ� �ȴ�. 0 frequency�� 1:10������ ������ �ǹ�..

% load tjrwuv2yr.mat

[ns, nt] = size(u);

df = (nt-1)/(max(t)-min(t)); %unit : cycles per day

if mod(nt,2), f=df*[-(nt-1)/2:(nt-1)/2]/nt;
else, f = df*[-nt/2:nt/2-1]/nt;
end
deltaf=diff(f(1:2));

u=u-nanmean(u);      % 0 frequency�� ���Գ��� �׸����°� ����� ���� �־� ���ش�. 

ii=isnan(u);u(ii)=0; %addhoc
%ii = find(isnan(u) | isnan(v));

d=u;
fc=fftshift(fft(d))/nt;
phase = angle(fc);
s=abs(fc).^2/deltaf;

f_final=f(0<f);
s_final=s(0<f);



%paseval's theorem

% sum(s)*deltaf;  %spetrum�� �� ����(������)�� data�� �л갪�� ��ġ�Ѵ�. ��, s���� ��հ��� ����� ����.
% nanvar(u);      %�л��� ��հ��� �������� �ʱ� ������
end

