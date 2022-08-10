% plot the adcp data as direction arrows.  Loads the .mat files created by
% adcpdata2array.m.  Must be run from the same directory the files are in.
% Bill Scuba   Oct. 15, 2003 and edit by Peter(2004-4-23)
% INPUT data: ���� ���丮�� �ִ� ��� *.mat ����
% m-file: plotadcp1.m
% array���� ����� *.mat������ �о� �迭 ping�� �����͸� �ط����ͷ� ����
% ���� adcpdata2array.m�� ���� �� �� ���������� �����ϸ� �����ط��������� �ۼ�
% 

clear dirstruct            % �迭 dirstruct�� ����� ������ ����

figure                     % �׸� â (1��)�� ����
hold on                    % ��� �׸��� ������ â�� �׸���
dirstruct=dir('*.mat')     % ���� ���丮�� Ȯ���ڰ� mat�� ��� ������ �迭 dirstruct�� ����

% ���� ���丮�� Ȯ���ڰ� mat�� ��� ������ �ݺ� ó���ϴ� ��ƾ
	for yi=1:length(dirstruct)   % yi�� 1���� �迭 dirstruct�� ���ϰ��� ��ȣ����
    clear ping                       %�迭 ping�� ���� ������ ���� 
    load(dirstruct(yi).name) 
	for xi=1:50:length(ping)     % xi�� 1����� ping�� ��ü�� ��ȣ���� 50��° �� ��ȣ���� ó���ϴ� ����
		quiver(ping(xi).lon,ping(xi).lat,ping(xi).adcpdata(1,4),ping(xi).adcpdata(1,5),.002)
		% �迭 ping�� xi��° ���� �̸��� �迭 lon (x��: �浵��), 
		% �迭 ping�� xi��° ���� �̸��� �迭 lat (y��: ������),
		% �迭 ping�� xi��° ���� �̸��� �迭 adcpdata�� 1��, 4���� �� (���漺��)
		% �迭 ping�� xi��° ���� �̸��� �迭 adcpdata�� 1��, 5���� �� (�Ϲ漺��))
		% ȭ��ǥ�� ���� (length of arrow)
	end
    fprintf('u = %g \n',sqrt(ping(xi).adcpdata(1,4)^2 +ping(xi).adcpdata(1,5)^2))
		% u ���� ���� ���漺��(U)�� �Ϲ漺��(V)�� ������ ������ ���Ѱ��� �������� ���� ����
    end
hold off           % ��� �׸��� ������ â�� �׸��� ����