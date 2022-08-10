% movie_exam1.m
%  �� �ڵ�� AVI������ matlab figureâ���� display�ϸ鼭 
%  1-D sine,cosine�� �����ϴ� ������ ���� �ۼ�
%  

clc;clear all;close all;

% Frame ���� ������ �ð�����
numframes = 50;               % Frames�� ����
seg = 2 * pi / numframes;     % 0~2pi������ frame���� ���� segment
FramesPerSecond = 10;         % �ʴ� ������ frame �� (fps;fames per second)
totaltime = numframes/FramesPerSecond % ������ �ɸ��� �ð�  

% AVIFILE �Լ��� �̿��� AVI���� ����
aviobj = avifile ('sin&cos.avi','fps', FramesPerSecond,'quality',100); 

% X�� ���� ������ ����
i = 1 : numframes + 1;
x = ( i - 1 )*seg;

%  ���� ���ӵ� sine�� cosine ������ ������ ����
% for���� addframe����� ���� AVI ������������ �������
for j =  1 : numframes+1
    y1 = sin ( x - ( j - 1 ) * seg );    % sine wave�� red ������ ǥ��
    plot ( x, y1,'r');hold on
    y2 = cos( x - ( j - 1 ) * seg );    % cosine wave�� blue ǥ��
    plot ( x, y2,'b');hold on
    y3 = y1+y2;                        % sine+cosine wave�� green ǥ��
    plot ( x, y3,'g','LineWidth',4); 
    legend('sine','cos','sine+cosine')
    plot([0,7],[0,0],'k:');hold on;
    title('trigonometric function display');
    frame = getframe (gca ); hold off;
    aviobj = addframe (aviobj, frame );   % addframe�� �̿��ؼ� ��� �׸���
end
%  Movie �۾��� ����
aviobj = close (aviobj);
fprintf ( 1, 'You have a good JOB !\n' );  % �۾� ���Ῡ�� ǥ��