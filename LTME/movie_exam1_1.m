% movie_exam1.m
%  이 코드는 AVI파일을 matlab figure창에서 display하면서 
%  1-D sine,cosine을 재현하는 동영상 파일 작성
%  

clc;clear all;close all;

% Frame 수와 재현할 시간지정
numframes = 50;               % Frames의 갯수
seg = 2 * pi / numframes;     % 0~2pi범위를 frame수로 나눈 segment
FramesPerSecond = 10;         % 초당 재현할 frame 수 (fps;fames per second)
totaltime = numframes/FramesPerSecond % 재현에 걸리는 시간  

% AVIFILE 함수를 이용한 AVI파일 생성
aviobj = avifile ('sin&cos.avi','fps', FramesPerSecond,'quality',100); 

% X축 값과 간격을 조정
i = 1 : numframes + 1;
x = ( i - 1 )*seg;

%  이제 연속된 sine과 cosine 동영상 파일을 생성
% for문과 addframe명령을 통해 AVI 동영상파일이 만들어짐
for j =  1 : numframes+1
    y1 = sin ( x - ( j - 1 ) * seg );    % sine wave를 red 색으로 표현
    plot ( x, y1,'r');hold on
    y2 = cos( x - ( j - 1 ) * seg );    % cosine wave를 blue 표현
    plot ( x, y2,'b');hold on
    y3 = y1+y2;                        % sine+cosine wave를 green 표현
    plot ( x, y3,'g','LineWidth',4); 
    legend('sine','cos','sine+cosine')
    plot([0,7],[0,0],'k:');hold on;
    title('trigonometric function display');
    frame = getframe (gca ); hold off;
    aviobj = addframe (aviobj, frame );   % addframe를 이용해서 계속 그리기
end
%  Movie 작업을 종료
aviobj = close (aviobj);
fprintf ( 1, 'You have a good JOB !\n' );  % 작업 종료여부 표현