% text_add.m
x=0:.5:3;                   % x는 0부터 5까지 0.5씩 증가한 벡터
y=[x;sin(x)];                 % x가 첫 열 sin(x)가 둘째 열인 벡터를 생성
fid=fopen('exp_out.txt','w++');    % exp_out.txt를 읽고 쓰기 위해 open
%'Input and output test' 입력 후 줄바꿈
fprintf(fid,'Input and output test\n');  
% 개수를 표시, 서식은 x를 소수점 이하 2 자리 전부 6문자의 고정 소수점표현 
% 공백 삽입, sin(x)의 값을 소수점 이하 8자리로 전부 12문자로 표현
count=fprintf(fid,'%6.2f%12.8f\n',y)
fclose('all');
fid=fopen('exp_out.txt','r+');   %비로소 읽고 쓰기 위해 주어진 파일을 열고 저장