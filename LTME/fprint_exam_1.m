% fprint_exam.m
x=0:.2:1;         % x는 0부터 1까지 0.1씩 증가
y=[x;exp(x)];     % y는 x값에 따른 지수함수
fid=fopen('outdata.txt','w++');
fprintf(fid,'%6.2f%12.8f\n',y);
fclose(fid)