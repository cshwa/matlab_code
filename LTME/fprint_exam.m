% fprint_exam.m
x=0:.2:1;         % x�� 0���� 1���� 0.1�� ����
y=[x;exp(x)];     % y�� x���� ���� �����Լ�
fid=fopen('outdata.txt','w++');
fprintf(fid,'%6.2f%12.8f\n',y);
fclose(fid)