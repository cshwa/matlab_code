% text_add.m
x=0:.5:3;                   % x�� 0���� 5���� 0.5�� ������ ����
y=[x;sin(x)];                 % x�� ù �� sin(x)�� ��° ���� ���͸� ����
fid=fopen('exp_out.txt','w++');    % exp_out.txt�� �а� ���� ���� open
%'Input and output test' �Է� �� �ٹٲ�
fprintf(fid,'Input and output test\n');  
% ������ ǥ��, ������ x�� �Ҽ��� ���� 2 �ڸ� ���� 6������ ���� �Ҽ���ǥ�� 
% ���� ����, sin(x)�� ���� �Ҽ��� ���� 8�ڸ��� ���� 12���ڷ� ǥ��
count=fprintf(fid,'%6.2f%12.8f\n',y)
fclose('all');
fid=fopen('exp_out.txt','r+');   %��μ� �а� ���� ���� �־��� ������ ���� ����