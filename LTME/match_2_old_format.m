clc; clear all; close all;

% ���� �ڷ�� ������ ���߱� ���� �ڵ�

filename ='�����ؾ��������-20150101~20161231.xls';
sheet = 1;
[num,txt,raw] = xlsread(filename,sheet);

data = str2double(txt(3:end,[2 3 5 6 8 9 11]));
% �������� txt�� �������Ѿ� ��

% ������ȣ ó��
A = data(:,1);
B = str2num([num2char(data(:,2),2)]);
st_no = horzcat(A,B);

po_no = [num2char(data(:,2),2)]; % 2�ڸ� ���� ����� 
po_no = str2double(po_no);

DateString = txt(:,7);
formatIn ='yyyy-mm-dd';
re = datenum(DateString,formatIn);

t = datetime(txt(:,7),'InputFormat','yyyy-MM-dd');