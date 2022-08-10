clc; clear all; close all;

% 과거 자료와 포맷을 맞추기 위한 코드

filename ='정선해양관측정보-20150101~20161231.xls';
sheet = 1;
[num,txt,raw] = xlsread(filename,sheet);

data = str2double(txt(3:end,[2 3 5 6 8 9 11]));
% 관측날은 txt로 유지시켜야 함

% 정점번호 처리
A = data(:,1);
B = str2num([num2char(data(:,2),2)]);
st_no = horzcat(A,B);

po_no = [num2char(data(:,2),2)]; % 2자리 수로 만들기 
po_no = str2double(po_no);

DateString = txt(:,7);
formatIn ='yyyy-mm-dd';
re = datenum(DateString,formatIn);

t = datetime(txt(:,7),'InputFormat','yyyy-MM-dd');