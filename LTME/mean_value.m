% M-file edit���� ������ �Է��� 
% m-file: mean_value.m�� �����Ѵ�.
function mean_value=average(x)
[m,n]=size(x);
    if(~((m==1) | (n==1)) | (m==1& n==1))
    error('input must be a vector')
    end
mean_value=sum(x)/length(x);
