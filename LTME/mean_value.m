% M-file edit에서 다음을 입력후 
% m-file: mean_value.m로 저장한다.
function mean_value=average(x)
[m,n]=size(x);
    if(~((m==1) | (n==1)) | (m==1& n==1))
    error('input must be a vector')
    end
mean_value=sum(x)/length(x);
