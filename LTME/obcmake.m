clear M
for j=1:16
    for i=1:16
        M(j+(i-1)*(16),2)=j;
        M(j+(i-1)*(16),1)=i;
        M(j+(i-1)*(16),3)=0;
    end
end

%land
%strait
M(33,3)=3,0;
M(34,3)=3.0;
M(147,3)=3.0;
M(148,3)=3.0;
dlmwrite('Vel_U.dat', M)