clear M
for j=1:16
    for i=1:16
        M(j+(i-1)*(16),2)=j+30;
        M(j+(i-1)*(16),1)=i+124;
        M(j+(i-1)*(16),3)=2000;
    end
end

%land
for i=1:32
    M(i,3)=0;
end
for i=225:256
    M(i,3)=0;
end
for i=33:16:241
    M(i,3)=0;
end
for i=34:16:242
    M(i,3)=0;
end
for i=15:16:127
    M(i,3)=0;
end
for i=16:16:128
    M(i,3)=0;
end
for i=159:16:255
    M(i,3)=0;
end
for i=160:16:256
    M(i,3)=0;
end

%strait
M(33,3)=200;
M(34,3)=200;
M(147,3)=200;
M(148,3)=200;
dlmwrite('XYZ_BathymetryData_Box.xyz', M)