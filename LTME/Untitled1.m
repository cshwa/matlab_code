temp1 = textread('temp1.txt','');
temp2 = textread('temp2.txt','');
salp1 = textread('salp1.txt','');
salp2 = textread('salp2.txt','');

temp1(:,2) = -1 * abs(temp1(:,2));
temp2(:,2) = -1 * abs(temp2(:,2));
salp1(:,2) = -1 * abs(salp1(:,2));
salp2(:,2) = -1 * abs(salp2(:,2));


fid = fopen('temp1.txt','w');
fprintf(fid, '%f %f %f %f\n', temp1');
fclose(fid);
fid = fopen('temp2.txt','w');
fprintf(fid, '%f %f %f %f\n', temp2');
fclose(fid);
fid = fopen('salp1.txt','w');
fprintf(fid, '%f %f %f %f\n', salp1');
fclose(fid);
fid = fopen('salp2.txt','w');
fprintf(fid, '%f %f %f %f\n', salp2');
fclose(fid);
