
close all
clear all
clc

kw2 = [200304; 200307; 200310; 200312];
for dd = 1:length(kw2)
    st1 = [3, 5, 7, 11, 12, 14, 18, 19, 21, 24, 25, 28, 29, 31, 32, 34, 35, 36, 38, 39, 41];
    if dd == 4
        st1 = [11, 12, 14, 18, 19, 21, 24, 25, 28, 29, 31, 32, 34, 35, 36, 38, 39, 41];
    end
    ll1 = length(st1);
    for i = 1:ll1
        aa = num2char(st1(i),2);
        if dd == 1
            sheet = (['kw',aa]);
        elseif dd == 2
            sheet = (['HD',aa]);
        elseif dd == 3
            if i < 4
                aa = num2char(st1(i),1);
            end
            sheet = (['hd',aa]);
        else
            sheet = (['hd',aa]);
        end
        file = ([num2char(kw2(dd),6),'하동.xls']);
        output = (['./',num2char(kw2(dd),6),'./',sheet,'.dat']);
        
%         [num,txt,raw] = xlsread(file,sheet,'a4');
        data = xlsread(file,sheet,'b5:i100');
        ind = find(isnan(data));
        data(ind) = 0;
        fid = fopen(output,'w');
%         fprintf(fid,'%20s\n',nam);
        fprintf(fid,'%5d%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n',data');
        fclose(fid);
    
    end
end

for dd = 1:length(kw2)
    st2 = [42, 43, 44, 46, 47, 49, 54, 56, 57, 58, 59, 60, 62, 64, 66, 68, 70, 72, 74, 78, 81, 84, 86, 87];
    if dd == 3
        st2 = [42, 43, 44, 47, 49, 56, 57, 58, 59, 62, 64, 70];
    elseif dd == 4
        st2 = [42, 43, 44, 47, 49, 54, 56, 57, 58, 59, 60, 62, 64, 66, 68, 70, 72, 74, 77, 78, 81, 84, 86];
    end
    ll2 = length(st2);
    for i = 1:ll2
        aa = num2char(st2(i),2);
        if dd == 1
            sheet = (['gj',aa]);
        elseif dd == 2
            sheet = (['HD',aa]);
        elseif dd == 3
            sheet = (['hd',aa]);
        else
            sheet = (['hd',aa]);
        end
        file = ([num2char(kw2(dd),6),'하동.xls']);
        output = (['./',num2char(kw2(dd),6),'./',sheet,'.dat']);
        
%         [num,txt,raw] = xlsread(file,sheet,'a4');
        data = xlsread(file,sheet,'b5:i100');
        ind = find(isnan(data));
        data(ind) = 0;
        fid = fopen(output,'w');
%         fprintf(fid,'%20s\n',nam);
        fprintf(fid,'%5d%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n',data');
        fclose(fid);
    
    end
end
