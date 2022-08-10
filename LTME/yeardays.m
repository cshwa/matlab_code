function d = yeardays(year);

if mod(year,4) == 0
    if mod(year,100) == 0
        d = 365;
    else
        d = 366;
    end
else
    d = 365;
end

return;