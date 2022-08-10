
function f = fSSR(a, xm, ym)
    yp = a(1)*xm.^a(2);
    F = sum((ym-yp).^2);
