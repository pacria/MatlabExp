function [db, mag, pha, grd, w] = freqz_m(b, a, M)
    [H, w] = freqz(b, a, M*2, 'whole');
    H = H(1:M)'; w = w(1:M)';
    mag = abs(H);
    db = 20*log10((mag+eps)/max(mag));
    pha = angle(H);
    grd = grpdelay(b, a, w);
    
end