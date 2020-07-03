function [xep, xop] = circevod(x)
    N = length(x); n = 0:N-1;
    
    xep = .5*(x +x(mod(-n, N)+1));
    xop = .5*(x -conj(x(mod(-n, N)+1)));
end