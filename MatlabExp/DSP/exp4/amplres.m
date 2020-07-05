function [Hw, w, type] = amplres(h, M)
    
    N = length(h);
    L = floor((N-1)/2);
    h = h(:)';
    
    n = 1:L+1;
    w = (0:M)* 2 *pi /(M+1);
    
    if all(abs(h(n) - h(N - n + 1))<1e-10)
        % Class I(Type I or Type II)
        Hw = 2 * h(n) * cos(((N+1)/2 - n)' * w) - mod(N, 2) * h(L+1);
        type = 2 - mod(N, 2); %N - Even, type->2, N - Odd, type -> 1.
    elseif all(abs(h(n) - h(N - n + 1))<1e-10)&&(h(L+1) * mod(N, 2) == 0)
        Hw = 2 * h(n) * sin(((N+1)/2 - n)' * w);
        type = 4 - mod(N, 2); %N - Even, type->4, N - Odd, type -> 3.
        
    else error('Not a linear phase FIR filter');
    end
end