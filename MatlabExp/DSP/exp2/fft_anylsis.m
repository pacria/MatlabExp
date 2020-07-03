function [k, x_n, mag_X_k, ang_X_k, re_X_k, im_X_k, X_k] = fft_anylsis(x_n, N, anylsis)
    M = length(x_n);
    x_n = [x_n zeros(1, N-M)];
    
    X_k = fftshift(fft(x_n));
    X_k(N+1) = X_k(1); % 多补充一位以利于体现对称性（在n=(N+1)/2 处补上 n=-(N-1)/2的值 相当于将序列在圆周上\pi处排列的值重复一遍）
    
    k = floor(-(N-1)/2:(N+1)/2);
    k = k .* (2 * pi / N);
    
    mag_X_k = abs(X_k);
    ang_X_k = angle(X_k);
    
    re_X_k = real(X_k);
    im_X_k = imag(X_k);
    
    if anylsis
        subplot(1, 2, 1)
        plot(k, mag_X_k, '--');
        axis([-pi, pi, -inf, inf]);
        xlabel('\omega')
        ylabel('|X(k)|')
        title(['mag   ', 'N=', num2str(N)])
        hold on 
        plot(k, mag_X_k, '*');
        hold off 
        
        subplot(1, 2, 2)
        plot(k, ang_X_k, '--');
        axis([-pi, pi, -inf, inf]);
        xlabel('\omega')
        ylabel('arg[X(k)]')
        title(['ang   ', 'N=', num2str(N)])
        hold on 
        plot(k, ang_X_k, '*');
        hold off
    end
end
