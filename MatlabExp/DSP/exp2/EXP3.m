[k, x4a, mag_x4, ~] = fft_anylsis(x4, N, 0);
[k, x5a, mag_x5, ~] = fft_anylsis(x5, N, 0);

[xep, xop] = circevod(X_k);
xep = fftshift(xep);
xop = fftshift(xop);

subplot(2, 2, 1)
plot(n, abs(xep))
title('|X_ep (k)|')

subplot(2, 2, 2)
plot(n, abs(xop))
title('|X_op (k)|')

subplot(2, 2, 3)
plot(k, mag_x4)
title('|DFS(Re[x(n)])|(|X_4(k)|)')

subplot(2, 2, 4)
plot(k, mag_x5)
title('|DFS(jIm[x(n)])|(j*|X_5(k)|)')
