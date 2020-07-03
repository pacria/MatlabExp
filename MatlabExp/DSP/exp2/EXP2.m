% EXP-2
N = 16;
x_n = x4+x5;

[k, x_n, mag_X_k, ang_X_k, re_X_k, im_X_k] = fft_anylsis(x_n, N, 0);
[k, x4a, mag_x4, ~] = fft_anylsis(x4, N, 0);
[k, x5a, mag_x5, ~] = fft_anylsis(x5, N, 0);

subplot(2, 2, 1)
plot(k, re_X_k)
title('|Re[X(k)]|')

subplot(2, 2, 2)
plot(k, abs(im_X_k))
title('|Im[X(k)]|')

subplot(2, 2, 3)
plot(k, mag_x4)
title('|X_(ep)(k)|(|X_4(k)|)')

subplot(2, 2, 4)
plot(k, mag_x5)
title('|X_(op)(k)|(|X_5(k)|)')
