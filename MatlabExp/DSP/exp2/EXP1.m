clc;
clear

N1 = 8;
N2 = 16;
N = N2; % Adjust N here.
n = 0:N-1;
x1 = ones(1, 4);
x2 = [1:4 4:-1:1];
x3 = [4:-1:1 1:4];
x4 = cos(pi/4 * n);
x5 = sin(pi/8 * n);


x_n = x3; % Adjust x_n here.
[k, x_n, ~] = fft_anylsis(x_n, N, 1);
figure
stem(n, x_n);
axis([0 N -inf inf])


N = 16;
n = 0:N-1;
D = 2*pi/N;
n = D*n;
x_n = x4+1j* x5;

X_k = fft(x_n);


% x_6
clc;
clear
f =@(t) cos(8*pi*t) + cos(16*pi*t) + cos(20*pi*t);
Fs = 64;
T = 1/Fs;
N = 64; % 32, 64

n = 0:N-1;
x_6 = f(n*T);

[k, x_6, mag_x, ang_x, ~] = fft_anylsis(x_6, N, 0);
mag_x = T*mag_x; ang_x = T*ang_x;
subplot(2, 2, 1)
fplot(f, [0, N*T])
title('f')

subplot(2, 2, 2)
stem(n, x_6)
title('x_6')

subplot(2, 2, 3)
plot(k*Fs, mag_x, '--');
axis([-Fs*pi, Fs*pi, -inf, inf])
xlabel('\Omega')
ylabel('|X(j\Omega)|')
title(['Fs=', num2str(Fs), '    N=', num2str(N)])
hold on
plot(k*Fs, mag_x, '*')
hold off

subplot(2, 2, 4)
plot(k*Fs, ang_x, '--');
axis([-Fs*pi, Fs*pi, -inf, inf])
xlabel('\Omega')
ylabel('arg[X(j\Omega)]')
title(['Fs=', num2str(Fs), '    N=', num2str(N)])
hold on
plot(k*Fs, ang_x, '*')
hold off





