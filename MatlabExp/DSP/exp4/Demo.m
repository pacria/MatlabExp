% Pre-defined
N1 = 15;
N2 = 33;
N = N2;
omegaC = pi/4;
hanning = 'hn';
hamming = 'hm';
rect = 'r'; % Do nothing(Because it's already fixed to N point)
triang = 'tr'
win_type = hamming; % Choose one type (DEFAULT: hanning)

bianer = 3;
M = 2^ceil(log2(N)+bianer);
n = 0:N-1;

[hd, tau] = ideallp(omegaC, N);

[hdw, w, ~] = amplres(hd, M);

figure 
subplot(2, 1, 1)
stem(n, hd)
xlabel('n')
ylabel('h(n)')

subplot(2, 1, 2)
plot(w/pi, hdw)

grid on

hold on
plot([omegaC/pi omegaC/pi], ylim)
hold off
xlabel('\omega/\pi')

legend('|H_d(\omega)|(Rect)', 'Cut-off frenquency')

% Get window
wn = choose_win_type(win_type, N)';
[Hww, w, ~] = amplres(wn, M);

figure
subplot(2, 1, 1)
stem(n, wn)
xlabel('n')
ylabel('w(n)')

subplot(2, 1, 2)
plot(w/pi, Hww)

grid on

xlabel('\omega/\pi')
ylabel('|W(\omega)|')

% Get -h 
h = wn.*hd; % Get your window


[hw, w, ~] = amplres(h, M);

figure
subplot(2, 1, 1)
stem(n, hd)

hold on
stem(n, h)
hold off

xlabel('n')

legend('H_d(n)', 'H(n)')

subplot(2, 1, 2)
plot(w/pi, hw)

grid on

hold on
plot([omegaC/pi omegaC/pi], ylim)
plot(w/pi, hdw)
hold off

xlabel('\omega/\pi')

legend('|H(\omega)|', '|H_d(\omega)|(Rect)', 'Cut-off frenquency')

[db, mag, pha, grd, w] = freqz_m(h, [1], M);
    
figure 
subplot(2, 2, 1)
stem(n, h)
xlabel('n')
ylabel('H(n)')

subplot(2, 2, 2)
plot(w/pi, db)
grid on 
xlabel('\omega/\pi');
ylabel('20log(|H(w)|)/dB')

subplot(2, 2, 3)
plot(w/pi, mag)
grid on
xlabel('\omega/\pi')
ylabel('|H(w)|')

subplot(2, 2, 4)
plot(w/pi, pha)
xlabel('\omega/\pi')
ylabel('arg(|H(w)|)')







