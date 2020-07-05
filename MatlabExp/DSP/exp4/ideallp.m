function [hd, tau] = ideallp(wc, N)
% Ideal low-pass digital filter(with cut-off angular frenquency wc)
    tau = (N-1)/2;
    n = 0:N-1;
    m = n - tau + eps; % Avoid ZeroDivision Error
    hd = sin(wc * m)./(pi * m);
end