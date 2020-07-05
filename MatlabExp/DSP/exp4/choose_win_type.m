function w = choose_win_type(win_type, N)
    switch win_type
        case 'hn'
            w = hanning(N);
        case 'hm'
            w = hamming(N);
        case 'r'
            w = rectwin(N);
        case 'tr'
            w = triang(N);
    end
end