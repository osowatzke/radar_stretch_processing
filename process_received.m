function rdm = process_received(data_cube, en_range_win, en_amp_comp, en_dopp_win)

    % array of ranges
    R = flip(CONSTANTS.Rs*(0:(CONSTANTS.M-1)).');

    % apply range window
    if en_range_win
        data_cube = data_cube.*hamming(CONSTANTS.M);
    end

    % apply range FFT
    data_cube = fft(data_cube,[],1);

    % amply amplitude compensation
    % accounts for loss in received signal amplitude
    if en_amp_comp
        data_cube = data_cube.*(R.^2);
    end

    % apply doppler window
    if en_dopp_win
        data_cube = data_cube.*(hamming(CONSTANTS.K).');
    end

    % apply doppler FFT
    rdm = fft(data_cube,[],2);
end