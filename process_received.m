% Author: O. Sowatzke
%
% Updated: 11/08/2022
%
% Subject: Function processes the received signal for a radar system using
% stretch processing. Function output is an RDM.
%
% Inputs:
%       data_cube           Received data cube [fast time samples x 
%                           slow time samples]
%       en_range_win        Enable range window (Fast Time FFT Window)
%       en_amp_comp         Enable R.^2 Amplitude Compensation
%       en_dopp_win         Enable doppler window (Slow Time FFT Window)
%
% Outputs:
%       rdm                 Range Doppler Matrix [range gate x doppler bin]
%
function rdm = process_received(data_cube, en_range_win, en_amp_comp, en_dopp_win)

    % Apply Range Window
    if en_range_win
        data_cube = data_cube.*hamming(CONSTANTS.M);
    end

    % Apply range FFT
    data_cube = fft(data_cube,[],1);

    % Array of ranges for amplitude correction
    R = CONSTANTS.Rs*(0:(CONSTANTS.M-1)).';

    % Range axis is flipped for upchirp
    % R > 0 => f_r < 0
    if CONSTANTS.f_start < CONSTANTS.f_stop
        R = flip(R);
    end

    % Amply amplitude compensation
    % Accounts for loss in received signal amplitude
    if en_amp_comp
        data_cube = data_cube.*(R.^2);
    end

    % Apply doppler window
    if en_dopp_win
        data_cube = data_cube.*(hamming(CONSTANTS.K).');
    end

    % Apply doppler FFT
    rdm = fft(data_cube,[],2);
end