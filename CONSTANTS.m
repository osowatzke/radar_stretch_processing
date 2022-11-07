classdef CONSTANTS
    properties(Constant)
        % center frequency (Hz)
        f_c = 77e9;
        % chirp starting frequency (Hz)
        f_start = 76.8e9;
        % chirp ending frequency (Hz)
        f_stop = 77.2e9;
        % pulse repetition interval (s)
        PRI = 30e-6;
        % Chip pulse width (s)
        P_T = 20e-6;
        % fast time samples
        M = 400;
        % number of chips
        K = 500;
        % ADC sample frequency (Hz)
        f_s = 20e6;
        % light speed (m/s)
        c = 3e8;
    end
    methods(Static)

        % chirp bandwidth (Hz)
        function y = B
            y = abs(CONSTANTS.f_stop - CONSTANTS.f_start);
        end

        % range bin size (m)
        function y = Rs

            % time difference corresponding to Fast Time FFT bin
            dt = (CONSTANTS.f_s/CONSTANTS.M)*CONSTANTS.P_T/CONSTANTS.B;

            % convert to range
            y = CONSTANTS.c*dt/2;
        end

        % unambiguous range (m)
        function y = Rua

            % Compute time corresponding to unambiguous range. Limiting 
            % factor is ADC rate instead of PRI in stretch processing
            t = CONSTANTS.f_s*CONSTANTS.P_T/CONSTANTS.B;

            % compute unambiguous range
            y = CONSTANTS.c*t/2;
        end
    end
end