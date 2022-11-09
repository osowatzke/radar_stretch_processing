% Author: O. Sowatzke
%
% Modified: 11/25/2022
%
% Subject: Class generates RDMs for a radar system using stretch processing
%
classdef radar < handle

    % public properties
    properties

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

    % private properties
    properties(Access=protected)

        % Nx1 vector of target ranges (m)
        R;

        % Nx1 vector of target velocities (m/s)
        v;

        % Nx1 vector of targets RCS (m^2)
        RCS;

        % Normalized thermal noise power expressed in units of 
        % ((4*pi)^3)/(Pt*G^2*lambda^2). No normalization factor if
        % Pt*G^2*lambda^2=(4*pi)^3
        thermal_noise_power;

        % Received data cube [fast time samples x slow time samples]
        data_cube;

        % Enable FFT window for range FFT
        en_range_win;

        % Enable R^2 Amplitude Compensation
        en_amp_comp;

        % Enable FFT window for doppler FFT
        en_dopp_win;

        % Range Doppler Matrix [range gate x doppler bin]
        rdm;

        % Nx1 vector of range gates
        range_gate;

        % Nx1 vector of doppler bins
        dopp_bin;
    end

    % Public methods
    methods

        % Class constructor. Can assign property values with comma 
        % separated list. Ex: "radar('f_c', 10e9)" creates a radar object 
        % that uses a 10GHz center frequency 
        function self = radar(varargin)
            for i = 1:2:nargin
                self.(varargin{i}) = varargin{i+1};
            end
        end

        % chirp bandwidth (Hz)
        function y = B(self)
            y = abs(self.f_stop - self.f_start);
        end

        % range bin size (m)
        function y = Rs(self)

            % time difference corresponding to Fast Time FFT bin
            dt = (self.f_s/self.M)*self.P_T/self.B;

            % convert to range
            y = self.c*dt/2;
        end

        % unambiguous range (m)
        function y = Rua(self)

            % Compute time corresponding to unambiguous range. Limiting 
            % factor is ADC rate instead of PRI in stretch processing
            t = self.f_s*self.P_T/self.B;

            % compute unambiguous range
            y = self.c*t/2;
        end

        % Function generates and plots the RDM for a radar system using
        % stretch processing. Windows and amplitude correction can be 
        % configured for different scenarios.
        %
        function generate_rdm(self,R,v,RCS,thermal_noise_power,...
                en_range_win,en_amp_corr,en_dopp_win)
    
            % save function inputs as private class properties
            % can be accessed in each of the class methods
            self.R = R;
            self.v = v;
            self.RCS = RCS;
            self.thermal_noise_power = thermal_noise_power;
            self.en_range_win = en_range_win;
            self.en_amp_comp = en_amp_corr;
            self.en_dopp_win = en_dopp_win;

            % generate received data cube
            self.simulate_received;
    
            % form rdm from received data cube
            self.process_received;
    
            % get location of each target in RDM
            self.find_rdm_peaks;

            % plot rdm
            figure;
            m = mesh(db(self.rdm),'FaceColor','flat');
            view(2);

            % mark each of the peaks
            for i = 1:length(self.range_gate)
                datatip(m, self.dopp_bin(i), self.range_gate(i));
            end

            % scale axis
            colorbar;
        
            % label plot
            xlabel('Doppler Bin');
            ylabel('Range Gate');
        end
    end

    % Private methods
    methods(Access=protected)

        % Function generates the received data cube 
        % [fast time samples x slow time samples] for a system employing
        % stretch processing
        function simulate_received(self)

            % timestep (s)
            dt = 1/self.f_s;
        
            % generate thermal noise for each data cube sample
            thermal_noise = normrnd(0,1,self.M,self.K) + ...
                1i.*normrnd(0,1,self.M,self.K);
        
            thermal_noise = thermal_noise*sqrt(self.thermal_noise_power);
        
            % compute frequency offet due to target range
            f_r = 2*self.R/self.c*(self.f_start - self.f_stop)/self.P_T;
        
            % frequency offset due to doppler shift
            f_d = 2*self.v*self.f_c/self.c;
        
            % array of times in (s)
            t = 0:dt:dt*(self.PRI*self.f_s*self.K-1);
        
            % generate received signal
            received = (sqrt(self.RCS)./(self.R.^2))'...
                *exp(1i.*2.*pi.*(f_r + f_d)*t);
        
            % generate datacube
            reshape_received = reshape(received,self.PRI*self.f_s,self.K);
            received = reshape_received(1:self.M,:);
            self.data_cube = received+thermal_noise;
        end

        % Function processes the received signal for a radar system using
        % stretch processing. Function output is an RDM.
        function process_received(self)

            % apply range window
            if self.en_range_win
                range_fft_in = self.data_cube.*hamming(self.M);
            else
                range_fft_in = self.data_cube;
            end

            % Apply range FFT
            range_fft_out = fft(range_fft_in,[],1);

            % FFT along range axis outputs frequencies from [0, 2*pi) for 
            % downchirp and [-2*pi, 0) for upchirp. 
            % 
            % If upchirp, must flip after FFT to place small ranges at 
            % bottom of data cube. This produces frequencies from 
            % (0, -2*pi]. Note that a range of zero is at end of RDM. 
            % Shift is needed to produce frequencies from [0, -2*pi) 
            % and place R=0 at bottom of RDM
            if self.f_start < self.f_stop
                range_fft_out = circshift(flip(range_fft_out,1),1,1);
            end

            % Array of ranges for amplitude correction
            R_vec = self.Rs*(0:(self.M-1)).';
        
            % Amply amplitude compensation
            % Accounts for loss in received signal amplitude
            if self.en_amp_comp
                dopp_fft_in = range_fft_out.*(R_vec.^2);
            else
                dopp_fft_in = range_fft_out;
            end
        
            % Apply doppler window
            if self.en_dopp_win
                dopp_fft_in = dopp_fft_in.*(hamming(self.K).');
            end
        
            % Apply doppler FFT
            self.rdm = fft(dopp_fft_in,[],2);
        end

        % Function computes the expected range gate and doppler bin for a
        % vectors of targets with range R and velocity V
        function find_rdm_peaks(self)
            
            % compute doppler shift
            fd = 2*self.v*self.f_c/self.c;

            % compute range shift due to target velocity
            % upchirp   => deltaR < 0
            % downchirp => deltaR > 0
            deltaR = self.c*fd*self.P_T/(2*self.B);
            if self.f_start < self.f_stop
                deltaR = -deltaR;
            end
    
            % compute adjusted range
            R_adj = self.R + deltaR;

            % compute ambiguous range
            Ra = mod(R_adj, self.Rua);

            % quantize to a range gate
            % add one to produce MATLAB indexing
            self.range_gate = mod(round(Ra/self.Rs),self.M) + 1;

            % normalized frequency fd/PRF
            fnorm = mod(fd*CONSTANTS.PRI,1);

            % doppler shift due to frequency offset
            % add one to produce MATLAB indexing
            self.dopp_bin = mod(round(fnorm*CONSTANTS.K),CONSTANTS.K) + 1;
        end
    end
end
