% Author: O. Sowatzke
%
% Modified: 12/01/2022
%
% Subject: Class generates RDMs for a radar system using stretch processing
%
classdef radar < handle

    % public properties
    properties

        % center frequency (Hz)
        f_c = 77e9;

        % chirp starting frequency (Hz)
        f_start = 77.2e9;

        % chirp ending frequency (Hz)
        f_stop = 76.8e9;

        % pulse repetition interval (s)
        PRI = 30e-6;

        % Chip pulse width (s)
        P_T = 20e-6;

        % number of chips
        K = 500;

        % ADC sample frequency (Hz)
        f_s = 20e6;

        % light speed (m/s)
        c = 3e8;

        % unambiguous velocity interval select
        % 0 => vua_int = [0 vua)
        % 1 => vua_int = [-vua/2 vua)
        vua_sel = 0;
    end

    % private properties
    properties(Access=protected)

        % Nx1 vector of target ranges (m)
        R;

        % Nx1 vector of target velocities (m/s)
        v;

        % Nx1 vector of targets RCS (m^2)
        RCS;

        % Thermal noise power (W)
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

        % velocity step size corresponding to a doppler bin (m/s)
        function y = vs(self)
            
            % compute width of each doppler bin
            Fd = 1/(self.PRI*self.K);

            % convert to velocity
            y = Fd*self.lambda/2;
        end

        % wavelength (m)
        function y = lambda(self)
            y = self.c/self.f_c;
        end

        % unambiguous range (m)
        function y = Rua(self)

            % Compute time corresponding to unambiguous range. Limiting 
            % factor is ADC rate instead of PRI in stretch processing
            t = self.f_s*self.P_T/self.B;

            % compute unambiguous range
            y = self.c*t/2;
        end

        % unambiguous velocity (m/s)
        function y = vua(self)
            y = self.lambda/(2*self.PRI);
        end

        % size of data cube along fast time axis
        function y = M(self)

            y = round(self.P_T*self.f_s);
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

            % plot resulting rdm
            self.plot_rdm;
        end
    end

    % Private methods
    methods(Access=protected)

        % Function generates the received data cube 
        % [fast time samples x slow time samples] for a system employing
        % stretch processing
        function simulate_received(self)

            % generate thermal noise for each data cube sample
            thermal_noise = normrnd(0,1,self.M,self.K) + ...
                1i.*normrnd(0,1,self.M,self.K);
            thermal_noise = thermal_noise*sqrt(self.thermal_noise_power);
            
            % compute frequency offet due to target range (Hz)
            f_r = 2*self.R/self.c*(self.f_start-self.f_stop)/self.P_T;

            % frequency offset due to doppler shift (Hz)
            f_d = 2*self.v/self.lambda;
            
            % PRI length in samples
            pri_len = round(self.PRI*self.f_s);

            % compute amplitude of each target return
            %   Pr = (Pt*G^2*lambda^2*RCS)/((4*pi)^3*R^4)
            %
            % equation assumes:
            %   (Pt*G^2*lambda^2)/((4*pi)^3) = 1
            %
            amp = ((sqrt(self.RCS))./((self.R).^2));

            % timestep (s)
            dt = 1/self.f_s;

            % time axis for a pulse (s)
            t_pulse=dt*(0:(pri_len-1));

            % time axis for a cpi (s)
            t_cpi=dt*(0:(pri_len*self.K-1));

            % compute return of target due to range effects only
            % repeats each pulse due to indentical transmit signal
            rec1 = repmat(exp(1i*2*pi*f_r.*t_pulse),1,self.K);

            % compute doppler shift of each target return
            % complex exponential operates over full CPI
            % produce doppler shift and range-doppler coupling
            rec2 = exp(1i*2*pi*f_d.*t_cpi);

            % compute overall return using superposition of each target return
            rec3 = sum(amp.*rec1.*rec2,1);

            % compute data cube due to signal returns only
            received = reshape(rec3,pri_len,self.K);
            received = received(1:self.M,:);

            % compute total received data cube (signal + noise)
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
                range_fft_out = circshift(flip(range_fft_out,1),1);
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
            self.range_gate = mod(round(Ra/self.Rs),self.M);

            % normalized frequency fd/PRF
            fnorm = mod(fd*self.PRI,1);

            % doppler shift due to frequency offset
            self.dopp_bin = mod(round(fnorm*self.K),self.K);
        end

        % plot rdm and mark target returns
        function plot_rdm(self)

            % range axis
            raxis = self.Rs*(0:(self.M-1));

            % velocity axis
            % [0, vua) is default unambiguous velocity interval
            vaxis = self.vs*(0:(self.K-1));

            % for [-vua/2, vua/2) unambiguous velocity interval
            if self.vua_sel
                vaxis = vaxis - self.vua/2;
            end

            % generate RDM for plotting
            % [0, vua) is default unambiguous velocity interval
            rdm_plot = 20*log10(abs(self.rdm));

            % for [-vua/2, vua/2) unambiguous velocity interval
            if self.vua_sel
                rdm_plot = fftshift(rdm_plot,2);
            end

            % plot rdm
            figure;
            m = mesh(vaxis,raxis,rdm_plot,'FaceColor','flat');
            view(2);

            % compute measured range and velocity of each target
            % [0, vua) is default unambiguous velocity interval
            vm = self.vs*self.dopp_bin;
            Rm = self.Rs*self.range_gate;

            % for [-vua/2, vua/2) unambiguous velocity interval
            if self.vua_sel
                vm = vm - self.vua*(vm >= (self.vua/2));
            end

            % mark each of the peaks
            for i = 1:length(self.range_gate)
                datatip(m, vm(i), Rm(i));
            end

            xlim([vaxis(1) vaxis(end)])
            ylim([raxis(1) raxis(end)])

            % scale axis
            colorbar;
        
            % label plot
            xlabel('Velocity (m/s)')
            ylabel('Range (m)')
        end
    end
end
