% Author: O. Sowatzke
%
% Modified: 11/25/2022
%
% Subject: Function generates and plots the RDM for a radar system using
% stretch processing. Windows and amplitude correction can be configured
% for different scenarios.
%
% Inputs:
%       R                       Nx1 vector of target ranges (m)
%       v                       Nx1 vector of target velocities (m/s)
%       RCS                     Nx1 vector of target RCS (m^2)
%       thermal_noise_power     Normalized thermal noise power expressed in
%                               units of ((4*pi)^3)/(Pt*G^2*lambda^2)
%       en_range_win            Enable FFT window for range FFT
%       en_amp_corr             Enable R^2 Amplitude Compensation
%       en_dopp_win             Enable FFT window for doppler FFT
%
function generate_rdm(R,v,RCS,thermal_noise_power,en_range_win,en_amp_corr,en_dopp_win)
    
    % generate received data cube
    data_cube = simulate_received(R,v,RCS,thermal_noise_power);
    
    % form rdm from received data cube
    rdm = process_received(data_cube,en_range_win,en_amp_corr,en_dopp_win);
    
    % get location of each target in RDM
    [range_gate,dopp_bin] = find_rdm_peaks(R,v);

    % plot rdm
    figure;
    m = mesh(db(rdm),'FaceColor','flat');
    view(2);

    % mark each of the peaks
    for i = 1:length(range_gate)
        datatip(m, dopp_bin(i), range_gate(i));
    end

    % scale axis
    colorbar;

    % label plot
    xlabel('Doppler Bin');
    ylabel('Range Gate');
end