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