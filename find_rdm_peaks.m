function [range_gate, dopp_bin] = find_rdm_peaks(R, v)
    
    % compute doppler shift
    fd = 2*v*CONSTANTS.f_c/CONSTANTS.c;

    % compute range shift
    % fd = B/tau*delta_t => delta_t = fd*tau/B
    % deltaR = c*delta_t/2 => delta_R =2* c*fd*tau/(2*B)
    deltaR = -CONSTANTS.c*fd*CONSTANTS.P_T/(2*CONSTANTS.B);

    % adjusted range
    R_adj = R + deltaR;

    % compute ambiguous range
    Ra = mod(R_adj, CONSTANTS.Rua);

    % quantize to a range gate
    % add 1 for MATLAB indexing
    range_gate = mod(round(-Ra/CONSTANTS.Rs),CONSTANTS.M) + 1;

    % normalized frequency fd/PRF
    fnorm = mod(fd*CONSTANTS.PRI,1);

    % compute doppler bin
    dopp_bin = fnorm*CONSTANTS.K;

    % doppler shift due to frequency offset
    dopp_bin = mod(round(dopp_bin),CONSTANTS.K)+ 1;
end