% Author: O. Sowatzke
%
% Modified: 11/25/2022
%
% Subject: Function computes the expected range gate and doppler bin for a
% vectors of targets with range R and velocity V
%
% Inputs:
%           R               Nx1 vector of target ranges (m)
%           v               Nx1 vector of target velocities (m/s)
%
% Outputs:
%           range_gate      Nx1 vector of range gates
%           dopp_bin        Nx1 vector of doppler bins
%
function [range_gate, dopp_bin] = find_rdm_peaks(R, v)
    
    % compute doppler shift
    fd = 2*v*CONSTANTS.f_c/CONSTANTS.c;

    % compute range shift due to target velocity
    % upchirp   => deltaR < 0
    % downchirp => deltaR > 0
    deltaR = CONSTANTS.c*fd*CONSTANTS.P_T/(2*CONSTANTS.B);
    if CONSTANTS.f_start < CONSTANTS.f_stop
        deltaR = -deltaR;
    end
    
    % compute adjusted range
    R_adj = R + deltaR;

    % compute ambiguous range
    Ra = mod(R_adj, CONSTANTS.Rua);

    % quantize to a range gate
    % add one to produce MATLAB indexing
    range_gate = mod(round(Ra/CONSTANTS.Rs),CONSTANTS.M) + 1;

    % normalized frequency fd/PRF
    fnorm = mod(fd*CONSTANTS.PRI,1);

    % compute doppler bin
    dopp_bin = fnorm*CONSTANTS.K;

    % doppler shift due to frequency offset
    % add one to produce MATLAB indexing
    dopp_bin = mod(round(dopp_bin),CONSTANTS.K) + 1;
end