% Author: O. Sowatzke
%
% Updated: 11/25/2022
% 
% Subject: Function generates the received data cube 
% [fast time samples x slow time samples] for a system employing
% stretch processing
%
% Inputs:
%       R                       Nx1 vector of target ranges (m)
%       v                       Nx1 vector of target velocities (m/s)
%       RCS                     Nx1 vector of targets RCS (m^2)
%       thermal_noise_power     Normalized thermal noise power expressed in
%                               units of ((4*pi)^3)/(Pt*G^2*lambda^2)
%
% Outputs:
%       data_cube               Received data cube [fast time samples x
%                               slow time samples]
%
function data_cube = simulate_received(R, v, RCS, thermal_noise_power)
    
    % timestep (s)
    dt=1/CONSTANTS.f_s;

    % generate thermal noise for each data cube sample
    thermal_noise=normrnd(0,1,CONSTANTS.M,CONSTANTS.K) + ...
        1i.*normrnd(0,1,CONSTANTS.M,CONSTANTS.K);

    thermal_noise=thermal_noise*sqrt(thermal_noise_power);

    % compute frequency offet due to target range
    f_r=2*R/CONSTANTS.c*(CONSTANTS.f_start - CONSTANTS.f_stop)/CONSTANTS.P_T;

    % frequency offset due to doppler shift
    f_d=2*v*CONSTANTS.f_c/CONSTANTS.c;

    % array of times in (s)
    t=0:dt:dt*(CONSTANTS.PRI*CONSTANTS.f_s*CONSTANTS.K-1);

    % generate received signal
    RCS=sqrt(RCS);
    received=((RCS)./R.^2)'*exp(1i.*2.*pi.*(f_r + f_d)*t);

    % generate datacube
    reshape_received=reshape(received,CONSTANTS.PRI*CONSTANTS.f_s,CONSTANTS.K);
    data_cube=reshape_received(1:CONSTANTS.M,:);
    data_cube=data_cube+thermal_noise;
end