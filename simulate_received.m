function data_cube = simulate_received(R, v, RCS, thermal_noise_power)
    
    % timestep (s)
    dt=1/CONSTANTS.f_s;

    % generate thermal noise
    thermal_noise=normrnd(0,1,CONSTANTS.M,CONSTANTS.K) + ...
        1i.*normrnd(0,1,CONSTANTS.M,CONSTANTS.K);
    thermal_noise=thermal_noise*sqrt(thermal_noise_power);

    % compute frequency offet due to target range
    f_r=2*R/CONSTANTS.c*CONSTANTS.B/CONSTANTS.P_T;

    % frequency offset due to doppler shift
    f_d=2*v*CONSTANTS.f_c/CONSTANTS.c;

    % array of times in (s)
    t=0:dt:dt*(CONSTANTS.PRI*CONSTANTS.f_s*CONSTANTS.K-1);

    % generate received signal
    %
    % Pr = Pt*G^2*lambda^2*sigma/((4*pi)^3*R^4)
    % Vr = (sqrt(Pt/(4*pi)^3)*G*lambda)*sqrt(sigma)/R^2
    %
    % Since only the sqrt(sigma)/R^2 term is accounted for, noise standard
    % deviation is given in units of (sqrt(Pt/(4*pi)^3)*G*lambda)
    %
    RCS=sqrt(RCS);
    received=((RCS)./R.^2)'*exp(1i.*2.*pi.*(-f_r+f_d)*t);

    % generate datacube
    reshape_received=reshape(received,CONSTANTS.PRI*CONSTANTS.f_s,CONSTANTS.K);
    data_cube=reshape_received(1:CONSTANTS.M,:);
    data_cube=data_cube+thermal_noise;
end