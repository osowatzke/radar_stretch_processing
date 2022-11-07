% close any open figures
close all;

%% Problem 1
% parameters
R = [10;100];
v = [15;10];
RCS = [1;10000];
thermal_noise_power = 0;
en_range_win = 0;
en_amp_corr = 0;
en_dopp_win = 0;

% generate rdm
generate_rdm(R,v,RCS,thermal_noise_power,en_range_win,en_amp_corr,en_dopp_win);

%% Problem 2a)
% parameters
R = [10;10];
v = [15;14];
RCS = [100;1];
thermal_noise_power = 1e-4;
en_range_win = 0;
en_amp_corr = 1;
en_dopp_win = 0;

% generate rdm
generate_rdm(R,v,RCS,thermal_noise_power,en_range_win,en_amp_corr,en_dopp_win);

%% Problem 2b)
% parameters
R = [10;10];
v = [15;14];
RCS = [100;1];
thermal_noise_power = 1e-4;
en_range_win = 1;
en_amp_corr = 1;
en_dopp_win = 1;

% generate rdm
generate_rdm(R,v,RCS,thermal_noise_power,en_range_win,en_amp_corr,en_dopp_win);

%% Problem 3)
% parameters
R = 200;
v = 80;
RCS = 10000;
thermal_noise_power = 1e-6;
en_range_win = 0;
en_amp_corr = 0;
en_dopp_win = 0;

% generate rdm
generate_rdm(R,v,RCS,thermal_noise_power,en_range_win,en_amp_corr,en_dopp_win);