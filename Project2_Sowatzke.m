% Author: O. Sowatzke
% 
% Date: 12/01/2022
%
% Subject: Script generates input data cube for a stretch processor given 
% sets of input parameters. Data cube is processed to form RDM. Resulting
% RDM is then plotted.
%

% close any open figures
close all;

%% Problem 1
% input parameters
R = [10;100];
v = [15;10];
RCS = [1;10000];
thermal_noise_power = 0;
en_range_win = 0;
en_amp_corr = 0;
en_dopp_win = 0;

% create radar object with default properties
r = radar;

% generate rdm
r.generate_rdm(R,v,RCS,thermal_noise_power,en_range_win,...
    en_amp_corr,en_dopp_win);

%% Problem 2a)
% input parameters
R = [10;10];
v = [15;14];
RCS = [100;1];
thermal_noise_power = 1e-4;
en_range_win = 0;
en_amp_corr = 1;
en_dopp_win = 0;

% create radar object with default properties
r = radar;

% generate rdm
r.generate_rdm(R,v,RCS,thermal_noise_power,en_range_win,...
    en_amp_corr,en_dopp_win);

% create zoomed-in figure
ax1 = gca;
f = figure;
copyobj(ax1,f);
xlim([12 17]);
ylim([6 14]);

%% Problem 2b)
% input parameters
R = [10;10];
v = [15;14];
RCS = [100;1];
thermal_noise_power = 1e-4;
en_range_win = 1;
en_amp_corr = 1;
en_dopp_win = 1;

% create radar object with default properties
r = radar;

% generate rdm
r.generate_rdm(R,v,RCS,thermal_noise_power,en_range_win,...
    en_amp_corr,en_dopp_win);

% create zoomed-in figure
ax1 = gca;
f = figure;
copyobj(ax1,f);
xlim([12 17]);
ylim([6 14]);

%% Problem 3b)
% input parameters
R = 200;
v = 80;
RCS = 10000;
thermal_noise_power = 1e-8;
en_range_win = 0;
en_amp_corr = 0;
en_dopp_win = 0;

% create radar object with default properties
r = radar;

% generate rdm
r.generate_rdm(R,v,RCS,thermal_noise_power,en_range_win,...
    en_amp_corr,en_dopp_win);

%% Problem 3c)
% input parameters
R = 200;
v = 80;
RCS = 10000;
thermal_noise_power = 1e-8;
en_range_win = 0;
en_amp_corr = 0;
en_dopp_win = 0;

% create radar object with modified configuration
% increased ADC rate increases size of range window
% increased PRF increases unambiguous velocity
r = radar(...
    'K',500,...
    'f_s',80e6,...
    'PRI',20e-6,...
    'P_T',10e-6);

% generate rdm
r.generate_rdm(R,v,RCS,thermal_noise_power,en_range_win,...
    en_amp_corr,en_dopp_win);

% create radar object with modified configuration
% use for 2nd CPI in multiple CPI solution
r = radar(...
    'K',500,...
    'f_s',15e6,...
    'PRI',40e-6,...
    'P_T',20e-6);

% generate rdm
r.generate_rdm(R,v,RCS,thermal_noise_power,en_range_win,...
    en_amp_corr,en_dopp_win);