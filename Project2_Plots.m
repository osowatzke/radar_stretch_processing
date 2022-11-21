% Author: O. Sowatzke
%
% Updated: 12/01/2022
%
% Subject : Function generates helper plots used in report

%% Instantaneous Frequencies for Transmit and Receive Signals
% array of frequencies and times for plotting (arbitrarily chosen)
f = 0:0.01:1;
t = 0:0.01:1;

% time offset
t0 = 0.4;

% frequency offset
fd = 0.2;

% create figure
figure(1)
clf;

% plot instantaneous frequency of transmit signal
plot(t,f,'LineWidth',1.5);
hold on;

% compute times of received signals (delayed by t0)
t = t + t0;

% plot instantaneous frequency of receive signal (static target)
plot(t,f,'LineWidth',1.5);

% apply doppler shift for moving target
f = f + fd;

% plot instantaneous frequency of receive signal (moving target)
plot(t,f,'LineWidth',1.5);

% label time axis
xticks([0, t0, 1, 1+t0])
xticklabels({'0','t_0','\tau','\tau+t_0'});

% label frequency axis
yticks([0, fd, 1, 1+fd])
yticklabels({'0','F_d','\beta','\beta + F_d'});

% add labels to plot
grid on;
xlabel('Time');
ylabel('Frequency')
legend({'Transmit Signal','Receive Signal (Static Target)','Receive Signal (Moving Target)'},'location','northwest');
title('Instantaneous Frequency of Transmit and Receive Signals')

%% Beat Frequencies for Target Returns
% create figure
figure(2)
clf;

% array of times
t = 0:0.01:1;

% zero beat frequency for transmit signal (mixes down to DC)
f = zeros(1,length(t));

% plot beat frequency of transmit signal
plot(t,f,'LineWidth',1.5);
hold on;

% add time offset for received signal
t = t + t0;

% compute beat frequency of received signal (negative for upchirp)
f = -t0*ones(1,length(t));

% plot beat frequency of received signal (static target)
plot(t,f,'LineWidth',1.5);

% apply frequency offset for moving target
f = f + fd;

% plot beat frequency of received signal (static target)
plot(t,f,'LineWidth',1.5);

% scale frequency axis
ylim([-1 1]);

% label time axis
xticks([0, t0, 1, 1+t0]);
xticklabels({'0','t_0','\tau','\tau+t_0'});

% label frequency axis
yticks([-t0 -t0+0.2 0]);
yticklabels({'F_r','F_r + F_d','0'});

% add labels to plot
grid on;
xlabel('Time');
ylabel('Frequency')
title('Beat Frequency of Each Target Return')
legend({'Transmit Signal','Receive Signal (Static Target)','Receive Signal (Moving Target)'},'location','northwest');

%% Beat Frequencies of Target Returns after Group Delay Correction
% create new figure
figure(3)
clf;

% array of times
t = 0:0.01:1;

% zero beat frequency for transmit signal (mixes down to DC)
f = zeros(1,length(t));

% plot beat frequency of transmit signal
plot(t,f,'LineWidth',1.5);
hold on;

% compute beat frequency of received signal (negative for upchirp)
f = -t0*ones(1,length(t));

% plot beat frequency of received signal (static target)
plot(t,f,'LineWidth',1.5);

% apply frequency offset for moving target
f = f + fd;

% plot beat frequency of received signal (static target)
plot(t,f,'LineWidth',1.5);

% scale frequency axis
ylim([-1 1]);

% scale time axis
xlim([0 1+t0]);

% label time axis
xticks([0, t0, 1, 1+t0]);
xticklabels({'0','t_0','\tau','\tau+t_0'});

% label frequency axis
yticks([-t0 -t0+0.2 0]);
yticklabels({'F_r','F_r + F_d','0'});

% add labels to plot
xlabel('Time');
ylabel('Frequency')
title('Beat Frequencies after Analog Filtering')
legend({'Transmit Signal','Receive Signal (Static Target)','Receive Signal (Moving Target)'},'location','northwest');
grid on;