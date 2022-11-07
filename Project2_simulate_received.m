clear all;
close all;

f_c=77e9; % center frequency
f_start=76.8e9;% chirp starting frequency
f_stop=77.2e9;% chirp end frequency
PRI=30e-6;%Pulse Repetition Interval
P_T=20e-6; %chirp pulse width
M=400; % fast time samples
K=500; % number of chirps
f_s=20e6; % ADC sample frequency
c=3e8; % light speed
dt=1/f_s;
theormal_noise_power=0; % thermal noise power

R=[10, 100]'; % target range
v=[15, 10]'; % target velocity
RCS=[1, 10000]'; % target RCS


theormal_noise=normrnd(0,1,M,K)+1i.*normrnd(0,1,M,K);
theormal_noise=theormal_noise*sqrt(theormal_noise_power);

f_r=2*R/c*(f_stop-f_start)/P_T;
f_d=2*v*f_c/c;

t=0:dt:dt*(PRI*f_s*K-1);

RCS=sqrt(RCS);
received=((RCS)./R.^2)'*exp(1i.*2.*pi.*(f_r+f_d)*t);

reshape_received=reshape(received,PRI*f_s,K);
data_cube=reshape_received(1:round(P_T*f_s),:);
data_cube=data_cube+theormal_noise;
