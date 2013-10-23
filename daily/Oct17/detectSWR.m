function [evt, t, pwr] = detectSWR(csc,varargin)
%%
clear all
close all
% function evt = detectSWR(csc,varargin)
%
% detect putative sharp wave-ripple events in csc input signal
%
% INPUTS
% csc: [1 x 1 mytsd]
%
% varargins (with defaults):
% ripple_band = [140 180]; % frequency band to use
% threshold = 5; % number of SDs above mean to use for detection
%
% OUTPUTS
% evt: [1 x 1 struct] with fields
% .t: [1 x nEvents double] times (in s) of events
% .pwr: [1 x nEvents double] power of events (in SDs above mean)
%%
cd('D:\data\promoted\R042\R042-2013-08-18');
csc = LoadCSC('R042-2013-08-18-CSC03a.ncs');

%restrict the csc to 6000s to 6300s
cscR = Restrict(csc, 6000, 6900);

ripple_band = [140 180]; 
threshold = 5; 
evt = [];
t = [];
pwr = [];

%%
% extract varargins
extract_varargin;
% load signal

Fs = 2000;%specify sampling frequency

%downsample the data
dsf = 10;
x = decimate(Data(cscR),dsf);
tvec = downsample(Range(cscR),dsf);
Fs = Fs/dsf;
% filter in the ripple band
%use cheby filter because it has a sharper rolloff good for avoiding 
%confusion with chewing artifacts. this will give you +/-5db distortion 
Wp = ripple_band * 2 /(Fs*dsf);
Ws = [ripple_band(1)-2, ripple_band(2)+2] * 2 / (Fs*dsf);
[N,Wn] = cheb1ord(Wp, Ws, 3, 20);
[b_c,a_c] = cheby1(N,0.1,Wn);%suppress the ripples to less than 0.1dB
%fvtool(b_c,a_c);
y = filtfilt(b_c,a_c,x);
%plot(tvec,x,'b',tvec,y,'r');
% convert to power envelope
power_env = y.^2;
figure
subplot(3,1,1);
plot(tvec,x,'b',tvec,y,'r');
legend('original','filtered');
ylabel('LFP');
%%
% convert to z-score (SDs from the mean; use nanmean() and nanstd())
power_filtered = medfilt1(power_env,401);% filter window is specified in samples, so this is ~200ms
%plot(tvec,x,'b',tvec,power_filtered,'r');

mu = nanmean(power_filtered);
theta = nanstd(power_filtered);
z = (power_filtered - mu)/theta;

% find times when above threshold
temp = z>threshold;
% find crossings from below to above threshold and vice versa (can use diff())
temp_diff = diff(temp);
zero_xing_start = find(temp_diff == 1);
zero_xing_end = find(temp_diff == -1);

% get center time and power, return
centre = round((zero_xing_start+zero_xing_end)/2);
%t = tvec(centre)*10^-6;
t = tvec(centre);
pwr = power_env(centre);
evt = struct('time', t, 'power', pwr);
subplot(3,1,2)
plot(tvec,z,[6000 6600],[threshold threshold]);
legend('z-score','threshold');
ylabel('z-score plot')
subplot(3,1,3)
stem(t,pwr); legend('power');
ylabel('power of events');
end
