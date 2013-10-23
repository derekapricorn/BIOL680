%Generate a 10-second long white noise signal, sampled at 500Hz. 
%Filter it using a bandpass Butterworth filter between 50 and 100 Hz of order 4. 
%Plot the Welch spectrum, in dB, of the original and filtered signal, 
%using a 512-sample Hanning window. Evaluate the FFT over 2^14 points. 
% set up time axis
Fs = 500;
t0 = 0; t1 = 10;
tvec = t0:1/Fs:t1-(1/Fs);
 
% generate white noise
x = rand(1,length(tvec));
 
% get PSD
nP =  2^14;
wSize = 512;
[Porig,Forig] = pwelch(x, hanning(wSize),wSize/2,nP,Fs);
plot(Forig,10*log10(Porig),'r'); xlabel('Frequency (Hz)');
hold on

% design filter
Fc = Fs/2; %cut-off frequency
W1 = 50/Fc;
W2 = 100/Fc;
[b,a] = butter(4,[W1 W2]);
y = filter(b,a,x);%filter the original signal, not the processed coefficients
%after pwelch
 
% get PSD
[Pfilt,Ffilt] = pwelch(y, hanning(wSize),wSize/2, nP,Fs);
 
% plot the resulting PSDs
%subplot(121)
plot(Ffilt,10*log10(Pfilt),'b');
grid on
%%
Wp = [ 50 100] * 2 / Fs; % passband - between 50 and 100 Hz
Ws = [ 45 105] * 2 / Fs; % stopband
[N,Wn] = buttord( Wp, Ws, 3, 20); % determine filter parameters
[b2,a2] = butter(N,Wn); % builds filter
%%
fvtool(b,a,b2,a2)

%%
Wp = [ 50 100] * 2 / Fs; 
Ws = [ 48 102] * 2 / Fs;
[N,Wn] = cheb1ord( Wp, Ws, 3, 20); 
[b_c1,a_c1] = cheby1(N,0.5,Wn);
fvtool(b2,a2,b_c1,a_c1)

%%
Fs = 500; dt = 1./Fs;
t = [0 10];
tvec = t(1):dt:t(2)-dt;
 
s1 = sin(2*pi*80*tvec+pi/6);
s2 = sin(2*pi*40*tvec);
s = s1 + s2;
 
sf = filter(b_c1,a_c1,s);
 
plot(tvec,s,'k',tvec,sf,'r--'); hold on;
legend({'original','filtered'});
xlim([0 0.2]);

%%
sf = filtfilt(b_c1,a_c1,s);
 
plot(tvec,s,'k',tvec,sf,'r--'); hold on;
legend({'original','filtered'});
xlim([0 0.2]);

%% compare freq responses
Fs = 500; dt = 1./Fs;
t = [0 10];
tvec = t(1):dt:t(2)-dt;
 
x = rand(size(tvec)); % white noise input
[P,F] = pwelch(x,hanning(512),256,2^14,Fs);
 
y1 = filter(b_c1,a_c1,x);
[P1,F1] = pwelch(y1,hanning(512),256,2^14,Fs);
 
y2 = filtfilt(b_c1,a_c1,x);
[P2,F2] = pwelch(y2,hanning(512),256,2^14,Fs);
 
plot(F,10*log10(P),F,10*log10(P1),F,10*log10(P2));
legend({'original','filter','filtfilt'});

%%
[b,a] = butter(10, [59 61] * 2 / Fs, 'stop');
fvtool(b,a);

%% test the notch filter on white noise
[z,p,k] = butter(10, [59 61] * 2 / Fs, 'stop'); % note, we ask for 3 outputs instead of 2
[sos,g] = zp2sos(z,p,k); % convert to SOS format
h = dfilt.df2sos(sos,g); % create filter object
%fvtool(h);
y3 = filtfilt(b,a,x);
[P3,F3] = pwelch(y3,hanning(512),256,2^14,Fs);
plot(F,10*log10(P3));

%%
%% cd to R016-2012-10-08 folder first
cd('D:\data\promoted\R016\R016-2012-10-08');
csc = LoadCSC('R016-2012-10-08-CSC02b.ncs');
cscR = Restrict(csc,1270,1272);
plot(cscR)

%%
x = Data(cscR);
tvec = Range(cscR);
 
Fs = 2000;
Wp = [ 180 220] * 2 / Fs;
Ws = [ 178 222] * 2 / Fs;
[N,Wn] = cheb1ord( Wp, Ws, 3, 20); % determine filter parameters
[b_c1,a_c1] = cheby1(N,0.5,Wn); % builds filter
 
%fvtool(b_c1,a_c1); % remember to check your filter!
 
y = filtfilt(b_c1,a_c1,x);
plot(tvec,x,'b',tvec,y,'r');
%%
chew_power = y.^2;
chew_power_filtered = medfilt1(chew_power,101); % filter window is specified in samples, so this is ~50ms
plot(tvec,x,'b',tvec,chew_power_filtered,'r');


%%

function [evt, t, pwr] = detectSWR(csc,varargin)
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

cd('D:\data\promoted\R042\R042-2013-08-18');
csc = LoadCSC('R042-2013-08-18-CSC03a.ncs');

ripple_band = [140 180]; 
threshold = 5; 
evt = [];
t = [];
pwr = [];


% extract varargins
extract_varargin;
% load signal
x = Data(csc);
tvec = Range(csc);
% filter in the ripple band
Fs = 2000;
%use cheby filter because it has a sharper rolloff good for avoiding 
%confusion with chewing artifacts. this will give you +/-5db distortion 
Wp = ripple_band * 2 /Fs;
Ws = [ripple_band(1)-2, ripple_band(2)+2] * 2 / Fs;
[N,Wn] = cheb1ord(Wp, Ws, 3, 20);
[b_c,a_c] = cheby1(N,0.1,Wn);
fvtool(b_c,a_c);
y = filtfilt(b_c,a_c,x);
plot(tvec,x,'b',tvec,y,'r');

% convert to power envelope
power_env = y.^2;
% convert to z-score (SDs from the mean; use nanmean() and nanstd())
chew_power_filtered = medfilt1(power_env,404);% filter window is specified in samples, so this is ~200ms
plot(tvec,x,'b',tvec,chew_power_filtered,'r');
% find times when above threshold
 
% find crossings from below to above threshold and vice versa (can use diff())
 
% get center time and power, return












