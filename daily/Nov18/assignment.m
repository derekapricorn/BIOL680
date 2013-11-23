%%
%*************assignment***********%

%% loading
cd 'D:\data\promoted\R042\R042-2013-08-18'
 
fc = FindFiles('*.t');
S = LoadSpikes(fc);

% for each of the two neurons, restrict the data to [3200 5650] (the time interval when the rat was running on the track)
cell1_id = 5; cell2_id = 42;
s1 = Restrict(S{cell1_id},3200,5650); % restrict to on-track times only
s2 = Restrict(S{cell2_id},3200,5650);

% compute the spike density function for each, making sure that your tvec runs from 3200 to 5650 also, and that you have a 50ms SD for the Gaussian convolution kernel
dt = 0.001;
t = [3200 5650]; % time interval (length) of spike train to generate
tvec = t(1):dt:t(2);
tvec = tvec(1:end-1);
binsize = 0.001; % in seconds, so everything else should be seconds too
tbin_edges = t(1):binsize:t(2);
tbin_centers = tbin_edges(1:end-1)+binsize/2;

gauss_window = 1./binsize; % 1 second window
gauss_SD = 0.05./binsize; % 0.05 seconds (50ms) SD
gk = gausskernel(gauss_window,gauss_SD); gk = gk./binsize; % normalize by binsize

spk_t1 = data(s1);
spk_t2 = data(s2);
spk_count1 = histc(spk_t1,tbin_edges);
spk_count1 = spk_count1(1:end-1);
spk_count2 = histc(spk_t2,tbin_edges);
spk_count2 = spk_count2(1:end-1);

gau_sdf1 = conv2(spk_count1,gk,'same'); % convolve with gaussian window
gau_sdf2 = conv2(spk_count2,gk,'same'); % convolve with gaussian window

% to use these SDFs to generate Poisson spike trains, convert the firing rates 
%given by the SDF to a probability of emitting a spike in a given bin. (As you did above for a 0.47 Hz constant firing rate.)
 
pspike1 = gau_sdf1 * 10^-3; % probability of generating a spike in bin
pspike2 = gau_sdf2 * 10^-3; % probability of generating a spike in bin
%rng default; % reset random number generator to reproducible state
spk_poiss = rand(size(tvec)); % random numbers between 0 and 1
spk_poiss_idx1 = find(spk_poiss < pspike1'); % index of bins with spike
spk_poiss_t1 = tvec(spk_poiss_idx1)'; % use idxs to get corresponding spike time
spk_poiss_idx2 = find(spk_poiss < pspike2'); % index of bins with spike
spk_poiss_t2 = tvec(spk_poiss_idx2)'; % use idxs to get corresponding spike time

% generate Poisson spike trains, making sure to use the same tvec
 line([spk_poiss_t1 spk_poiss_t1],[-1 -0.5],'Color',[0 0 0]); % note, plots all spikes in one command
 title('spike train for cell 1');
 figure
 line([spk_poiss_t2 spk_poiss_t2],[-1 -0.5],'Color',[0 0 0]); % note, plots all spikes in one command
 title('spike train for cell 2');


% convert Poisson spike trains to ts objects and compute the ccf
ts_new1 = ts(spk_poiss_t1);
ts_new2 = ts(spk_poiss_t2);

[xcorr,xbin] = ccf(ts_new1,ts_new2,binsize,1);

plot(xbin,xcorr);
set(gca,'FontSize',20); xlabel('lag (s)'); ylabel('xcorr');
title(sprintf('%d-%d',cell1_id,cell2_id));