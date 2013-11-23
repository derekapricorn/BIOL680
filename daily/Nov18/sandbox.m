%% loading
cd 'D:\data\promoted\R042\R042-2013-08-18'
 
fc = FindFiles('*.t');
S = LoadSpikes(fc);
 
%% plot
iC = 47;
t = [5801 5801.7];
 
spk_t = Data(Restrict(S{iC},t(1),t(2))); % get spike times
 
line([spk_t spk_t],[-1 -0.5],'Color',[0 0 0]); % note, plots all spikes in one command

%%
binsize = 0.1;
tbin_edges = t(1):binsize:t(2); % vector of time bin edges (for histogram)
tbin_centers = tbin_edges(1:end-1)+binsize/2; % vector of time bin centers (for plotting)
 
spk_count = histc(spk_t,tbin_edges); % get spike counts for each bin
spk_count = spk_count(1:end-1); % ignore spikes falling exactly on edge of last bin.

%%
hold on;
h = bar(tbin_centers,spk_count);
set(h,'BarWidth',1,'EdgeColor','none','FaceColor',[1 0 0]); % reformat bar appearance
 
yt = get(gca,'YTick'); ytl = get(gca,'YTickLabel');
set(gca,'YLim',[-1.5 10],'YTick',yt(2:end),'YTickLabel',ytl(2:end,:)); xlabel('time (s)');

%%
binsize = 0.001; % select a small bin size for good time resolution
tbin_edges = t(1):binsize:t(2);
tbin_centers = tbin_edges(1:end-1)+binsize/2;
 
spk_count = histc(spk_t,tbin_edges);
spk_count = spk_count(1:end-1);
 
rw_sdf = conv2(spk_count,rectwin(50),'same'); % convolve with rectangular window
plot(tbin_centers,rw_sdf,'b');
 
gau_sdf = conv2(spk_count,gausswin(50),'same'); % convolve with gaussian window
plot(tbin_centers,gau_sdf,'g');

%%
binsize = 0.001; % in seconds, so everything else should be seconds too
gauss_window = 1./binsize; % 1 second window
gauss_SD = 0.02./binsize; % 0.02 seconds (20ms) SD
gk = gausskernel(gauss_window,gauss_SD); gk = gk./binsize; % normalize by binsize
gau_sdf = conv2(spk_count,gk,'same'); % convolve with gaussian window
plot(tbin_centers,gau_sdf,'g');

%% real isi
iC = 47;
spk_t = Data(S{iC}); % spike times
isi = diff(spk_t); % interspike intervals
 
dt = 0.001; % in s, because spike times are in s
isi_edges = 0:dt:0.25; % bin edges for ISI histogram
isi_centers = isi_edges(1:end-1)+dt/2; % for plotting
 
isih = histc(isi,isi_edges);
 
bar(isi_centers,isih(1:end-1)); % remember to ignore the last bin of histc() output
set(gca,'FontSize',20,'XLim',[0 0.25]); xlabel('ISI (s)'); ylabel('count'); grid on;


%% return plot
isi_next = isi(1:(end-1));
isi_next = [0;isi_next];
plot(isi,isi_next,'.');xlim([0 0.25]);ylim([0 0.25]);
%% synthetic plot
dt = 0.001;
t = [0 4500]; % time interval (length) of spike train to generate
tvec = t(1):dt:t(2);

pspike = 4.7*10^-4; % probability of generating a spike in bin
%rng default; % reset random number generator to reproducible state
spk_poiss = rand(size(tvec)); % random numbers between 0 and 1
spk_poiss_idx = find(spk_poiss < pspike); % index of bins with spike
spk_poiss_t = tvec(spk_poiss_idx)'; % use idxs to get corresponding spike time
 
line([spk_poiss_t spk_poiss_t],[-1 -0.5],'Color',[0 0 0]); % note, plots all spikes in one command
axis([0 0.1 -1.5 5]); %set(gca,'YTick','');
%%
isi = diff(spk_poiss_t);
dt = 0.001; % in s, because spike times are in s
isi_edges = 0:dt:0.25; % bin edges for ISI histogram
isi_centers = isi_edges(1:end-1)+dt/2; % for plotting
 
isih = histc(isi,isi_edges);
figure
bar(isi_centers,isih(1:end-1)); % remember to ignore the last bin of histc() output
set(gca,'FontSize',20,'XLim',[0 0.25]); xlabel('ISI (s)'); ylabel('count'); grid on;

%% synthetic plot varying pspike
dt = 0.001;
t = [0 10]; % time interval (length) of spike train to generate
tvec = t(1):dt:t(2);

pspike = zeros(1,length(tvec)); % probability of generating a spike in bin
for ipspike = 2:2:length(pspike)
    pspike(ipspike) = 0.5;
end

%rng default; % reset random number generator to reproducible state
spk_poiss = rand(size(tvec)); % random numbers between 0 and 1
spk_poiss_idx = find(spk_poiss < pspike); % index of bins with spike
spk_poiss_t = tvec(spk_poiss_idx)'; % use idxs to get corresponding spike time
 
line([spk_poiss_t spk_poiss_t],[-1 -0.5],'Color',[0 0 0]); % note, plots all spikes in one command
axis([0 0.1 -1.5 5]); %set(gca,'YTick','');

%% comparing the autocorrelation of the real ici and poisson ici
binsize = 0.001; % select a small bin size for good time resolution

dt = 0.001;
t = [0 4500]; % time interval (length) of spike train to generate
tvec = t(1):dt:t(2);

pspike = 4.7*10^-4; % probability of generating a spike in bin
%rng default; % reset random number generator to reproducible state
spk_poiss = rand(size(tvec)); % random numbers between 0 and 1
spk_poiss_idx = find(spk_poiss < pspike); % index of bins with spike
spk_poiss_t = tvec(spk_poiss_idx)'; % use idxs to get corresponding spike time
line([spk_poiss_t spk_poiss_t],[-1 -0.5],'Color',[0 0 0]); % note, plots all spikes in one command
axis([0 0.1 -1.5 5]); %set(gca,'YTick','');

% ts_new = ts(spk_poiss_t);
% [ac,xbin] = acf(ts_new,binsize, 1);
% plot(xbin,ac);
% hold on
% [acorr,xbin] = acf(S{47},0.01,1);
% plot(xbin, acorr,'g')

%%
binsize = 0.01; % select a small bin size for good time resolution

dt = 0.001;
t = [0 4500]; % time interval (length) of spike train to generate
tvec = t(1):dt:t(2);

pspike = 4.7*10^-4; % probability of generating a spike in bin
%rng default; % reset random number generator to reproducible state
spk_poiss = rand(size(tvec)); % random numbers between 0 and 1
spk_poiss_idx = find(spk_poiss < pspike); % index of bins with spike
spk_poiss_t = tvec(spk_poiss_idx)'; % use idxs to get corresponding spike time
ts_new = ts(spk_poiss_t);
[cc,xbin] = ccf(ts_new, ts_new, binsize, 1);
plot(xbin, cc);

%%
cell1_id = 5; cell2_id = 42;
 
s1 = Restrict(S{cell1_id},3200,5650); % restrict to on-track times only
s2 = Restrict(S{cell2_id},3200,5650);
 
[xcorr,xbin] = ccf(s1,s2,0.01,1);
 
plot(xbin,xcorr);
set(gca,'FontSize',20); xlabel('lag (s)'); ylabel('xcorr');
title(sprintf('%d-%d',iX,iY));

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