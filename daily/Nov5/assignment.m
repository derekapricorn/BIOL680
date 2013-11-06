%%
%Assignment 7 

%% load and restrict the data
cd('D:\data\promoted\R016\R016-2012-10-03');
fname = 'R016-2012-10-03-CSC04a.Ncs';
csc = LoadCSC(fname);
 
cscR = Restrict(csc,2700,3300); % risk session only
Fs = 2000;
x = Data(cscR);
tvec = Range(cscR);

%% define frequency ranges of interest
gamma = [45 65];
delta = [3 4];
%% design filters for frequency ranges
Wp_d =  delta * 2 /Fs;
Ws_d = [delta(1)-2, delta(2)+ 2] * 2 / Fs;
[N,Wn] = cheb1ord(Wp_d, Ws_d, 3, 20);
[b_d,a_d] = cheby1(N,0.1,Wn); %suppress the ripples to less than 0.1dB
%fvtool(b_d,a_d); %check the filter shape

Wp_g =  gamma * 2 /Fs;
Ws_g = [gamma(1)-2, gamma(2)+ 2] * 2 / Fs;
[N_g,Wn_g] = cheb1ord(Wp_g, Ws_g, 5, 20);
[b_g,a_g] = cheby1(N_g,0.1,Wn_g); %suppress the ripples to less than 0.1dB
%fvtool(b_g,a_g); %checke the filter shape
%% filter the data (remember to use filtfilt!)
y = filtfilt(b_d,a_d,x);
y_g = filtfilt(b_g,a_g,x);
%% extract delta phase and low gamma power
threshold = 5;
power_env = y.^2;
power_filtered = medfilt1(power_env,101);
mu = nanmean(power_filtered);
theta = nanstd(power_filtered);
z = (power_filtered - mu)/theta;

% find times when above threshold
temp = z>threshold;
% find crossings from below to above threshold and vice versa (can use diff())
temp_diff = diff(temp);
zero_xing_start = find(temp_diff == 1);
zero_xing_end = find(temp_diff == -1);

% get center time and delta phase
centre = round((zero_xing_start+zero_xing_end)/2);
t = tvec(centre);
sgnl= y(centre);
std = z(centre);
h = hilbert(sgnl);
phi = angle(h); %extract delta phase

%extract low gamma power
low_gamma_pwr = y_g.^2;
low_gamma_pwr = low_gamma_pwr(centre);

%% use averageXbyYbin to plot relationship (ideally with standard deviations)
phi_edges = -pi:pi/8:pi;
pow_bin = averageXbyYbin(low_gamma_pwr,phi,phi_edges);
 
pow_bin(end-1) = pow_bin(end-1)+pow_bin(end); % add counts on last edge to preceding bin
pow_bin = pow_bin(1:end-1); % trim
 
phi_centers = phi_edges(1:end-1)+pi/16; % convert edges to centers
plot(phi_centers,pow_bin);
% hold on
% plot(phi_centers, std);
% legend('low gamma power','std');
title('phase-power correlation');
xlabel('phase');
ylabel('power bins');
