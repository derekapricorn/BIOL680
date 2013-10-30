function eventLFPplot(csc,event_times,varargin)
%
% INPUTS
%
% csc: [1 x 1] mytsd, LFP signal to be plotted
% event_times: [nEvents x 1] double with event times to align LFP on
%
% varargins (with defaults):
%
% t_window: [2 x 1] double indicating time window to use, e.g. [-1 3] for 1 second before to 3 seconds after event times

% load the data, for testing purpose
% remember to cd to the correct folder here, may need to get this file from the lab database
% cd('D:\data\promoted\R016\R016-2012-10-03');
% fname = 'R016-2012-10-03-CSC04a.Ncs';
% csc = MyLoadCSC(fname);
%get sampling frequency 


t_window = [-1 3];
dsf = 2;
color = 'b';
hdr = getCSCHeader(csc);
Fs = (hdr.SamplingFrequency)./dsf;
%filt_coeff = [];
%filt_coeff = [1 10 100]; %by default, create a butt passband filter with passband[0 100]
extract_varargin


 %event_times = [2013,2088,2100];
 %t_window = [-1 3];

%error checking 
max_t = max(event_times);
if(max_t > max(Range(csc)))
    exception = MException('VerifyOutput:OutOfBounds', ...
       'Results are outside the allowable limits');
    throw(exception);
end


figure
hold on

% for each event:
for ii = 1:length(event_times)
   
% *extract the corresponding piece of LFP
csc_R = Restrict(csc, event_times(ii)+t_window(1), event_times(ii)+t_window(2));
    
% *replace the original time axis with a new one based on the time window asked for
tvec = downsample(Range(csc_R),dsf);
data = decimate(Data(csc_R),dsf);
Fs_new = Fs./dsf;

% filter if necessary; couldn't get the filtering work
% if(~isempty(filt_coeff))
%     data = filt(data,filt_coeff(1),[filt_coeff(2),filt_coeff(3)], Fs_new);
% end

% *optionally, rescale the LFP
lfp_minmax = 25;  % range of LFP plotting
tvec0 = tvec - tvec(1)+ t_window(1); % align LFP 
data = rescale(data,-lfp_minmax,lfp_minmax); 

% *add a y-offset to the LFP to plot one above the other
lfp_cent = lfp_minmax * (2*ii-1);
data = data + lfp_cent;
% *plot the current piece of LFP
plot(tvec0, data, 'color', color);
 
end

% add a vertical line to the plot indicating time zero 
plot([0,0],[0, 2*ii*lfp_minmax],'-.r');
title('event LPF');
xlabel('time(s)');
ylabel('LFP');

end
% further niceties: add an option to decimate, to color, to filter before plotting