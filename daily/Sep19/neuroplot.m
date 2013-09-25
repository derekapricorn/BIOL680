function neuroplot(spikes,csc_sample,varargin)

% inputs:
%
% spikes: {nCells x 1} cell array of ts objects (spike train)
% csc: {nCSCs x 1} cell array of tsd objects (LFPs)
%
% varargins:
%
% cscColor: [nCSC x 3] RGB values to plot CSCs, default []
% spikeColor: [nCells x 3] RGB values to plot spikes, default []
%
% evt: [nEvents x 1] event times to plot, default []
% evtColor: [nEvents x 3] RGB values to plot events, default []
%
% interactiveMode: boolean to enable/disable arrow key navigation

cscColor = [];
spikeColor = [];
evt = [];
evtColor = [];
interactiveMode = 1;
extract_varargin; %override variables if needed



figure
hold on
% the shift_limits.m is located in the same repository
if(interactiveMode)
    f_hdl = figure('KeyPressFcn',@shift_limits);
end

% raster plot
for i = 1:length(spikes)
    temp = Data(spikes{i});
    for j = 1:length(temp);
        if isempty(spikeColor)
            line([temp(j),temp(j)],[i-1,i], 'Color', 'g');
        else
            line([temp(j),temp(j)],[i-1,i],'Color', spikeColor(i));
        end
    end
end

% LFP plot
DC_offset = length(spikes)+5;
for i = 1:length(csc_sample)
     temp = csc_sample(i);
     
     if isempty(cscColor)
         plot(Range(temp), Data(temp)*0.001+ DC_offset, 'Color', 'm');
     else
         plot(Range(temp), Data(temp)*0.001+ DC_offset, 'Color', cscColor(i));
     end
     
     DC_offset = DC_offset + 1;
end



%event plot
if ~isempty(evt)
      
    Y_limits = get(gca, 'YLim');
    if isempty(evtColor)
        for i = 1:length(evt)
            line([evt(i),evt(i)],[0, Y_limits(2)], 'Color', 'k');
        end
    else
        for i = 1:length(evt)
            line([evt(i),evt(i)],[0, Y_limits(2)], 'Color', evtColor(i));
        end
    end 
end









