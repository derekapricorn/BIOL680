function shift_limits(src,event)
xlimits = get(gca, 'XLim');

xRange = xlimits(2) - xlimits(1);

% Callback to parse keypress event data to print a figure
if strcmp(event.Key,'rightarrow')
    disp('ra');
    xlimits = xlimits + xRange / 2;
    set(gca, 'XLim', xlimits);
end

if strcmp(event.Key,'leftarrow')
    disp('la');
    xlimits = xlimits - xRange / 2;
    set(gca, 'XLim', xlimits);
end

end
