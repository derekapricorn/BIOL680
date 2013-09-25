%% load the data (note, may need to unzip position data first)
cd D:\data\promoted\R042\R042-2013-08-18;
fc = FindFiles('*.t');
S = LoadSpikes(fc);
 
[csc,csc_info] = LoadCSC('R042-2013-08-18-CSC03a.ncs');
 
[Timestamps, X, Y, Angles, Targets, Points, Header] = Nlx2MatVT('VT1.nvt', [1 1 1 1 1 1], 1, 1, [] );

%% verify csc object indeed has timestamps and data of the same length
length(Range(csc)) - length(Data(csc))
    

%% use Restrict() to extract data between 5950 and 6050s

csc_short = Restrict(csc, 5950, 6050);

%% normalize tsd for X and Y
X_sec = tsd(Timestamps*10^-6, X', 'sec');
Y_sec = tsd(Timestamps*10^-6, Y', 'sec');
X_short = Restrict(X_sec, 5950, 6050);

%% verify that the new truncated tsd's have the same length
EndTime(csc_short) - EndTime(X_short) 

%% plot the LFP for the segment between 5950 and 6050
plot(Range(csc_short), Data(csc_short));
hold on
%%  change the figure background to black
set(gcf,'Color',[0 0 0]);
%% properties that correspond to the axes color
get(gca)
%% set axis color to white
set(gca, 'Color', [1 1 1]);

%% save the image
print(gcf,'-dpng','-r300','R042-2013-08-18-LFPsnippet.png');

%% set invertHardcopy off
set(gcf, 'InvertHardcopy', 'off');

%%
hold on; box off;

csc_mean = nanmean(Data(csc));
xr = get(gca,'XLim');
%plot a red, dashed line with width 2
mean_hdl = plot(xr,[csc_mean csc_mean],'r:','LineWidth',2);

%%
set(gca, 'XLim', [5989, 5990], 'FontSize',[24]);

%% anonymous function
sqr_fn = @(x) x.^2;
sqr_fn(2)

%% the shift_limits.m is located in the same repository
f_hdl = figure('KeyPressFcn',@shift_limits);

%% test fun
test_fun(1, 'b', 0, 'c', 0)
a is 1, b is 0, c is 0



