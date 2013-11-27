%% sample code from the course website%%

cd 'D:\data\promoted\R042\R042-2013-08-18';
 
sd.fc = FindFiles('*.t');
sd.fc = cat(1,sd.fc,FindFiles('*._t')); % also load poorly isolated cells (if you don't have these, get from database)
sd.S = LoadSpikes(sd.fc);

% remember to unzip VT1.zip if not yet done
[Timestamps, X, Y, Angles, Targets, Points, Header] = Nlx2MatVT('VT1.nvt', [1 1 1 1 1 1], 1, 1, []);
Timestamps = Timestamps*10^-6;%convert to seconds

% remove zero samples, which happen when no LED is detected
toRemove = (X == 0 & Y == 0);
X = X(~toRemove); Y = Y(~toRemove); Timestamps = Timestamps(~toRemove);

sd.x = tsd(Timestamps,X');
sd.y = tsd(Timestamps,Y');
% find out the sampling frequency from the timestamps
S_freq = mode(diff(Timestamps));

%before looking at the data we first ignore the pre and post run segments
t = [3250 5650];
for iC = 1:length(sd.S)
    sd.S{iC} = Restrict(sd.S{iC},t(1),t(2));
end
sd.x = Restrict(sd.x,t(1),t(2));
sd.y = Restrict(sd.y,t(1),t(2));

%plot the position data
plot(Data(sd.x),Data(sd.y),'.','Color',[0.5 0.5 0.5],'MarkerSize',1);
axis off; hold on;



%plot the spikes of a single cell
% get x and y coordinate for times of spike
iC = 5;
spk_x = interp1(Range(sd.x),Data(sd.x),Data(sd.S{iC}),'linear');
spk_y = interp1(Range(sd.y),Data(sd.y),Data(sd.S{iC}),'linear');
 
h = plot(spk_x,spk_y,'.r');


%tuning curves
SET_xmin = 10; SET_ymin = 10; SET_xmax = 640; SET_ymax = 480;
SET_nxBins = 63; SET_nyBins = 47;
 
spk_binned = ndhist(cat(1,spk_x',spk_y'),[SET_nxBins; SET_nyBins],[SET_xmin; SET_ymin],[SET_xmax; SET_ymax]);
 
imagesc(spk_binned')
axis xy; colorbar


%spike count depending on the bin sizes
occ_binned = ndhist(cat(1,Data(sd.x)',Data(sd.y)'),[SET_nxBins; SET_nyBins],[SET_xmin; SET_ymin],[SET_xmax; SET_ymax]);
 
% this is a sample count, so need to convert to seconds (1/30s per sample) to get time
VT_Fs = 30;
tc = spk_binned./(occ_binned .* (1./VT_Fs)); % firing rate is spike count divided by time
 
pcolor(tc'); shading flat
axis xy; colorbar; axis off


%improving image quality by bining with finer scale
SET_nxBins = 630; SET_nyBins = 470; % more bins
 
kernel = gausskernel([30 30],8); % 2-D gaussian for smoothing: 30 points in each direction, SD of 8 bins
 
% spikes
spk_binned = ndhist(cat(1,spk_x',spk_y'),[SET_nxBins; SET_nyBins],[SET_xmin; SET_ymin],[SET_xmax; SET_ymax]);
spk_binned = conv2(spk_binned,kernel,'same');
 
% occupancy
occ_binned = ndhist(cat(1,Data(sd.x)',Data(sd.y)'),[SET_nxBins; SET_nyBins],[SET_xmin; SET_ymin],[SET_xmax; SET_ymax]);
occ_binned = conv2(occ_binned,kernel,'same');
 
occ_binned(occ_binned < 0.01) = 0;%improves the contrast
tc = spk_binned./(occ_binned .* (1 / VT_Fs));
tc(isinf(tc)) = NaN;%without this statement, the occ_binned == 1 case will give you inf for tc, creating useless regions
 
figure;
pcolor(tc'); shading flat; axis off
axis xy; colorbar


%preparing tuning curves for decoding
kernel = gausskernel([4 4],2); % 2-D gaussian, width 4 bins, SD 2
 
SET_xmin = 10; SET_ymin = 10; SET_xmax = 640; SET_ymax = 480;
SET_nxBins = 63; SET_nyBins = 47;
 
spk_binned = ndhist(cat(1,spk_x',spk_y'),[SET_nxBins; SET_nyBins],[SET_xmin; SET_ymin],[SET_xmax; SET_ymax]);
spk_binned = conv2(spk_binned,kernel,'same'); % smoothing
 
occ_binned = ndhist(cat(1,Data(sd.x)',Data(sd.y)'),[SET_nxBins; SET_nyBins],[SET_xmin; SET_ymin],[SET_xmax; SET_ymax]);
occ_mask = (occ_binned < 5);
occ_binned = conv2(occ_binned,kernel,'same'); % smoothing
 
occ_binned(occ_mask) = 0; % don't include bins with less than 5 samples
 
VT_Fs = 30;
tc = spk_binned./(occ_binned .* (1 / VT_Fs));
tc(isinf(tc)) = NaN;
 
pcolor(tc'); shading flat
axis xy; colorbar; axis off

% do the same for all the cells 
clear tc
nCells = length(sd.S);
for iC = 1:nCells
    spk_x = interp1(Range(sd.x),Data(sd.x),Data(sd.S{iC}),'linear');
    spk_y = interp1(Range(sd.y),Data(sd.y),Data(sd.S{iC}),'linear');
 
    spk_binned = ndhist(cat(1,spk_x',spk_y'),[SET_nxBins; SET_nyBins],[SET_xmin; SET_ymin],[SET_xmax; SET_ymax]);
    spk_binned = conv2(spk_binned,kernel,'same');
 
    tc = spk_binned./(occ_binned .* (1 / VT_Fs));
    tc(isinf(tc)) = NaN;
 
    sd.tc{iC} = tc;
 
end


%inspect the resulting tuning curves
ppf = 25; % plots per figure
for iC = 1:length(sd.S)
    nFigure = ceil(iC/ppf);
    figure(nFigure);
 
    subplot(5,5,iC-(nFigure-1)*ppf);
    pcolor(sd.tc{iC}); shading flat; axis off;
 
end

%if one cell has completel green place map, it probably means that cell is
%not associated with the cue


%preparing firing rates for decoding
clear Q;
binsize = 0.1;
tvec = t(1):binsize:t(2);
for iC = length(sd.S):-1:1
 
    spk_t = Data(sd.S{iC});
    Q(iC,:) = histc(spk_t,tvec);
 
end
nActiveNeurons = sum(Q > 0);


% look at it
imagesc(tvec,1:nCells,Q)
set(gca,'FontSize',16); xlabel('time(s)'); ylabel('cell #');


%final step
clear tc
nBins = numel(occ_binned);
nCells = length(sd.S);
for iC = nCells:-1:1
    tc(:,:,iC) = sd.tc{iC};
end
tc = reshape(tc,[size(tc,1)*size(tc,2) size(tc,3)]);
occUniform = repmat(1/nBins,[nBins 1]);


%decoding
len = length(tvec);
p = nan(length(tvec),nBins);
for iB = 1:nBins
    tempProd = nansum(log(repmat(tc(iB,:)',1,len).^Q));
    tempSum = exp(-binsize*nansum(tc(iB,:),2));
    p(:,iB) = exp(tempProd)*tempSum*occUniform(iB);
end

p = p./repmat(sum(p,2),1,nBins);
p(nActiveNeurons < 1,:) = 0;


%convert the actual position into binned form
xBinEdges = linspace(SET_xmin,SET_xmax,SET_nxBins+1);
yBinEdges = linspace(SET_ymin,SET_ymax,SET_nyBins+1);
 
xTempD = Data(sd.x); xTempR = Range(sd.x);
yTempD = Data(sd.y);
 
gS = find(~isnan(xTempD) & ~isnan(yTempD));
xi = interp1(xTempR(gS),xTempD(gS),tvec,'linear');
yi = interp1(xTempR(gS),yTempD(gS),tvec,'linear');
 
xBinned = (xi-xBinEdges(1))./median(diff(xBinEdges));
yBinned = (yi-yBinEdges(1))./median(diff(yBinEdges));

goodOccInd = find(occ_binned > 0);
for iT = 1:size(p,1)
    cla;
    temp = reshape(p(iT,:),[SET_nxBins SET_nyBins]);
    toPlot = nan(SET_nxBins,SET_nyBins);
    toPlot(goodOccInd) = temp(goodOccInd);
 
    pcolor(toPlot); axis xy; hold on; caxis([0 0.5]);
    shading flat; axis off;
 
    hold on; plot(yBinned(iT),xBinned(iT),'ow','MarkerSize',15);
 
    h = title(sprintf('t %.2f, nCells %d',tvec(iT),nActiveNeurons(iT))); 
    if nActiveNeurons(iT) == 0
        set(h,'Color',[1 0 0]);
    else
        set(h,'Color',[0 0 0]);
    end
    drawnow; pause(0.1);
end


%export to movie
h = figure; set(h,'Position',[100 100 320 240]);


%%
%assingment%
%find the error distance 
temp = (nActiveNeurons~=0);%find inactive bins
xBinned(~temp)=0;
xi(~temp)=0;
yBinned(~temp)=0;
yi(~temp)=0;

dist = sqrt((xBinned-xi).^2 + (yBinned - yi).^2);
plot(tvec,dist,'.');
title('decoding error');xlabel('time(s)');ylabel('distance(pixels)');
%by changing the binsize from 0.25 to 0.1, we see a more condensed plot
%with more scattering points, meaning more decoding errors to occur. That's
%probably because as the binsize decreases, the estimation become less
%accurate, resulting in more computational errors
%
