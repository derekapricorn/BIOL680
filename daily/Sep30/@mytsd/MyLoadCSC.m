function csc_new = MyLoadCSC(fname)

% function csc = LoadCSC(fname)
%
% load a Neuralynx .ncs file correctly

% INPUTS:
% fname: [1 x n] char, i.e. a filename such as 'R042-2013-08-18-CSC01a.ncs'
%cd('D:\data\promoted\R016\R016-2012-10-08');
%fname = 'R016-2012-10-08-CSC03b.Ncs';
% OUTPUTS:
%
% csc: [1 x 1] mytsd

% cd to your location here
[Timestamps, ~, SampleFrequencies, NumberOfValidSamples, Samples, Header] = Nlx2MatCSC(fname, [1 1 1 1 1], 1, 1, []);

%convert unit of timestamps to seconds
Timestamps = Timestamps / 10^6;

%
%convert samples to milivolts
Samples = Samples * str2double(Header{15}(12:end));
%

% return if there are negative diffs
if any(diff(Timestamps) < 0)
    fprintf('negative diff, invalid input');
end  

%find invalid entries
is_valid = zeros(size(Samples,1),size(Samples,2));
for i = 1:size(Samples,2)
    if(NumberOfValidSamples(i) == 512)
        is_valid(:,i) = 1;
    else
        for j=1:NumberOfValidSamples(i)
            is_valid(j:i) = 1;
        end
    end
end
clear i;
clear j;
is_valid = logical(reshape(is_valid, [size(is_valid, 1)*size(is_valid, 2) 1]));

% unwrap samples, and remove the invalid points
csc_data = reshape(Samples,[size(Samples,1)*size(Samples,2) 1]);
csc_data = csc_data(is_valid);   

% construct clean timestamps
csc_timestamps = repmat(Timestamps,[size(Samples,1) 1]);
dtvec = (0:size(Samples,1)-1)*(1/2000);
dtmat = repmat(dtvec',[1 size(Samples,2)]);

csc_timestamps = csc_timestamps+dtmat;
csc_timestamps = reshape(csc_timestamps,[size(csc_timestamps,1)*size(csc_timestamps,2) 1]);

csc_timestamps = csc_timestamps(is_valid);

%done
csc_new = mytsd(csc_timestamps, csc_data, Header);

% old method
[csc_old,csc_info_old] = LoadCSC('R016-2012-10-08-CSC03b.Ncs');
tvec = Range(csc_old);
raw_LFP = Data(csc_old);

%plot. The middle plot shows the outliers
outliers = find(is_valid == 0);
figure

ax(1) = subplot(311);
plot(tvec,raw_LFP,'b')

ax(2) = subplot(312);
plot(tvec(outliers),raw_LFP(outliers),'k')

ax(3) = subplot(313);
plot(Range(csc_new),Data(csc_new),'r')

set(ax, 'xlim', [500, 3600]);

end
