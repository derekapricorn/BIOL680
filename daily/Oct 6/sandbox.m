%%
%Assignment
%*************************************************************************
%part 1
clear all
close all
%Compute the PSD of “white noise”, i.e. a random signal where each sample is
%drawn independently from the open interval (0,1) with equal probability. 
%Is it 1/f? How would you describe its shape? Hint: use the MATLAB function rand(). 
endpoint = 1000;
x = rand(1000,1);
n = 0:endpoint-1;
Fs = 2000;

wSize = 512;

[Pxx_w,F] = pwelch(x, hamming(wSize),wSize/2,length(x),Fs);

[Pxx_p,F] = periodogram(x,hamming(length(x)),length(x),Fs);
plot(F,10*log10(Pxx_w),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
hold on
plot(F,10*log10(Pxx_p),'b'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 150]);
%**************************************************************************
%%
%Part 2
close all
cd('D:\data\promoted\R016\R016-2012-10-08');
%csc_vs = MyLoadCSC('R016-2012-10-08-CSC04d.ncs');
run(FindFile('*keys.m'));
hippo_LFP = getfield(ExpKeys, 'goodSWR');
hippo_LFP = strrep(hippo_LFP, '3', '8');
%csc_hippo = MyLoadCSC(hippo_LFP);  //why does it not work!!!!!?????
csc_hippo = MyLoadCSC('R016-2012-10-08-CSC02b.ncs');
csc_hippoR = Range(csc_hippo);
csc_hippoD = Data(csc_hippo); 
csc_vs = MyLoadCSC('R016-2012-10-08-CSC04d.ncs');
csc_vsR = Range(csc_vs);
csc_vsD = Data(csc_vs);

wSize = 1024;
[Pxx_hippo,F_hippo] = periodogram(csc_hippoD,hamming(length(csc_hippoD)),length(csc_hippoD),Fs);
plot(F_hippo,10*log10(Pxx_hippo),'r'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 150]);
hold on
[Pxx_vs,F_vs] = periodogram(csc_vsD,hamming(length(csc_vsD)),length(csc_vsD),Fs);
plot(F_vs,10*log10(Pxx_vs),'b'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 150]);
%the hippocampal LFP has bigger variation throughout the spectrum. At the 
%low freq end, the hippocampal LFP expresses higher intensity.
%**************************************************************************
%%
%Part 3
%hippocampal data
Fs = 2000;
nP = 1024;
csc_pre = Restrict(csc_hippo,0,ExpKeys.TimeOnTrack(1)-10);
csc_preR = Range(csc_pre);
csc_preD = Data(csc_pre);

figure
    iw = 40;%set window size
    [Pxx,F] = pwelch(csc_preD,hamming(iw),iw/2,nP,Fs);
    plot(F,Pxx,'r'); xlabel('Frequency (Hz)');
    hold on
    iw = 200;
    [Pxx,F] = pwelch(csc_preD,hamming(iw),iw/2,nP,Fs);
    plot(F,Pxx,'g'); xlabel('Frequency (Hz)');
    iw = 400;
    [Pxx,F] = pwelch(csc_preD,hamming(iw),iw/2,nP,Fs);
    plot(F,Pxx,'b'); 
    xlim([0, 200]);
    set(gcf,'Name','hippocampal graph')
    hleg1 = legend('iw = 40','iw = 200','iw = 400');
%ventral striatal data
    csc_pre = Restrict(csc_vs,0,ExpKeys.TimeOnTrack(1)-10);
    csc_preR = Range(csc_pre);
    csc_preD = Data(csc_pre);
    figure
    iw = 40;%set window size
    [Pxx,F] = pwelch(csc_preD,hamming(iw),iw/2,nP,Fs);
    plot(F,Pxx,'r'); xlabel('Frequency (Hz)');
    hold on
    iw = 200;
    [Pxx,F] = pwelch(csc_preD,hamming(iw),iw/2,nP,Fs);
    plot(F,Pxx,'g'); xlabel('Frequency (Hz)');
    iw = 400;
    [Pxx,F] = pwelch(csc_preD,hamming(iw),iw/2,nP,Fs);
    plot(F,Pxx,'b'); 
    xlim([0, 200]);
    set(gcf,'Name','ventral striatal graph')
    hleg1 = legend('iw = 40','iw = 200','iw = 400');
    
 %the larger the window size, the closer the plot is to the orignal
 %spectrum. The small window gives a smooth and steady curve
 %**************************************************************************
%%
%Part 4
dsf = 10;
%method 1
csc_hippoD_1 = decimate(csc_hippoD,dsf);
csc_hippoR_1 = decimate(csc_hippoR,dsf);

%method 2
csc_hippoR_2 = csc_hippoR(1:dsf:end);
csc_hippoD_2 = csc_hippoD(1:dsf:end);

figure
[Pxx,F] = periodogram(csc_hippoD,hamming(length(csc_hippoD)),length(csc_hippoD),Fs);
ax(1) = subplot(311);
plot(F,10*log10(Pxx),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 100]);
[Pxx,F] = periodogram(csc_hippoD_1,hamming(length(csc_hippoD_1)),length(csc_hippoD_1),Fs/dsf);
ax(2) = subplot(312);
plot(F,10*log10(Pxx),'g'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 100]);
[Pxx,F] = periodogram(csc_hippoD_2,hamming(length(csc_hippoD_2)),length(csc_hippoD_2),Fs/dsf);
ax(3) = subplot(313);
plot(F,10*log10(Pxx),'b'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 100]);

%at higher frequency, the decimate method(1) loses fidelity, whereas the
%fundamental downsampling method (2) still works okay




%%
%practice code
%**************************************************************************
% x = round(rand(1,8)*10); % generate a length 8 vector of integers between 0 and 10
% xlen = length(x);
%  
% % get magnitudes and phases of Fourier series
% X = fft(x);
% Xmag = abs(X); % magnitudes, a_n
% Xphase = angle(X); % phases, phi_n
%  
% n = 0:xlen-1;
% t = 0:0.05:xlen-1; % a finer timescale to show the smooth signal later
%  
% for iH = xlen-1:-1:0 % reconstruct each harmonic
%     s(iH+1,:) = Xmag(iH+1)*cos(2*pi*iH*n/xlen + Xphase(iH+1))/xlen;
%     sm(iH+1,:) = Xmag(iH+1)*cos(2*pi*iH*t/xlen + Xphase(iH+1))/xlen;
%     % detail: xlen appears here because the fundamental frequency used by fft() depends on this
% end
%  
% ssum = sum(s);
% smsum = sum(sm);
%  
% figure;
% plot(n, x, 'go', t, smsum, 'b', n, ssum, 'r*');
% legend({'original','sum - all','sum - points only'});
% 
% %%
% Fs = 100; % in samples per second (Hz)
% t0 = 0; t1 = 1; % start and end times
% tvec = t0:1/Fs:t1-(1/Fs); % construct time axis; generate exactly 20 samples
%  
% f = 2; % signal frequency
% y = sin(2*pi*f*tvec); % construct signal, a 2Hz sine wave sampled at 20Hz for 1s
%  
% yfft = fft(y,length(y));
% yfft_mag = abs(yfft); yfft_ph = angle(yfft);
% stem(yfft_mag)
% 
% 
% Npoints = length(y);
% F = [-Npoints/2:Npoints/2-1]./Npoints; % construct frequency axis
%  
% yfft_mag = fftshift(yfft_mag); % align output, see note below
% stem(F,yfft_mag);
%  
% xlabel('Frequency (Fs^{-1})');
% 
% %%
% Fs = 20;
% tvec = t0:1/Fs:t1;
% nPoints = [length(tvec) 64 256 1024];
%  
% for iP = 1:length(nPoints) % repeat fft with different numbers of points
%  
%     nP = nPoints(iP);
%     subplot(2,2,iP);
%  
%     y = sin(2*pi*f*tvec);
%     yfft = fft(y,nP);
%     yfft_mag = abs(yfft); yfft_ph = angle(yfft);
%  
%     F = [-nP/2:nP/2-1]./nP;
%     yfft_mag = fftshift(yfft_mag);
%     plot(F,yfft_mag,'kx',F,yfft_mag,'k');
%  
%     title(sprintf('%d point FFT',nP));
%     xlabel('Frequency (Fs^{-1})');
%  
% end
% 
% %%
% nP = 25;
% nPFFT = 1024;
%  
% windows = {'rectwin','triang','hamming','hanning','blackman'};
% cols = 'rgbcmyk';
%  
% for iW = 1:length(windows)
%  
%     eval(cat(2,'wn = ',windows{iW},'(nP);')); % make sure you understand this
%     wn = wn./sum(wn);
%  
%     subplot(211); % plot the window
%     plot(wn,cols(iW),'LineWidth',2); hold on;
%  
%     subplot(212);
%     yfft = fft(wn,nPFFT);
%     yfft_mag = abs(yfft); yfft_ph = angle(yfft);
%  
%     F = [-nPFFT/2:nPFFT/2-1]./nPFFT;
%     yfft_mag = fftshift(yfft_mag);
%  
%     h(iW) = plot(F,yfft_mag,cols(iW),'LineWidth',2); hold on;
%  
% end
%  
% xlabel('Frequency (Fs^{-1})');
% legend(h,windows);
% set(gca, 'XScale', 'log')
% %%
% tvec = t0:1/Fs:t1-(1/Fs);
% nRepeats = [1 2 4 8];
%  
% nP =  1024;
%  
% for iP = 1:length(nRepeats)
%  
%     subplot(2,2,iP);
%  
%     y = sin(2*pi*f*tvec);
%     y = repmat(y,[1 nRepeats(iP)]); % repeat the signal a number of times
%  
%     yfft = fft(y,nP);
%     yfft_mag = abs(yfft); yfft_ph = angle(yfft);
%  
%     F = [-nP/2:nP/2-1]./nP;
%     yfft_mag = fftshift(yfft_mag);
%     plot(F,yfft_mag,'kx',F,yfft_mag,'k');
%  
%     title(sprintf('%d repeats',nRepeats(iP)));
%     xlabel('Frequency (Fs^{-1})');
%  
% end
% %% successfully loaded
% tvec = t0:1/Fs:t1-(1/Fs);
% nRepeats = [1 2 4 8];
%  
% nP =  1024;
% figure
% for iP = 1:length(nRepeats)
%  
%     subplot(2,2,iP);
%  
%     y_h = sin(2*pi*f*tvec) .* rectwin(length(tvec))';
%     y_h = repmat(y_h,[1 nRepeats(iP)]); % repeat the signal a number of times
%  
%     yfft_h = fft(y_h,nP);
%     yfft_mag = abs(yfft_h); yfft_ph = angle(yfft_h);
%  
%     F = [-nP/2:nP/2-1]./nP;
%     yfft_mag = fftshift(yfft_mag);
%     plot(F,yfft_mag,'kx',F,yfft_mag,'k');
%  
%     title(sprintf('%d repeats',nRepeats(iP)));
%     xlabel('Frequency (Fs^{-1})');
%  
% end
%     
%     
% %%
% Fs = 20; % in samples per second (Hz)
% t0 = 0; t1 = 1;
% f = 2;
% nRepeats = 4;
%  
% tvec = t0:1/Fs:t1-(1/Fs);
%  
% nP =  1024;
% y = sin(2*pi*f*tvec);
% y = repmat(y,[1 nRepeats]);
%  
% [Pxx,F] = periodogram(y,rectwin(length(y)),nP,Fs);
% plot(F,Pxx);
%  
% hold on;
% wSize = 40;
% [Pxx,F] = pwelch(y,rectwin(wSize),wSize/2,nP,Fs);
% plot(F,Pxx,'r'); xlabel('Frequency (Hz)');
% 
% %%
% Fs = 20; % in samples per second (Hz)
% t0 = 0; t1 = 1; f = 2;
% nP =  1024;
% gaps = [5 10 15]; % idx of samples to be removed
%  
% tvec = t0:1/Fs:t1;%-(1/Fs);
% y = sin(2*pi*f*tvec);
%  
% subplot(211)
% plot(tvec,y,'k*'); hold on;
%  
% yfft = fft(y,nP);
% yfft_mag = abs(yfft); yfft_ph = angle(yfft);
%  
% F = [-nP/2:nP/2-1]./nP;
% yfft_mag = fftshift(yfft_mag);
%  
% subplot(212);
% plot(F,yfft_mag,'kx',F,yfft_mag,'k'); hold on;
%  
% xlabel('Frequency (Fs^{-1})');
%  
% % signal with gaps
% y = sin(2*pi*f*tvec);
% y2(gaps) = []; tvec(gaps) = []; % remove
%  
% subplot(211);
% plot(tvec,y2,'bo'); hold on;
%  
% yfft = fft(y2,nP);
% yfft_mag = abs(yfft); yfft_ph = angle(yfft);
%  
% F = [-nP/2:nP/2-1]./nP;
% yfft_mag = fftshift(yfft_mag);
%  
% subplot(212);
% plot(F,yfft_mag,'bx',F,yfft_mag,'b');
% 
% 
% 
% %%
% % cd to your R016-2012-10-08 folder
% cd('D:\data\promoted\R016\R016-2012-10-08');
% csc = MyLoadCSC('R016-2012-10-08-CSC04d.ncs');
% run(FindFile('*keys.m'));
% 
% % restrict to prerecord, leaving some time (10s) before rat actually goes on track
% csc_pre = Restrict(csc,0,ExpKeys.TimeOnTrack(1)-10);
% csc_preR = Range(csc_pre);
% csc_preD = Data(csc_pre);
%  
% % check if sampling is ok
% plot(diff(csc_preR)); % only minimal differences
% Fs = 1./mean(diff(csc_preR));



%%
%Assignment
%*************************************************************************
%part 1
clear all
close all
%Compute the PSD of “white noise”, i.e. a random signal where each sample is
%drawn independently from the open interval (0,1) with equal probability. 
%Is it 1/f? How would you describe its shape? Hint: use the MATLAB function rand(). 
endpoint = 1000;
x = rand(1000,1);
n = 0:endpoint-1;
Fs = 2000;

wSize = 512;

[Pxx_w,F] = pwelch(x, hamming(wSize),wSize/2,length(x),Fs);

[Pxx_p,F] = periodogram(x,hamming(length(x)),length(x),Fs);
plot(F,10*log10(Pxx_w),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
hold on
plot(F,10*log10(Pxx_p),'b'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 150]);
%**************************************************************************
%%
%Part 2
close all
cd('D:\data\promoted\R016\R016-2012-10-08');
%csc_vs = MyLoadCSC('R016-2012-10-08-CSC04d.ncs');
run(FindFile('*keys.m'));
hippo_LFP = getfield(ExpKeys, 'goodSWR');
hippo_LFP = strrep(hippo_LFP, '3', '8');
%csc_hippo = MyLoadCSC(hippo_LFP);  //why does it not work!!!!!?????
csc_hippo = MyLoadCSC('R016-2012-10-08-CSC02b.ncs');
csc_hippoR = Range(csc_hippo);
csc_hippoD = Data(csc_hippo); 
csc_vs = MyLoadCSC('R016-2012-10-08-CSC04d.ncs');
csc_vsR = Range(csc_vs);
csc_vsD = Data(csc_vs);

wSize = 1024;
[Pxx_hippo,F_hippo] = periodogram(csc_hippoD,hamming(length(csc_hippoD)),length(csc_hippoD),Fs);
plot(F_hippo,10*log10(Pxx_hippo),'r'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 150]);
hold on
[Pxx_vs,F_vs] = periodogram(csc_vsD,hamming(length(csc_vsD)),length(csc_vsD),Fs);
plot(F_vs,10*log10(Pxx_vs),'b'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 150]);
%the hippocampal LFP has bigger variation throughout the spectrum. At the 
%low freq end, the hippocampal LFP expresses higher intensity.
%**************************************************************************
%%
%Part 3
%hippocampal data
Fs = 2000;
nP = 1024;
csc_pre = Restrict(csc_hippo,0,ExpKeys.TimeOnTrack(1)-10);
csc_preR = Range(csc_pre);
csc_preD = Data(csc_pre);

figure
    iw = 40;%set window size
    [Pxx,F] = pwelch(csc_preD,hamming(iw),iw/2,nP,Fs);
    plot(F,Pxx,'r'); xlabel('Frequency (Hz)');
    hold on
    iw = 200;
    [Pxx,F] = pwelch(csc_preD,hamming(iw),iw/2,nP,Fs);
    plot(F,Pxx,'g'); xlabel('Frequency (Hz)');
    iw = 400;
    [Pxx,F] = pwelch(csc_preD,hamming(iw),iw/2,nP,Fs);
    plot(F,Pxx,'b'); 
    xlim([0, 200]);
    set(gcf,'Name','hippocampal graph')
    hleg1 = legend('iw = 40','iw = 200','iw = 400');
%ventral striatal data
    csc_pre = Restrict(csc_vs,0,ExpKeys.TimeOnTrack(1)-10);
    csc_preR = Range(csc_pre);
    csc_preD = Data(csc_pre);
    figure
    iw = 40;%set window size
    [Pxx,F] = pwelch(csc_preD,hamming(iw),iw/2,nP,Fs);
    plot(F,Pxx,'r'); xlabel('Frequency (Hz)');
    hold on
    iw = 200;
    [Pxx,F] = pwelch(csc_preD,hamming(iw),iw/2,nP,Fs);
    plot(F,Pxx,'g'); xlabel('Frequency (Hz)');
    iw = 400;
    [Pxx,F] = pwelch(csc_preD,hamming(iw),iw/2,nP,Fs);
    plot(F,Pxx,'b'); 
    xlim([0, 200]);
    set(gcf,'Name','ventral striatal graph')
    hleg1 = legend('iw = 40','iw = 200','iw = 400');
    
 %the larger the window size, the closer the plot is to the orignal
 %spectrum. The small window gives a smooth and steady curve
 %**************************************************************************
%%
%Part 4
dsf = 10;
%method 1
csc_hippoD_1 = decimate(csc_hippoD,dsf);
csc_hippoR_1 = decimate(csc_hippoR,dsf);

%method 2
csc_hippoR_2 = csc_hippoR(1:dsf:end);
csc_hippoD_2 = csc_hippoD(1:dsf:end);

figure
[Pxx,F] = periodogram(csc_hippoD,hamming(length(csc_hippoD)),length(csc_hippoD),Fs);
ax(1) = subplot(311);
plot(F,10*log10(Pxx),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 100]);
[Pxx,F] = periodogram(csc_hippoD_1,hamming(length(csc_hippoD_1)),length(csc_hippoD_1),Fs/dsf);
ax(2) = subplot(312);
plot(F,10*log10(Pxx),'g'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 100]);
[Pxx,F] = periodogram(csc_hippoD_2,hamming(length(csc_hippoD_2)),length(csc_hippoD_2),Fs/dsf);
ax(3) = subplot(313);
plot(F,10*log10(Pxx),'b'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 100]);

%at higher frequency, the decimate method(1) loses fidelity, whereas the
%fundamental downsampling method (2) still works okay
