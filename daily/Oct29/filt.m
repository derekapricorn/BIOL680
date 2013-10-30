function [y] = filt(x, type, passband, Fs)
%constructs a filter
%type: 1 - butt; 2 - cheby
%stop band = passband +- 5
%shape: 1 - low-pass; 2 - band-pass
    

Wp = [passband(1),passband(2)] * 2 / Fs; % passband - between 50 and 100 Hz
Ws = [passband(1)-5,passband(2)+5] * 2 / Fs; % stopband
if(type == 1)
    [N,Wn] = buttord( Wp, Ws, 3, 20); % determine filter parameters
    [b,a] = butter(N,Wn); % builds filter
else
    [N,Wn] = cheb1ord( Wp, Ws, 3, 20); 
    [b,a] = cheby1(N,0.5,Wn);
end
    y = filter(b,a,x);
    
end