function [psd_est,psd_est_f,psd_orig,psd_orig_f] = psdcom (attz_est,attz_orig,h)

% Pxx_true = periodogram(wqzI(2,:));
fs=1/h; %HzÎªµ¥Î»
n = length(attz_est(1,:));
yushu = mod(n,2);

% if yushu==0
    nfft = n;
% else
%     nfft = n-1;
% end
% nfft=5380;

window=boxcar(n);
[Pxx_est1,f1] = periodogram(attz_est(1,:),window,nfft,fs);   
[Pxx_est2,f2] = periodogram(attz_est(2,:),window,nfft,fs);   
[Pxx_est3,f3]= periodogram(attz_est(3,:),window,nfft,fs);   
psd_est = [Pxx_est1 Pxx_est2 Pxx_est3];
psd_est_f = f1 ;

[Pxx_orig1,of1] = periodogram(attz_orig(1,:),window,nfft,fs);   
[Pxx_orig2,of2] = periodogram(attz_orig(2,:),window,nfft,fs);   
[Pxx_orig3,of3] = periodogram(attz_orig(3,:),window,nfft,fs);
psd_orig = [Pxx_orig1 Pxx_orig2 Pxx_orig3];
psd_orig_f = of1;
