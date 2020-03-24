function [S_sum, f,t] = spectr(specs, neededRangeBins,fs_slow ) 
% THE FUCTION computes the spectrogram from the fft of the range map. There
% are two frequency clusteres available -- either at 2GHz and or at 4GHz.
% The function uses the 2GHz frequency bins. 
% input: 
% fft matrix            --> specs
% selected bins         --> neededRangeBins
% sample frequency      --> fs_slow
% returns: 
% computed spectrogram  --> S_sum
% Frequency vector      --> f
% time vector           --> t

for ki = 1:numel(neededRangeBins)
    specs111 = specs(neededRangeBins(ki),:);
    [S,f,t,P] = spectrogram(detrend(specs111,2),hann(128),120,128,fs_slow);
    % suming up the spectrograms
    if ki == 1; S_sum = zeros(size(S)); end 
    S_sum = S_sum + abs(S); 
    fig = figure(1);
    sspec = fftshift(10*log10(P/max(P(:)))',2)';
    imagesc(t, f-f(end)/2, sspec);
    colormap jet; colorbar; caxis([max(max(sspec))-50, inf]); axis xy
    title(['Spectrogram: ', num2str(ki)])
%     pause(0.1);
end