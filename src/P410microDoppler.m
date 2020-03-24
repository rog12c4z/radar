%##########################################################################
%#####  Progrtram converts the EXCEL file into a *.mat file ###############
%##########################################################################
%
%--------------------------------------
% Author:       Ronny (Gerhard) Guendel
% Written by:   Microwave Sensing, Signals and Systems (MS3)
% University:   TU Delft
% Email:        r.guendel@tudelft.nl
% Created:      01/15/2020
% Updated:      03/24/2020

%% set clear
clc; clear all;   %close all;
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultFigureWindowStyle','docked');
set(0,'Defaultlinelinewidth',2);

%% for the EXCEL file
% [cfg,req,scn,det] = readMrmRetLog;

% for the matlab file *.mat
[file,path] = uigetfile('C:\Research\PulsON check\threeRadCol','*.mat');
load([path,file]);

%% parameter and vector scan
NScans = length(scn);
scans = [scn.scn];
NTS = scn.Nscn;

%% MTI Filter
h = [1 -2 1];
scansMTI= filter(h,1,scans,[],2);
% scansMTI= scans;
scansMTI = reshape(scansMTI,[],NScans);
dt_ps = (cfg.Tstp-cfg.Tstrt)/1e12;
fs = 1e12/dt_ps; %[Hz]

fs_slow = 1/7.5e-3;
T = 1/fs_slow*size(scansMTI,2);

c = 3e8; 
Rmin = 1; % 1meter
Rmax = (cfg.Tstp-cfg.Tstrt)/1e12*c/2 + Rmin;
%% Subtracting the mean form the in-phase and quatrature components
scansMTI = scansMTI-mean(scansMTI,2);
fig = figure;
imagesc([0 T],[Rmin Rmax],10*log10(abs(scansMTI))); axis xy; title(['Range map -- ']);%, file(1:3)])
colormap jet
colorbar;
caxis([max(max(10*log10(abs(scansMTI))))-30, inf]);

%% hilbert enveloping
% hil = hilbert(scans); % Raw data
hil = hilbert(scansMTI(:)); % after MTI and mean substraction

% plot(scans(NTS*20:NTS*22)); hold on;
% plot(real(hil(NTS*20:NTS*22))); hold on;
% plot(imag(hil(NTS*20:NTS*22))); hold on;
% plot(abs(hil(NTS*20:NTS*22))); hold on;
% legend('original', 'in-phase of hilbert', 'quadrature of hilbert', 'absoulute Value of hilbert');

hil_resha = reshape(hil,NTS,NScans);
RM_hil = 10*log10(abs(hil_resha));
fig = figure;
imagesc([0 T],[Rmin Rmax],RM_hil);title('Range map'); axis xy
colormap jet
colorbar;
caxis([max(max(RM_hil))-30, inf]);

%% spectrogram
specs = fft(scansMTI);
neededRangeBins = [70:110]; %[70:110]; % ,270:310 The best for 384 bins
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
figure;
sspec_sum = fftshift(10*log10(S_sum)',2)';
imagesc(t, f-f(end)/2,sspec_sum); axis xy
colormap jet; colorbar

%     savefig([name(ki,:),'.fig']); saveas(fig,[name(ki,:),'.png']);
%     saveas(fig,[name(ki,:),'.pdf']);

%% time stamp
figure;
plot(diff([scn.Tstmp]),'.');

return; 
%% kernel cleaning
% starting point is the range-map from the hilbert transformation: RM_hil

%% range map and spectorgram rezizing
psd_resize = imresize(RM_hil,[128 64*T]);  % RM image
sspec_resize = mat2gray(imresize(sspec_sum,[128 64*T]));  % spectrogram image


%% columnwise normilization
maxVec = max(psd_resize);           % determine the maximum value
maxVec = movmean(maxVec,10);        % Filter function
psd_norm = psd_resize./maxVec;      % columnwise normilization
plot_Image(psd_norm,T,'RM')

%% eClean alorithm
psd_eC = eCLEAN_Mod(psd_norm,0,8);
plot_Image(psd_eC,T,'RM')

sspec_eC = eCLEAN_Mod(sspec_resize,0,15);
plot_Image(sspec_eC,T,'MD')

%% outlier removal
OLR = bwareaopen(psd_eC,10);        % Function removes outliers logical operator
psd_OLR = zeros(size(psd_eC));      % creating a matrix of the same size
psd_OLR(OLR==1) = psd_eC(OLR==1);   % copying values from OLR into psd_OLR
plot_Image(psd_OLR,T,'RM')

OLR = bwareaopen(sspec_eC,100);        % Function removes outliers logical operator
sspec_OLR = zeros(size(sspec_eC));      % creating a matrix of the same size
sspec_OLR(OLR==1) = sspec_eC(OLR==1);   % copying values from OLR into psd_OLR
plot_Image(sspec_OLR,T,'MD')

%% write the range-map into a larger matrix; add paddings above and below
[m,n] = size(OLR);
padding = 16;
OLR_2 =[zeros(padding,n);psd_OLR;zeros(padding,n)];
plot_Image(OLR_2,T,'RM')

%% kernel cleaning method
[psd_kernel,psd_binary,psdRzVec] = RMkernalClean(OLR_2,1,12,12); % kernel cleaning algorithm

%% remove the paddings above and below
RM_clean = psd_kernel(padding+1:end-padding,:);
plot_Image(RM_clean,T,'RM');
RM_clean_bin = psd_binary(padding+1:end-padding,:);
plot_Image(RM_clean_bin,T,'RM');
% close all;
RM_clean_vec = psdRzVec(padding+1:end-padding,:);
plot_Image(RM_clean_vec,T,'RM');

%% discrete first and second derivative

[yPoints, xPoints] = find(psdRzVec ==1);
% smoothing the data
Methods = ["movmean" ; "movmedian" ; "gaussian" ; "lowess" ; "loess" ; "rlowess" ; "rloess" ; "sgolay"];
for ki = 3 % go for gaussian 1:length(Methods)
    yPointsSmo = smoothdata(yPoints,Methods(ki),100);
    figure(ki);
    subplot(511); plot(xPoints, yPointsSmo,'.b'); hold on; legend(strcat('f filtered--',Methods(ki)));
    subplot(512); plot( diff(yPointsSmo),'.r'); legend(strcat('df unfiltered--',Methods(ki)));
    % first derivative  filtered
    clearvars df
    df  = smoothdata(diff(yPointsSmo),'gaussian',50);
    df = [df(1);df];
    subplot(513); plot(df,'.r'); legend(strcat('df filtered--',Methods(ki)));
    subplot(514); plot(diff(diff(yPointsSmo)),'.g'); legend(strcat('d2f unfiltered--',Methods(ki)));
    % second derivative filtered
    clearvars d2f
    d2f = diff(smoothdata(diff(yPointsSmo),'gaussian',50));
    d2f = [d2f(1); d2f; d2f(end)];
    subplot(515); plot(d2f,'.g'); legend(strcat('d2f filtered--',Methods(ki)));
end

%% cropping according to the first derivative
marg = 0.25; % defines the margin for the first derivative
marg2 = 0.005; % defines the margin for the second derivative


plot_Image(RM_clean_vec,size(RM_clean_vec,2),'RM'); % plot the original RM
figure; plot(df,'.r'); legend(strcat('df filtered--',Methods(ki))); hold on; % first derivative
hline = refline([0 marg]); hline.Color = 'g';   % setting first upper threshold
hline = refline([0 -marg]); hline.Color = 'g';   % setting first lower threshold
xlim([0 size(RM_clean_vec,2)]);

figure; plot(d2f,'.r'); legend(strcat('d2f filtered--',Methods(ki))); hold on; % second derivative
hline = refline([0 marg2]); hline.Color = 'r';   % setting second upper threshold
hline = refline([0 -marg2]); hline.Color = 'r';   % setting second lower threshold
xlim([0 size(RM_clean_vec,2)]);

%% logic operator to connect the first and second derivative

% logical operator first derivative
abs_df = abs(df);
vec_df = logical(zeros(size(abs_df)));
vec_df(abs_df<= marg)  =  true;

% logical operator second derivative
abs_d2f = abs(d2f);
vec_d2f = logical(zeros(size(abs_d2f)));
vec_d2f(abs_d2f<= marg2)  =  true;

% AND connection
vec = logical(zeros(size(abs_df)));
vec((vec_df==1)&(vec_d2f==1)) = true;
figure; plot(vec); xlim([0 size(RM_clean_vec,2)]); hold on
ylim([-0.1 1.1]);


% sets a fluctuation short period to one
samp = 100;
[vec_filt] = setToOne(vec', samp);
plot(vec_filt,'g'); xlim([0 size(RM_clean_vec,2)]); hold on

% sets a fluctuation short period to zero
samp = 100;
[vec_filt] = setToZero(vec_filt, samp);
figure; title('final binary vector')
plot(vec_filt,'r'); xlim([0 size(RM_clean_vec,2)]); hold on; ylim([-0.1 1.1]);

%% separating the image in translationa and inplace

Im = [(repmat(vec_filt,16,1));RM_clean;repmat(vec_filt,16,1)];
plot_Image(Im,size(Im,2),'RM');
colormap('bone');

Im = [(repmat(vec_filt,16,1)); sspec_OLR; repmat(vec_filt,16,1)];
plot_Image(Im,size(Im,2),'MD');
colormap('bone');














return;




