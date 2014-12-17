
% visualization parameters

font_size=12; % self-explanatory
disp_band=[0 12e3]; % spectrogram band to display
clipping=0; % clipping for spectrogram
imscale=15; % image scaling parameter
regression_timescale=.005; % sets the smoothness of the derivative calculation

% Load the sample sound

[signal, fs] = wavread('finchdoublet.wav');
signal_t=[1:length(signal)]/fs;

tic
[feature_matrix,raw_features,feature_names,consensus,f,t]=acontrast_contour(signal,fs,'regression_timescale',regression_timescale);
[envelope_score,envelope,sonogram,sono_f,sono_t]=acontrast_envelope(signal,fs,'regression_timescale',regression_timescale);
contour_score=sum(abs(feature_matrix));
toc

%% Visualization

figure
h(1) = subplot(4,1,1:2);
imagesc(sono_t,sono_f,max(log(abs(sonogram)+imscale),clipping));
ylim([disp_band])
ylabel('Frequency [Hz]','FontSize', font_size);
set(gca,'XTick',[]);
axis xy;
title('Sonogram','FontSize', font_size);
colormap(hot)

h(2) = subplot(4,1,3);
plot(signal_t,contour_score);
ylabel('Cont.','FontSize', font_size);
set(gca,'XTick',[]);
axis xy;

h(3)=subplot(4,1,4);
plot(signal_t,abs(envelope_score));
ylabel('Env.','FontSize',font_size);
xlabel('Time (in s)','FontSize',font_size);
linkaxes(h,'x')




