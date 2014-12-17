function [ACONTRAST_SCORE,ENVELOPE,SONOGRAM,F,T]=aconstrast_envelope(SIGNAL,FS,varargin)
%	
%
%
%%

use_band=[3e3 10e3]; % band to compute features in
pow_scale='log'; % power weighting (log or linear)
len = 23; % window length (in ms)
overlap = 22.6; % window overlap (in ms)
timescale= 2; % window timescale (in ms)
clength_threshold = 95; % keep only contours longer than this length percentile 98 or 99 recommended.
norm_amp=1; % normalize amplitude to [-1,1]
filtering=300; % high pass filtering of mic signal
regression_timescale=.005; % in seconds, computes smooth derivative through windowed regression
zeropad=0;

%% END USER PARAMETERS

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'len'
			len=varargin{i+1};
		case 'overlap'
			overlap=varargin{i+1};
		case 'timescale'
			timescale=varargin{i+1};
		case 'filtering'
			filtering=varargin{i+1};
		case 'norm_amp'
			norm_amp=varargin{i+1};
		case 'use_band'
			use_band=varargin{i+1};
		case 'pow_scale'
			pow_scale=varargin{i+1};
		case 'regression_timescale'
			regression_timescale=varargin{i+1};
		case 'zeropad'
			zeropad=varargin{i+1};
	end
end

len=round((len/1e3)*FS);
overlap=round((overlap/1e3)*FS);
timescale=(timescale/1e3)*FS;

if zeropad==0
	zeropad=round(len/2);
end

if ~isempty(zeropad)
	SIGNAL=[zeros(zeropad,1);SIGNAL(:);zeros(zeropad,1)];
	disp(['Zero pad: ' num2str(zeropad/FS) ' S']);
end

win_t=-len/2+.5:len/2-.5;
window=exp(-(win_t/timescale).^2);

[SONOGRAM,F,T]=spectrogram(SIGNAL,len,overlap,[],FS);

if ~isempty(zeropad)
	T=T-zeropad/FS;
end

minf=max(find(F<=use_band(1)));
maxf=min(find(F>=use_band(2)));

if isempty(minf), minf=1; end
if isempty(maxf), maxf=length(F); end

switch lower(pow_scale)
	case 'log'
		ENVELOPE=sum(abs(log(SONOGRAM(minf:maxf,:))));
	case 'lin'	
		ENVELOPE=sum(abs(SONOGRAM(minf:maxf,:)));
	case 'z'
		ENVELOPE=sum(abs(zscore(SONOGRAM(minf:maxf,:),[],2)));
end

score_fs=1./(T(2)-T(1));
tmp=acontrast_deltacoef(ENVELOPE,round(score_fs*regression_timescale));

original_time=[1:length(SIGNAL)]/FS;
ACONTRAST_SCORE=interp1(T,tmp,original_time,'nearest');







