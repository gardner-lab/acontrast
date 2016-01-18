function [ACONTRAST_SCORE,ENVELOPE,SONOGRAM,F,T]=aconstrast_envelope(SIGNAL,FS,varargin)
%Computes acoustic contrast scores using the derivative of the amplitude envelope. The
%envelope is computed by taking the summed abs value of the Gabor transform (i.e. Gaussian
%windowed spectrogram) across a frequency range of interest (typically 1e3-12e3 for
%birdsong).  This function accepts a number of parameters, see either the Github wiki
%or edit the function itself for parameter descriptions.  Note that all parameters are passed
%as parameter/value pairs (see demo).
%
%	[ACONTRAST_SCORE,ENVELOPE,SONOGRAM,F,T]=aconstrast_envelope(SIGNAL,FS,varargin)
%
%	SIGNAL
%	audio signal (double)
%
%	FS
%	sampling rate
%
%	the function has the following outputs:
%
%	ACONTRAST_SCORE
%	abs derivative of the amplitude envelope
%
%	ENVELOPE
%	the amplitude envelope
%
%	SONOGRAM
%	sonogram used to compute the amplitude envelope
%
%	F
%	frequency vector for sonogram
%
%	T
%	time vector for sonogram
%
%See also acontrast_demo.m, acontrast_deltacoef.m, acontrast_contour.m

if nargin<1 | isempty(SIGNAL)
	error('Need signal to continue');
end

if nargin<2 |  isempty(FS)
	error('Need sampling rate to continue');
end

%% BEGIN USER PARAMETERS

use_band=[1e3 12e3]; % band to compute features in
pow_scale='log'; % power weighting (log or linear)
len = 23; % window length (in ms)
overlap = 22.9; % window overlap (in ms)
timescale= 1; % window timescale (in ms)
norm_amp=1; % normalize amplitude to [-1,1]
filtering=300; % high pass filtering of mic signal
regression_timescale=.005; % in seconds, computes smooth derivative through windowed regression
zeropad=0; % zeropad for the signal (in ms, set to 0 to automatically set to spectrogram len/2, empty for none)

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

% convert parameters to samples

len=round((len/1e3)*FS);
overlap=round((overlap/1e3)*FS);
timescale=(timescale/1e3)*FS;

if zeropad==0
	zeropad=round(len/2);
elseif ~isempty(zeropad)
	zeropad=round((zeropad/1e3)*FS);
end

original_time=[1:length(SIGNAL)]/FS;

% zeropad if necessary

if ~isempty(zeropad)
	SIGNAL=[zeros(zeropad,1);SIGNAL(:);zeros(zeropad,1)];
	disp(['Zero pad: ' num2str(zeropad/FS) ' S']);
end

% compute spectrogram w/ Gauss window

win_t=-len/2+.5:len/2-.5;
window=exp(-(win_t/timescale).^2);

[SONOGRAM,F,T]=spectrogram(SIGNAL,window,overlap,[],FS);

if ~isempty(zeropad)
	T=T-zeropad/FS;
end

minf=max(find(F<=use_band(1)));
maxf=min(find(F>=use_band(2)));

if isempty(minf), minf=1; end
if isempty(maxf), maxf=length(F); end

% select amplitude weighting (log, linear or z-score)

switch lower(pow_scale)
	case 'log'
		ENVELOPE=sum(abs(log(SONOGRAM(minf:maxf,:))));
	case 'lin'
		ENVELOPE=sum(abs(SONOGRAM(minf:maxf,:)));
	case 'z'
		ENVELOPE=sum(abs(zscore(SONOGRAM(minf:maxf,:),[],2)));
	otherwise
		error('Did not understand power scale argument (options are %s, %s and %s)','log','lin','z');
end

score_fs=1./(T(2)-T(1));
tmp=acontrast_deltacoef(ENVELOPE,round(score_fs*regression_timescale));

ACONTRAST_SCORE=interp1(T,tmp,original_time,'nearest');
