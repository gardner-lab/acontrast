function [FEATURE_MAT,RAW_FEATURES,FEATURE_NAME,CONSENSUS,F,T]=aconstrast_contour(SIGNAL,FS,varargin)
%Computes acoustic contrast scores using the contour transform.
%
%	[FEATURE_MAT,RAW_FEATURES,FEATURE_NAME,CONSENSUS,SONOGRAMS,F,T]=aconstrast_contour(SIGNAL,FS,varargin)
%
%	SIGNAL
%	audio signal (double)
%
%	FS
%	sampling rate
%
%	the function has the following outputs:
%
%	FEATURE_MAT
%	regression coefficient matrix (sum absolute values across rows to get score, see demo for example)
%
%	RAW_FEATURES
%	raw contour features (typically not used)
%
%	FEATURE_NAMES
%	parameters used for each row of the RAW_FEATURES matrix
%
%	CONSENSUS
%	consensus contour image used as the basis for RAW_FEATURES
%	
%	F
%	frequency vector for CONSENSUS
%
%	T
%	time vector for CONSENSUS
%
%See also acontrast_demo.m, acontrast_deltacoef.m, acontrast_envelope.m

if nargin<1 | isempty(SIGNAL)
	error('Need signal to continue');
end

if nargin<2 |  isempty(FS)
	error('Need sampling rate to continue');
end

%% PARAMETERS (PASSED AS PARAMETER/VALUE PAIRS) 

use_band=[3e3 10e3]; % band to compute features in
pow_weight=0; % enables power weighting
pow_scale='log'; % power weighting (log or linear)
len = 23.2; % window length (in ms)
overlap = 22.8; % window overlap (in ms)
angle_list = (pi/8:pi/8:pi) + pi/8; % compute contours in all these angle_list - recommended not to change.
timescale_list = 0.5:0.2:2.2; % time scales in milliseconds for this analysis
clength_threshold = 95; % keep only contours longer than this length percentile 98 or 99 recommended.
norm_amp=1; % normalize amplitude to [-1,1]
filtering=300; % high pass filtering of mic signal
regression_timescale=.005; % in seconds, computes smooth derivative through windowed regression

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
		case 'angle_list'
			angle_list=varargin{i+1};
		case 'timescale_list'
			timescale_list=varargin{i+1};
		case 'clength_threshold'
			clength_threshold=varargin{i+1};
		case 'filtering'
			filtering=varargin{i+1};
		case 'norm_amp'
			norm_amp=varargin{i+1};
		case 'pow_weight'
			pow_weight=varargin{i+1};
		case 'use_band'
			use_band=varargin{i+1};
		case 'pow_scale'
			pow_scale=varargin{i+1};
		case 'regression_timescale'
			regression_timescale=varargin{i+1};
	end
end

[CONSENSUS,F,T,auditory_contour,SONOGRAMS]=acontour(SIGNAL,FS,'pow_weight',pow_weight,'norm_amp',1,'filtering',filtering,...
	'clength_threshold',clength_threshold,'len',len,'overlap',overlap,'timescale_list',timescale_list,'angle_list',angle_list);

minf=max(find(F<=use_band(1)));
maxf=min(find(F>=use_band(2)));

if isempty(minf), minf=1; end
if isempty(maxf), maxf=length(F); end

% form feature matrix using the auditory contours

nangles=length(angle_list);
ntimescales=length(timescale_list);

disp('Tiling feature matrix...');

RAW_FEATURES=zeros(nangles*ntimescales,length(T));

weights=SONOGRAMS;

% get weighting from sonograms

if pow_weight & strcmp(lower(pow_scale(1:3)),'log')	
	for i=1:length(weights)
		weights{i}=log(abs(weights{i}));
	end
elseif pow_weight & strcmp(lower(pow_scale(1:3)),'lin')
	for i=1:length(weights)
		weights{i}=abs(weights);
	end
else
	for i=1:length(weights)
		weights{i}=ones(size(weights{i}));
	end
end

counter=1;
for i=1:ntimescales
	for j=1:nangles
		RAW_FEATURES(counter,:)=full(sum(auditory_contour{i}{j}(minf:maxf,:).*weights{i}(minf:maxf,:)));
		FEATURE_NAME{counter}=sprintf('Timescale %5.3f (ms) & angle %5.3f (rads)',timescale_list(i),angle_list(j));
		counter=counter+1;
	end
end

score_fs=1./(T(2)-T(1)); % get the frame rate of the feature vector
tmp=acontrast_deltacoef(RAW_FEATURES,round(score_fs.*regression_timescale)); % yoonseob used 10 msec regression window	

% interpolate back to original signal space

original_time=[1:length(SIGNAL)]/FS;

disp('Interpolating features...');

FEATURE_MAT=zeros(nangles*ntimescales,length(original_time));
for i=1:size(RAW_FEATURES,1)
	FEATURE_MAT(i,:)=interp1(T,tmp(i,:),original_time,'nearest');
end


