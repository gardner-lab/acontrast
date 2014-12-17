function DELTAC=acontrast_deltacoef(FEATURE,WIN,PAD)
% acontrast_deltacoef computes the delta coefficients for cepstral coefficient matrix
% FEATURE (coefficients x frames)
%
%	acontrast_deltacoef(FEATURE,WIN,PAD)
%
%	FEATURE
%	feature matrix (features x samples)
%
%	WIN
%	Window size 
%
%	PAD
%	Zero pad to return vector same size as number of frames in 
%	FEATURE (default: 1)
%
% simple function that computes the delta coefficients

% if padding enabled, zero pad to return vector of same size

if nargin<3 | isempty(PAD), PAD=1; end
if nargin<2 | isempty(WIN), WIN=2; end
if nargin<1, error('Need feature matrix to continue'); end

if PAD
	FEATURE=[ones(size(FEATURE,1),WIN).*repmat(FEATURE(:,1),[1 WIN]) ...
	       	FEATURE ...
		ones(size(FEATURE,1),WIN).*repmat(FEATURE(:,end),[1 WIN])];
end

WIN=round(WIN);

[rows,columns]=size(FEATURE);

% lose the edges via the window

DELTAC=zeros(rows,columns-(2*(WIN+1)));

for i=WIN+1:columns-(WIN)
	
	% compute the regression coefficient

	deltanum=sum(repmat(1:WIN,[rows 1]).*(FEATURE(:,i+1:i+WIN)-FEATURE(:,i-WIN:i-1)),2);
	deltaden=2*sum([1:WIN].^2);
	DELTAC(:,i-(WIN))=deltanum./deltaden;

end
