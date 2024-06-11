function cmap = polarmap(varargin)
%POLARMAP Polarized color map
%	POLARMAP applies a "polarized" blue-white-red colormap to current figure,
%	and adjusts the color axis limits to be centered to zero.
%
%	POLARMAP(M) fixes the number of colors to M (default is 64).
%
%	POLARMAP(MAP) applies linear shading to white to the center of colormap
%	MAP which can be any of existing colormaps (an Mx3 matrix of RGB).
%
%	POLARMAP(MAP,C) uses exponent C to modify the shading contrast. Default 
%	is C = 1 for linear shading. Use C = 2 to strengthen the shading, or 
%	C = 0.5 to attenuate it.
%
%	C=POLARMAP(...) returns an M-by-3 matrix containing the colormap, that 
%	can be used with COLORMAP function like other colormaps.
%
%	Examples:
%		pcolor(peaks), shading interp
%		polarmap, colorbar
%
%	then try the following
%		polarmap(jet,0.5)
%
%	Note the polar shading has no real interest with colormaps that include
%	white color as one of the extremes (like GRAY, BONE, HOT, ...).
%
%	See also JET, HSV, COPPER, SPRING, SUMMER, WINTER, COOL, COLORMAP, RGBPLOT.

% default parameters
m = 64;	% number of colors
c = 1;	% exponent of shading factor (1 = linear)
if nargin > 0
	if ~isnumeric(varargin{1}) | (size(varargin{1},2) ~= 3 & ~isscalar(varargin{1}))
		error('First argument must be numeric: scalar M or Mx3 color matrix');
	end
	if isscalar(varargin{1})
		m = varargin{1};
	end
end
if nargin > 0 & size(varargin{1},2) == 3
		map = varargin{1};
		m = size(map,1);
else
	map = bluered(m);
end
if nargin > 1 & isscalar(varargin{2})
	c = varargin{2};
end
% linear shading from min/max (colormap value) to center (white)
r = repmat(abs(linspace(1,-1,m)).^c,[3,1])';
map = map.*r + 1 - r;
if nargout > 0
	cmap = map;
else
	colormap(map)
	caxis([-1,1]*max(abs(caxis)))
	% Note: this fixes color axis to manual mode...
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function map = bluered(m)
if mod(m,2)
	z = [0,0,0];
	m2 = floor(m/2);
else
	z = zeros([0,3]);
	m2 = m/2;
end
map = [repmat([0,0,1],[m2,1]);z;repmat([1,0,0],[m2,1])];
