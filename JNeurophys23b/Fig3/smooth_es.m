function [ys] = smooth_es(y, sd, yx)

%SMOOTH Simple smoothing with a Gaussian kernel
%
%    [ys] = smooth(y, sd, [yx])
%
% Convolves the vector y with a Gaussian function of width sd, where
% the spacing of values along the x axis is assumed to be 1. If y is a
% matrix then each column is smoothed. The parameter yx controls what
% happens at the edges:
%
%    yx = 0 --> nothing is done
%
%    yx = 1 --> the last points are reflected across the boundary,
%               which produces a correction; this is the default

% Emilio Salinas, June 1998; rev July 2009

 if nargin < 3 | isempty(yx)
     yx = 1;
 else
     yx = 1*(yx > 0);
 end
    
 %
 % check input dimensions
 %
 if size(y,1) == 1
     y = y';
 end
 [Nr Nc] = size(y);

 %
 % make x-grid for Gaussian kernel and compute coefficients
 %
 sd2 = sd*sd;
 Nxs = 4*sd;
 if Nxs < 1
     Nxs = 1;
 else
     Nxs = ceil(Nxs);
 end
 xs = [-Nxs:Nxs]';
 g  = exp(-0.5*xs.^2/sd2); 
 g  = g/sum(g);            

 %
 % width of Gaussian is larger than data range; use spike counts
 % instead
 %
 if Nxs > Nr
     warning('smoothing width is too large')
 end

 %
 % to correct edge effects, add a mirror image of the first and 
 % last Nxs points, or zeroes if no correction is requested
 %
 if yx > 0
     lo_edge = yx*y(Nxs:-1:1,:);
     hi_edge = yx*y(end:-1:end-Nxs+1,:);
 else
     lo_edge = zeros(Nxs,Nc);
     hi_edge = zeros(Nxs,Nc);
 end
 y = [lo_edge; y; hi_edge];

 %
 % now convolve the vectors
 %
 [Nr Nc] = size(y);
 ys = zeros(Nr,Nc);
 for j=1:Nc
     yst = conv(y(:,j),g);        % convolve the 'spikes' with the kernel
     Nc1 = length(yst);           % length is size(g) + size(y) - 1
     yst = yst(Nxs+1:Nc1-Nxs);    % resize smoothed vector to length(y) 
     ys(:,j) = yst;    
 end

 %
 % finally, get rid of the extra points that were added to correct
 % for edge effects
 %
 ys(1:Nxs,:) = [];
 ys(end-Nxs+1:end,:) = [];





