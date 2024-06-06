function [spks, tlo, thi] = spkmtx(cdata, multi, Tlim)

%SPKMTX Produces a spike-train matrix from a cell array with spike times
%
%   [spks, tlo, thi] = spkmtx(cdata, [multi], [Tlim])
%
% spks has the same number of rows as cdata. Its numbers of columns is
% equal to the length of [tlo:thi].
%
% multi: 0 --> spks has only zeros and ones
%        1 --> spks may have values larger than 1, if more than 1
%              spike is found in a given time bin
%        default: 0
%
% Tlim: [tlo thi] time range. Values must be integers. Tlim ensures
%       two things, (1) that only spikes between tlo and thi are
%       included, and (2) that the length of spks is precisely
%       length(tlo:thi). 

% Emilio Salinas, July 2009, rev May 2014

 %
 % set default value of flag for multi-unit trains
 %
 if nargin < 2 | isempty(multi)
     multi = 0;
 elseif ~isscalar(multi)
     error('Parameter multi must be 0 or 1')
 elseif all(multi ~= [0 1])
     error('Parameter multi must be 0 or 1')
 end

 Nr = length(cdata);

 %
 % which spikes should be included?
 %
 if nargin < 3 | isempty(Tlim)
	 % find min and max spike times, to create a matrix that fits all
	 % the spikes
     tlo = Inf;
     thi = -Inf;
     for j=1:Nr
         svec = cdata{j};
         if ~isempty(svec)
             tlo = min(tlo, min(svec));
             thi = max(thi, max(svec));
         end
     end
 else
     % or set the size of the matrix as requested and eliminate
     % extraneous spikes
     tlo = Tlim(1);
     thi = Tlim(2);
     if thi <= tlo
         error('Limits Tlim = [Tlo Thi] must be increasing integers')
     elseif any(isinf(Tlim))
         error(['Limits Tlim = [Tlo Thi] cannot be infinite' char(10) ...
                'Leave Tlim empty to include all spikes'])
     end
     for j=1:Nr
         svec = cdata{j};
         ii = (svec < tlo) | (svec > thi);
         cdata{j} = svec(~ii);
     end
 end

 %
 % create the matrix and fill it assuming a 1 ms resolution
 %
 tlo = floor(tlo);
 thi = ceil(thi);
 ii = tlo:thi;
 Nc = length(ii);
 spks = zeros(Nr,Nc);
 
 if multi 
	 % data are from multiple neurons, so multiple spikes may appear
	 % in a given time bin; this is slow
     for j=1:Nr
         svec = cdata{j};
         N1 = length(svec);
         if N1 > 0
             ic = round(svec - tlo) + 1;
             for jj=1:N1
                 spks(j,ic(jj)) = spks(j,ic(jj)) + 1;
             end
         end
     end
 else
	 % data are from a single neuron, which either fires 1 spike or
	 % does not, in each bin
     for j=1:Nr
         svec = cdata{j};
         if ~isempty(svec)
             ic = round(svec - tlo) + 1;
             spks(j,ic) = 1;
         end
     end
 end
 
