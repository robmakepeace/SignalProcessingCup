

classdef PeakFinder < handle
    
    properties
        Kc;
        max_num_peaks;
        min_peak_loc;
    end       
        
    methods
        
        function peak_finder = PeakFinder(Kc, max_num_peaks, min_peak_loc)
            peak_finder.Kc = Kc;
            peak_finder.max_num_peaks = max_num_peaks;
            peak_finder.min_peak_loc = min_peak_loc;
        end
        
        function candidates = find_peaks(self, power_spectrum)
            
            [peak_values, peak_locs] = self.findpeaks(power_spectrum, 'minpeakheight', self.Kc, 'minpeakdistance', 1, 'sortstr', 'descend');
                        
            % remove peaks that are at low frequencies
            if numel(peak_locs) > 0
                for i = wrev(1:numel(peak_values))
                    if peak_locs(i) < self.min_peak_loc
                        peak_locs(i) = [];  % hacky matlab way of deleting vector element
                        peak_values(i) = []; 
                    end
                end
            end
            
            if numel(peak_locs) == 0
                candidates = [];
            end
            
            for i = 1:min(self.max_num_peaks, numel(peak_locs))
                candidate.power = peak_values(i);
                candidate.loc = peak_locs(i);
                
                candidates(i) = candidate;
            end
            
            % return candidates
        end
        
               
    end
    
    methods (Static)
        
        % copied from matlab system findpeaks function but with warnings disabled (I couldn't find a way to disable warnings..)
        function [pks,locs] = findpeaks(X,varargin)

            error(nargchk(1,11,nargin,'struct'));

            validateattributes(X,{'numeric'},{'nonempty','real','vector'},...
                'findpeaks','X');

            %#function dspopts.findpeaks
            hopts = uddpvparse('dspopts.findpeaks',varargin{:});

            Ph  = hopts.MinPeakHeight;
            Pd  = hopts.MinPeakDistance;
            Th  = hopts.Threshold;
            Np  = hopts.NPeaks;
            Str = hopts.SortStr;

            pks = [];
            locs = [];


            M = numel(X);

            if (M < 3)
                datamsgid = generatemsgid('emptyDataSet');
                error(datamsgid,'Data set must contain at least 3 samples.');
            else
                Indx = find(X > Ph);
                if(isempty(Indx))
                    %mphmsgid = generatemsgid('largeMinPeakHeight');
                    %warning(mphmsgid,'Invalid MinPeakHeight. There are no data points greater than MinPeakHeight.');
                else
                    % validate value of Pd and set default values for Pd and Np
                    [Pd,Np] = setvalues(Pd,Np,M);
                    if(Pd >= M)
                        pdmsgid = generatemsgid('largeMinPeakDistance');
                        error(pdmsgid,'Invalid MinPeakDistance. Set MinPeakDistance as an integer in the range between 1 and %s.',...
                            num2str(M));
                    else
                        [pks,locs] =getpeaks(X,Indx,Pd,Th,Np);
                    end
                end
            end

            if isempty(pks)
                %npmsgid = generatemsgid('noPeaks');
                %warning(npmsgid,'No peaks found.')
            elseif(~strcmp(Str,'none'))
                [pks,s]  = sort(pks,Str);
                if(~isempty(locs))
                    locs = locs(s);
                end
            end
        end

    end
end



% copied from matlab system findpeaks function but with warnings disabled (I couldn't find a way to disable warnings..)
%--------------------------------------------------------------------------
% private static functions
function [pks,locs] =getpeaks(Data,Indx,Pd,Th,Np)
    % This function finds peaks in data set Data whose index set Indx is
    % disjoint. Some data points were removed from the original set through
    % preprocessing

    m = 0;                  % counter for peaks found
    L = numel(Indx);
    M = numel(Data);

    % Pre-allocate for speed
    pks  = zeros(1,Np);
    locs = zeros(1,Np);

    endindx = M;      % last point in the index set

    j = 0;
    % First, the "Pd" neighbors, on either side, of the current data point are
    % found. Then the current data point is compared with these neighbors to
    % determine whether it is a peak.

    while (j < L) && (m < Np)
        j = j+1;

        currentIdx = Indx(j);

        % leftmost neighbor index
        endL = max(1,currentIdx - Pd);

        % Update the leftmost neighbor index if there is a peak within "Pd"
        % neighbors of leftmost index
        if(m > 0)
            prevPeakBoundL = min([locs(m)+Pd, endindx-1]);
            if currentIdx < prevPeakBoundL
                k = find(Indx(j+1:end),1,'first');
                if ~isempty(k)
                    j = j+k;
                    currentIdx = Indx(j);
                else
                    break;
                end
            end
        end

        % rightmost neighbor index
        endR = min(currentIdx + Pd,endindx);

        % create neighbor data set
        temp = Data(endL:endR);

        % set current data point to -Inf in the neighbor data set
        temp(currentIdx-endL+1) = -Inf;

        % Compare current data point with all neighbors
        if(all((Data(currentIdx) > temp+Th)) && (currentIdx ~=1)&& (currentIdx~=M))
            m = m+1;
            locs(m) = currentIdx;  % loctions of peaks
            pks(m)  = Data(currentIdx);  % peaks
        end
    end

    % return all peaks found
    if m~=0
        locs = locs(locs > 0);
        pks  = pks(1:length(locs));
    else
        locs = [];
        pks = [];
    end
end

function [Pd,Np] = setvalues(Pd,Np,L)

    if ~isempty(Pd) && (~isnumeric(Pd) || ~isscalar(Pd) ||any(rem(Pd,1)) || (Pd < 1))
        Nmsgid = generatemsgid('invalidMinPeakDistance');
        error(Nmsgid,'MinPeakDistance should be an integer greater than 0.');
    end

    if(isempty(Pd))
        Pd = 1;
    end

    if(isempty(Np))
        Np = L;
    end
end