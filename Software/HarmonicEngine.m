

classdef HarmonicEngine < handle

    methods
        
        % Ch will be 1 unless candidate is the 2nd/3rd/4th harmonic of another candidate,
        function candidates = compute(self, candidates)
            
            for i = 1:numel(candidates)
                candidates(i).Ch = 1;
            end
        
            for i = 1:numel(candidates)
                for j = (i+1):numel(candidates)
                    higher_loc_index = i;
                    lower_loc_index = j;
                    if candidates(j).loc > candidates(i).loc
                        higher_loc_index = j;
                        lower_loc_index = i;
                    end

                    harmonic_ratio = candidates(higher_loc_index).loc/candidates(lower_loc_index).loc;
                    diff = [2, 3, 4] - harmonic_ratio;
                    min_diff = min(abs(diff));

                    % give score of 0 at anything less than 0.05, then straight line increase to 1 at 0.1
                    score = max(0, min((min_diff/0.05 - 1), 1));
                    candidates(higher_loc_index).Ch = min(candidates(higher_loc_index).Ch, score);
                end
            end
            
            % return candidates;
        end
        
    end
    
end