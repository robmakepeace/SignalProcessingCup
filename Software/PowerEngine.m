

classdef PowerEngine < handle     
        
    methods
       
        function candidates = compute(self, candidates)
            
            if numel(candidates) > 0
                powers = [candidates(:).power];
                max_power = max(powers);

                for i = 1:numel(candidates)
                    candidates(i).Cp = candidates(i).power / max_power;                
                end
            end
            
            % return candidates
        end 
        
    end
end