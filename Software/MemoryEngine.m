

classdef MemoryEngine < handle 
    
    properties
        last_outcome;
    end
        
    methods
       
        function memory_engine = MemoryEngine()
            memory_engine.invalid_last_outcome();
        end
        
        function candidates = compute(self, candidates)
            
            if numel(candidates) > 0 
                if self.last_outcome.loc ~= -1
                
                    for i = 1:numel(candidates)
                        D = abs((self.last_outcome.loc - candidates(i).loc)/candidates(i).loc);
                        candidates(i).Cm = min(1, max((1-D/0.3), 0));
                    end
                else
                    for i = 1:numel(candidates)
                        candidates(i).Cm = 1;
                    end
                end
            end
            
            % return candidates
        end
        
        
        function update_last_outcome(self, candidate)
            self.last_outcome = candidate;
        end
        
        
        % not to be called if self.last_outcome is still invalid
        function candidate = repeat_last_outcome(self)
            candidate = self.last_outcome;
        end
        
        
        function loc = get_last_outcome_loc(self)
            loc = self.last_outcome.loc;
        end
        
        
        function invalid_last_outcome(self)
            self.last_outcome.loc = -1;
        end
        
    end
end