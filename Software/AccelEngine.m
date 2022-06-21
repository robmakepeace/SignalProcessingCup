

classdef AccelEngine < handle
    properties
        NFFT;
        FS;
        TOLERANCE_BPM;
    end       
        
    methods
        
        function accel_engine = AccelEngine(NFFT, FS, TOLERANCE_BPM)
            accel_engine.NFFT = NFFT;
            accel_engine.FS = FS;
            accel_engine.TOLERANCE_BPM = TOLERANCE_BPM;
        end
        
        
        function candidates = compute(self, candidates, acc_combined_peaks)
                   
            if numel(acc_combined_peaks) > 0

                for i = 1:numel(candidates)
                    diffs = acc_combined_peaks - candidates(i).loc;
                    min_diff = min(abs(diffs));
                    min_diff_bpm = min_diff*60*self.FS/self.NFFT;
                    
                    % straight line increase to 1 at BPM_DIFF_TOLERANCE diff
                    if min_diff_bpm > self.TOLERANCE_BPM
                        candidates(i).Ca = 1;
                    else
                        candidates(i).Ca = min_diff_bpm/self.TOLERANCE_BPM;
                    end
                    
                    
                end
            else
                for i = 1:numel(candidates)
                    candidates(i).Ca = 1;
                end
            end
            % return candidates
        end
        
        
    end
    
    methods (Static)
        
       
    end
end