

classdef AccPeakProcessor < handle
    properties
        FS;
        NFFT;
        Kc2_ACC;
    end       
        
    methods
        
        function acc_peak_processor = AccPeakProcessor(NFFT, FS, Kc2_ACC)
            acc_peak_processor.FS = FS;
            acc_peak_processor.NFFT = NFFT;
            acc_peak_processor.Kc2_ACC = Kc2_ACC;
        end
        
        % Combine into one list of peak locations. Also add in missing harmonics
        function processed_peaks = process(self, acc1_peaks, acc2_peaks, acc3_peaks)
           
            combined_peaks = [];
            combined_peak_powers = [];
            if numel(acc1_peaks) > 0
                combined_peaks = [combined_peaks, acc1_peaks(:).loc];
                combined_peak_powers = [combined_peak_powers, acc1_peaks(:).power];
            end
            if numel(acc2_peaks) > 0
                combined_peaks = [combined_peaks, acc2_peaks(:).loc];
                combined_peak_powers = [combined_peak_powers, acc2_peaks(:).power];
            end
            if numel(acc3_peaks) > 0
                combined_peaks = [combined_peaks, acc3_peaks(:).loc];
                combined_peak_powers = [combined_peak_powers, acc3_peaks(:).power];
            end
            
            
            % for each acc peak, determine whether another peak at its (2nd or 3rd) harmonic is also present. If not, invalidate the peak
            valid = zeros(1, numel(combined_peaks));
            for i = 1:numel(combined_peaks)
                for j = (i+1):numel(combined_peaks)
                    ratio = max(combined_peaks(i), combined_peaks(j)) / min(combined_peaks(i), combined_peaks(j));
                    diff = [2, 3] - ratio;
                    if min(abs(diff)) < 0.05
                        valid(i) = 1;
                        valid(j) = 1;
                    end
                end
            end
            
            % but if peak has very high power, do not invalidate it even if it has no harmonics
            % put these peaks into separate list, add them back later (ie do not add their harmonics into the mix..)
            add_back = [];
            for i = 1:numel(combined_peaks)
                if combined_peak_powers(i) > self.Kc2_ACC && ~valid(i)
                    add_back = [add_back, combined_peaks(i)];
                end
            end
            combined_peak_powers = []; % invalidate this vector to prevent potential bugs..
            
            temp = [];
            for i = 1:numel(valid)
                if valid(i) 
                    temp = [temp, combined_peaks(i)];
                end
            end
            combined_peaks = temp;
                        
            % for each peak, add in its 2nd, 3rd harmonic 
            TOLERANCE_HARM = ceil(3*self.NFFT/(60*self.FS));  % 3bpm tolerance
            added_harmonics = [];
            for i = 1:numel(combined_peaks)
                for j = 1:3
                    harmonic = combined_peaks(i) * j;
                    % check if harmonic is already in peak list
                    temp = combined_peaks - harmonic;
                    temp = sort(abs(temp));
                    if temp(1) > TOLERANCE_HARM
                        added_harmonics = [added_harmonics, harmonic];
                    end
                end
            end
            combined_peaks = [combined_peaks, added_harmonics];
            
            % add back peaks
            combined_peaks = [combined_peaks, add_back];
            
            % remove double-ups
            TOLERANCE_DOUBLE = ceil(1*self.NFFT/(60*self.FS));  % 1bpm tolerance
            combined_peaks = self.remove_double_ups(combined_peaks, TOLERANCE_DOUBLE);
            
            
            % remove peaks that are too low bpm or too high
            LOWER_LIM = 20*self.NFFT/(60*self.FS);  % 20 bpm
            UPPER_LIM = 400*self.NFFT/(60*self.FS);
            temp = [];
            for i = 1:numel(combined_peaks)
                if combined_peaks(i) > LOWER_LIM && combined_peaks(i) < UPPER_LIM
                    temp = [temp, combined_peaks(i)];
                end
            end
            
            processed_peaks = temp;
            % return processed_peaks
        end
        
    end
    
    methods (Static)

        function list = remove_double_ups(list, tolerance) 
            list = sort(list);    % sort in ascending order
            temp1 = zeros(1, numel(list));
            temp1(1) = -Inf;

            temp1(2:end) = list(1:end-1);
            temp2 = list - temp1;

            j = 1;
            temp3 = [];
            for i = 1:numel(temp2)
                if temp2(i) > tolerance
                    temp3(j) = list(i);
                    j = j + 1;
                end
            end
            list = temp3;
        end
       
    end
end