
classdef CandidateAnalyser < handle
    properties
        
        % default parameters for algorithm
        BPF_Fc1 = 1;
        BPF_Fc2 = 6;    
        
        NFFT = 50000;
        WINDOW_N = 1000;
        
        Kc_PPG = 200;
        MAX_NUM_PEAKS_PPG = 4;
        Kc1_ACC = 0.5;
        Kc2_ACC = 10;
        MAX_NUM_PEAKS_ACC = 4;
        
        ACCEL_TOLERANCE_BPM = 10;
        FUND_SHARPER_BETA = 0.5;
        MIN_C_OVERALL = 0.6;
        TRACK_ENGINE_N_FRAMES = 10;
        MEMORY_ENGINE_REPEATS_ALLOWED = 1;
        SEARCH_RANGE_BPM = 20;
        
        VERIFICATION_Ct = 0.7;
        VERIFICATION_Cta = 0.9;
        VERIFICATION_Cth = 0.9;
        VERIFICATION_K = 4;

        Ca_weight = 1;
        Cp_weight = 1;
        Ch_weight = 1;
        Ct_weight = 1;
        Cm_weight = 2;
        
        % state variables used by algorithm
        FS;
        H_bpf;
        peak_finder_ppg;
        peak_finder_acc;
        acc_peak_processor;
        power_engine;
        accel_engine;
        harmonic_engine;
        track_engine;
        memory_engine;
        
        output_i;
        memory_engine_repeated
        
        % debug variables
        candidate_frames;
        acc_peak_frames;
        %reference_output;
        
    end       
        
    methods
        
        
        function initialise(self, fs)
        
            self.FS = fs;
            
            % design 10th order butterworth filter
            d = fdesign.bandpass('N,F3dB1,F3dB2', 10, self.BPF_Fc1, self.BPF_Fc2, 125);
            self.H_bpf = design(d,'butter');
            
            % initialise component classes
            self.peak_finder_ppg = CandidateAnalyserPackage.PeakFinder(self.Kc_PPG, self.MAX_NUM_PEAKS_PPG, 0);
            self.peak_finder_acc = CandidateAnalyserPackage.PeakFinder(self.Kc1_ACC, self.MAX_NUM_PEAKS_ACC, 10*self.NFFT/(60*self.FS));
            self.acc_peak_processor = CandidateAnalyserPackage.AccPeakProcessor(self.NFFT, self.FS, self.Kc2_ACC);
            self.power_engine = CandidateAnalyserPackage.PowerEngine();
            self.accel_engine = CandidateAnalyserPackage.AccelEngine(self.NFFT, self.FS, self.ACCEL_TOLERANCE_BPM);
            self.harmonic_engine = CandidateAnalyserPackage.HarmonicEngine();
            self.track_engine = CandidateAnalyserPackage.TrackEngine(self.TRACK_ENGINE_N_FRAMES, self.NFFT, self.FS);
            self.memory_engine = CandidateAnalyserPackage.MemoryEngine();
            
            self.memory_engine_repeated = 0;
            self.output_i = 1;
            
            % reset debug variables
            self.candidate_frames = [];
            self.acc_peak_frames = [];

        end
        
        
        function output_bpm = compute_block(self, ppg_seg, acc_seg)
            
            % apply BPF to segments
            ppg_seg(1, :) = filter(self.H_bpf, ppg_seg(1, :));
            ppg_seg(2, :) = filter(self.H_bpf, ppg_seg(2, :));


            % computer power spectrum
            S_PPG1 = periodogram(ppg_seg(1, :), hanning(self.WINDOW_N), self.NFFT);
            S_PPG2 = periodogram(ppg_seg(2, :), hanning(self.WINDOW_N), self.NFFT);
            S_ACC1 = periodogram(acc_seg(1, :), hanning(self.WINDOW_N), self.NFFT);
            S_ACC2 = periodogram(acc_seg(2, :), hanning(self.WINDOW_N), self.NFFT);
            S_ACC3 = periodogram(acc_seg(3, :), hanning(self.WINDOW_N), self.NFFT);


            % find peaks
            candidates1 = self.peak_finder_ppg.find_peaks(S_PPG1);
            candidates2 = self.peak_finder_ppg.find_peaks(S_PPG2);
            acc1_peaks = self.peak_finder_acc.find_peaks(S_ACC1);
            acc2_peaks = self.peak_finder_acc.find_peaks(S_ACC2);
            acc3_peaks = self.peak_finder_acc.find_peaks(S_ACC3);

            % increase the power of each candidate by FUND_SHARPER_BETA * its 2nd harmonic
            candidates1 = self.sharpen_fund_freq_candidates(candidates1, S_PPG1, self.FUND_SHARPER_BETA);
            candidates2 = self.sharpen_fund_freq_candidates(candidates2, S_PPG2, self.FUND_SHARPER_BETA);


            % compute Cp - decouples Cp of the two candidate sets 
            candidates1 = self.power_engine.compute(candidates1);
            candidates2 = self.power_engine.compute(candidates2);

            % combine candidates from the two ppg channels
            candidates = [candidates1, candidates2];


            % combine acc peaks, add in missing harmonics and filters out bad peaks
            acc_combined_peaks = self.acc_peak_processor.process(acc1_peaks, acc2_peaks, acc3_peaks);


            % compute Ca (confidence accelerometer)
            candidates = self.accel_engine.compute(candidates, acc_combined_peaks);

            % compute Ch (condidence harmonic engine)
            candidates = self.harmonic_engine.compute(candidates);

            % compute Ct (confidence track engine)
            candidates = self.track_engine.compute(candidates);

            % compute Cm (confidence memory engine)
            candidates = self.memory_engine.compute(candidates);

            % compute C_overall (aggregate score)
            candidates = self.compute_C_overall(candidates);

            % decide the best candidate for output
            chosen_candidate = self.decide_candidate(candidates);

            % compute output bpm from candidate
            output_bpm = 60*self.FS*chosen_candidate.loc/self.NFFT;

            % for debug
            self.candidate_frames(self.output_i).loc = [candidates(:).loc];
            self.candidate_frames(self.output_i).bpm = [candidates(:).loc]*(60*self.FS/self.NFFT);
            self.candidate_frames(self.output_i).power = [candidates(:).power];
            self.candidate_frames(self.output_i).Cp = [candidates(:).Cp];
            self.candidate_frames(self.output_i).Ca = [candidates(:).Ca];
            self.candidate_frames(self.output_i).Ch = [candidates(:).Ch];
            self.candidate_frames(self.output_i).Ct = [candidates(:).Ct];
            self.candidate_frames(self.output_i).Ctp = [candidates(:).Ctp];
            self.candidate_frames(self.output_i).Cta = [candidates(:).Cta];
            self.candidate_frames(self.output_i).Cth = [candidates(:).Cth];
            self.candidate_frames(self.output_i).Cm = [candidates(:).Cm];
            self.candidate_frames(self.output_i).C_overall = [candidates(:).C_overall];

            self.acc_peak_frames(self.output_i).loc = acc_combined_peaks;
            self.acc_peak_frames(self.output_i).bpm = acc_combined_peaks*(60*self.FS/self.NFFT);


            % debug plots
            if self.output_i == 0
                figure()
                plot(S_PPG1(1:500))
                title('S_PPG1')

                figure()
                plot(S_PPG2(1:500))
                title('S_PPG2')

                figure()
                plot(S_ACC1(1:500))
                title('S_ACC1')

                figure()
                plot(S_ACC2(1:500))
                title('S_ACC2')

                figure()
                plot(S_ACC3(1:500))
                title('S_ACC3')

%                 if numel(self.reference_output) > 0
%                     self.reference_output(self.output_i)*self.NFFT/(60*self.FS)
%                 end

                temp_acc_peaks = [];
                if numel(acc1_peaks) > 0
                    temp_acc_peaks = [temp_acc_peaks, acc1_peaks(:).loc]
                end
                if numel(acc2_peaks) > 0
                    temp_acc_peaks = [temp_acc_peaks, acc2_peaks(:).loc]
                end
                if numel(acc3_peaks) > 0
                    temp_acc_peaks = [temp_acc_peaks, acc3_peaks(:).loc]
                end
                temp_acc_peaks

                acc_combined_peaks

                temp_candidates_locs = [candidates(:).loc]
                temp_candidates_Cp = [candidates(:).Cp]
                temp_candidates_Ca = [candidates(:).Ca]
                temp_candidates_Ch = [candidates(:).Ch]
                temp_candidates_Cm = [candidates(:).Cm]
                temp_candidates_Ct = [candidates(:).Ct]

            end
            
            % increment output_i as it is used by a few algorithm blocks
            self.output_i = self.output_i + 1;
        
        end
        
        
        % Debug function - not used for competition
        function output = compute_bpm(self, ppg_signal, acc_signal, fs)
            
            %ToolsPackage.do_spectrogram(ppg_signal(1,:), self.reference_output);
            %ToolsPackage.do_spectrogram(ppg_signal(2,:), self.reference_output);
            %ToolsPackage.do_spectrogram(ppg_signal(1,:)-ppg_signal(2,:), self.reference_output);
            %ToolsPackage.do_spectrogram(acc_signal(1,:), self.reference_output);
            %ToolsPackage.do_spectrogram(acc_signal(2,:), self.reference_output);
            %ToolsPackage.do_spectrogram(acc_signal(3,:), self.reference_output);
     
                        
            step = 2*fs;
            windowNb = (length(ppg_signal)-self.WINDOW_N)/step + 1; 
            self.initialise(fs);
            output = zeros(1, floor(windowNb));
            output_i = 1;
            
            for i = 1:windowNb
                curSegment = (i-1)*step+1 : (i-1)*step+self.WINDOW_N;
                output(output_i) = self.compute_block(ppg_signal(:, curSegment), acc_signal(:, curSegment));
                output_i = output_i + 1;
            end
                
                            
            self.plot_candidates(self.candidate_frames, self.acc_peak_frames, output);
            %self.plot_confidence(self.candidate_frames);
            %ToolsPackage.do_spectrogram(ppg_signal(1,:), output);
            %ToolsPackage.do_spectrogram(ppg_signal(2,:), output);
            %ToolsPackage.do_spectrogram(acc_signal(1,:), output);
            %ToolsPackage.do_spectrogram(acc_signal(2,:), output);
            %ToolsPackage.do_spectrogram(acc_signal(3,:), output);
            
            % Debug hack.
            temp = output;
            clear output;
            output.output = temp;
            output.analyser_optional_output = self.candidate_frames;
            
        end
        
        % computes C_overall for all candidates
        function candidates = compute_C_overall(self, candidates)
            for i = 1:numel(candidates)
                candidate = candidates(i);
                score = (self.Ca_weight*candidate.Ca + self.Cp_weight*candidate.Cp + self.Ch_weight*candidate.Ch + self.Cm_weight*candidate.Cm + self.Ct_weight*candidate.Ct)/(self.Ca_weight + self.Cp_weight + self.Ch_weight + self.Cm_weight + self.Ct_weight);
                candidates(i).C_overall = score;
            end
            % return candidates
        end
        
        
        function chosen_candidate = decide_candidate(self, candidates)
            
            % sort candidates on descending order for C_overall
            cand_cells = struct2cell(candidates);
            cand_size = size(cand_cells);
            cand_fields = fieldnames(candidates);
            cand_cells = reshape(cand_cells, cand_size(1), []);
            cand_cells = cand_cells';
            cand_cells = sortrows(cand_cells, -find(ismember(cand_fields, 'C_overall') ~= 0));
            cand_cells = reshape(cand_cells', cand_size);
            candidates = cell2struct(cand_cells, cand_fields, 1);
        
            if candidates(1).C_overall < self.MIN_C_OVERALL && self.memory_engine_repeated < self.MEMORY_ENGINE_REPEATS_ALLOWED && self.output_i ~= 1
                chosen_candidate = self.memory_engine.repeat_last_outcome();
                self.memory_engine_repeated = self.memory_engine_repeated + 1;
            else
               
                % invalidate all candidates that are outside search range based on previous output location
                if self.output_i ~= 1
                    
                    search_range_valid = ones(1, numel(candidates));
                    
                    last_output_loc = self.memory_engine.get_last_outcome_loc();
                    %closest_candidate = [];
                    smallest_diff = Inf;
                    for i = wrev(1:numel(candidates))
                        abs_diff = abs(candidates(i).loc - last_output_loc);
                        if abs_diff < smallest_diff
                            smallest_diff = abs_diff;
                            %closest_candidate = candidates(i);
                        end
                        if abs_diff > self.SEARCH_RANGE_BPM*self.NFFT/(60*self.FS)
                            search_range_valid(i) = 0;
                        end
                    end
                    
                    chosen_candidate = [];
                    if max(search_range_valid) == 0 % ie no single valid
                        chosen_candidate = candidates(1);
                    else
                        for i = 1:numel(candidates)
                            if search_range_valid(i) == 1
                                chosen_candidate = candidates(i);
                                break;
                            end
                        end
                    end
                    
                    % verify we aren't following the wrong track
                    for i = 1:numel(candidates)
                        if candidates(i).Ct > self.VERIFICATION_Ct && candidates(i).Cta > self.VERIFICATION_Cta && candidates(i).Cth > self.VERIFICATION_Cth && candidates(i).Cta > self.VERIFICATION_K*chosen_candidate.Cta
                            chosen_candidate = candidates(i);
                            %disp('override!!!')
                            %a = self.output_i
                            break;
                        end
                    end

                else
                    chosen_candidate = candidates(1);
                end
                                
                self.memory_engine.update_last_outcome(chosen_candidate);
                self.memory_engine_repeated = 0;
            end
        end
        
        
        % Debug plots
        function plot_candidates(self, candidate_frames, acc_peak_frames, output)
            figure()
            hold on
            for i = 1:numel(candidate_frames)
                for j = 1:numel(candidate_frames(i).loc)
                    plot(i, candidate_frames(i).bpm(j), '*g');
                end
            end
            
            for i = 1:numel(acc_peak_frames)
                for j = 1:numel(acc_peak_frames(i).loc)
                    plot(i, acc_peak_frames(i).bpm(j), '*y');
                end
            end
            
%             if numel(self.reference_output) > 0
%                 plot(self.reference_output, '-r');
%             end
            plot(output, '-b')
            
            hold off
        end
        
        
        % Debug plots
        function plot_confidence(self, candidate_frames)
            
            indices = [];
            bpms = [];
            Cp = [];
            Ca = [];
            Ch = [];
            Ct = [];
            
            indices1 = [];
            bpms1 = [];
            Ctp = [];
            Cta = [];
            Cth = [];
            
            for i = 1:numel(candidate_frames)
                indices = [indices, i*ones(1, numel(candidate_frames(i).loc))];
                bpms = [bpms, candidate_frames(i).bpm];
                Cp = [Cp, candidate_frames(i).Cp];
                Ca = [Ca, candidate_frames(i).Ca];
                Ch = [Ch, candidate_frames(i).Ch];
                Ct = [Ct, candidate_frames(i).Ct];
                
                if candidate_frames(i).Ctp(1) ~= -1
                    indices1 = [indices1, i*ones(1, numel(candidate_frames(i).loc))];
                    bpms1 = [bpms1, candidate_frames(i).bpm];
                    Ctp = [Ctp, candidate_frames(i).Ctp];
                    Cta = [Cta, candidate_frames(i).Cta];
                    Cth = [Cth, candidate_frames(i).Cth];
                end
            end
            
            figure()
            scatter(indices, bpms, [], Cp, '*');
            colorbar;
            title('Cp')
            xlabel('Output i');
            ylabel('Frequency (bpm)');
            
            figure()
            scatter(indices, bpms, [], Ca, '*');
            colorbar;
            title('Ca')
            xlabel('Output i');
            ylabel('Frequency (bpm)');
            
            figure()
            scatter(indices, bpms, [], Ch, '*');
            colorbar;
            title('Ch')
            xlabel('Output i');
            ylabel('Frequency (bpm)');
            
            figure()
            scatter(indices, bpms, [], Ct, '*');
            colorbar;
            title('Ct')
            xlabel('Output i');
            ylabel('Frequency (bpm)');
            
            figure()
            scatter(indices1, bpms1, [], Ctp, '*');
            colorbar;
            title('Ctp')
            xlabel('Output i');
            ylabel('Frequency (bpm)');
            
            figure()
            scatter(indices1, bpms1, [], Cta, '*');
            colorbar;
            title('Cta')
            xlabel('Output i');
            ylabel('Frequency (bpm)');
            
            figure()
            scatter(indices1, bpms1, [], Cth, '*');
            colorbar;
            title('Cth')
            xlabel('Output i');
            ylabel('Frequency (bpm)');
        end
        
    end
    
    methods (Static)
        
        function candidates = sharpen_fund_freq_candidates(candidates, power_spectrum, beta)
            for i = 1:numel(candidates)
                if candidates(i).loc*2 < numel(power_spectrum)
                    candidates(i).power = candidates(i).power + beta*power_spectrum(candidates(i).loc*2);
                end
            end
            % return candidates
        end
                
    end
    
end
