

classdef TrackEngine < handle
    properties
        N_FRAMES;
        NFFT;
        FS;
        frames;
        dp_parent;
        
        % used for edge score
        m;
        b;
        
        % constants for internal use
        PARENT_NULL = 0;
        
        
        % define some constants
        NUM_SKIPS_ALLOWED = 1;
        ZERO_SKIP = 1;
        ONE_SKIP = 2;
    end       
        
    methods
        
        function track_engine = TrackEngine(N_FRAMES, NFFT, FS)
            track_engine.N_FRAMES = N_FRAMES;
            track_engine.NFFT = NFFT;
            track_engine.FS = FS;
            
            % if diff is smaller than 3bpm, give it full score. if > 10 give 0. straight line in between
            loc1 = 3*NFFT/(60*FS);
            loc2 = 10*NFFT/(60*FS);
            track_engine.m = -1/(loc2-loc1);
            track_engine.b = 1 - track_engine.m*loc1;
        end
        
        
        function candidates = compute(self, candidates)
            
            % construct a frame
            frame.locs = [candidates(:).loc];
            frame.Cp = [candidates(:).Cp];
            frame.Ca = [candidates(:).Ca];
            frame.Ch = [candidates(:).Ch];
            
        
            % fill in the frames array
            if numel(self.frames) >= self.N_FRAMES
                self.frames = [self.frames(2:end), frame];
            elseif numel(self.frames) == 0
                self.frames(1).locs = frame.locs;
                self.frames(1).Cp = frame.Cp;
                self.frames(1).Ca = frame.Ca;
                self.frames(1).Ch = frame.Ch;
            else
                self.frames = [self.frames, frame];
            end
            
            scores = self.compute_iter(self.frames);
            for i = 1:numel(candidates)
                candidates(i).Ct = scores(i);
                [candidates(i).Ctp, candidates(i).Cta, candidates(i).Cth] = self.retrace(i);
            end
            
            % return candidates
        end
        
        
        function scores = compute_iter(self, frames)
                           
            dp_mem = cell(self.NUM_SKIPS_ALLOWED+1, numel(frames));

            % handle base case of frames(1) and frames(2)
            dp_mem{self.ZERO_SKIP}{1} = ones(1, numel(frames(1).locs));
            dp_mem{self.ONE_SKIP}{1} = ones(1, numel(frames(1).locs));
            if numel(frames) >= 2
                dp_mem{self.ONE_SKIP}{2} = ones(1, numel(frames(2).locs));
            end
            
            % initialise base case parent
            self.dp_parent = cell(self.NUM_SKIPS_ALLOWED+1, numel(frames));
            self.dp_parent{self.ZERO_SKIP}{1} = self.PARENT_NULL * ones(1, numel(frames(1).locs));
            self.dp_parent{self.ONE_SKIP}{1} = self.PARENT_NULL * ones(1, numel(frames(1).locs));
            if numel(frames) >= 2
                self.dp_parent{self.ONE_SKIP}{2} = self.PARENT_NULL * ones(1, numel(frames(2).locs));
            end
            
            % iterate through each following frame. each frame corresponds to an output_i
            for frame_i = 2:numel(frames)
                
                curr_frame = frames(frame_i);
                prev_frame = frames(frame_i-1);
                prev_frame_i = frame_i - 1;
                if frame_i >= 3
                    prev_two_frame = frames(frame_i-2);
                end
                prev_two_frame_i = frame_i - 2;
                
                for i = 1:numel(curr_frame.locs)
                    curr_loc = curr_frame.locs(i);
                    
                    % compute no skip score
                    score = 0;
                    max_j = 0;
                    for j = 1:numel(prev_frame.locs)
                        prev_loc = prev_frame.locs(j);
                        edge_score = self.compute_edge_score(curr_loc, prev_loc);
                        temp = edge_score*dp_mem{self.ZERO_SKIP}{prev_frame_i}(j);
                        if temp > score
                            score = temp;
                            max_j = j;
                        end
                    end
                    dp_mem{self.ZERO_SKIP}{frame_i}(i) = score;
                    self.dp_parent{self.ZERO_SKIP}{frame_i}(i) = max_j;
                    
                    % compute one skip score
                    if frame_i >= 3
                        score = 0;
                        max_j = 0;
                        for j = 1:numel(prev_two_frame.locs)
                            prev_loc = prev_two_frame.locs(j);
                            edge_score = self.compute_edge_score(curr_loc, prev_loc);
                            temp = edge_score*dp_mem{self.ZERO_SKIP}{prev_two_frame_i}(j);
                            if temp > score
                                score = temp;
                                max_j = -j; % note the hack where negative means parent is 2 frames back
                            end
                        end
                        
                        for j = 1:numel(prev_frame.locs)
                            prev_loc = prev_frame.locs(j);
                            edge_score = self.compute_edge_score(curr_loc, prev_loc);
                            temp = edge_score*dp_mem{self.ONE_SKIP}{prev_frame_i}(j);
                            if temp > score
                                score = temp;
                                max_j = j;
                            end
                        end
                        dp_mem{self.ONE_SKIP}{frame_i}(i) = score;
                        self.dp_parent{self.ONE_SKIP}{frame_i}(i) = max_j;
                    end
                    
                end
            end
            
            
            scores = dp_mem{self.ONE_SKIP}{numel(frames)};
            % return scores
        end
        
        function score = compute_edge_score(self, loc1, loc2)
            loc_diff = abs(loc1 - loc2);
            score = self.m*loc_diff + self.b;
            score = min(max(score , 0), 1);            
        end
        
        
        function [Ctp, Cta, Cth] = retrace(self, index)
            
            if numel(self.frames) == self.N_FRAMES
            
                trace = self.PARENT_NULL * ones(1, numel(self.frames)-1);
                trace_loc = zeros(1, numel(self.frames)-1);
                trace_Cp = zeros(1, numel(self.frames)-1);
                trace_Ca = zeros(1, numel(self.frames)-1);
                trace_Ch = zeros(1, numel(self.frames)-1);

                parent = index;
                trace_i = numel(self.frames)-1;

                % keep going along the self.ONE_SKIP parent trace until find the skip, then go along self.ZERO_SKIP parent trace
                for frame_i = wrev(1:numel(self.frames))
                    if parent < 0 || parent == self.PARENT_NULL
                        break;
                    end
                    
                    % fill in trace info
                    trace(trace_i) = parent;
                    if parent > 0
                        trace_loc(trace_i) = self.frames(frame_i).locs(parent);
                        trace_Cp(trace_i) = self.frames(frame_i).Cp(parent);
                        trace_Ca(trace_i) = self.frames(frame_i).Ca(parent);
                        trace_Ch(trace_i) = self.frames(frame_i).Ch(parent);
                    else
                        trace_loc(trace_i) = self.frames(frame_i-1).locs(abs(parent));
                        trace_Cp(trace_i) = self.frames(frame_i-1).Cp(abs(parent));
                        trace_Ca(trace_i) = self.frames(frame_i-1).Ca(abs(parent));
                        trace_Ch(trace_i) = self.frames(frame_i-1).Ch(abs(parent));
                    end
                    trace_i = trace_i - 1;
                    
                    % update parent by reading dp_parent
                    parent = self.dp_parent{self.ONE_SKIP}{frame_i}(parent);

                end

                % decrease frame_i because we've hit the skip. If we got here because of hitting a dead track head, it's fine also
                frame_i = frame_i - 1;
                for frame_i = wrev(1:frame_i)
                    if parent == self.PARENT_NULL
                        break;
                    end

                    % fill in trace info
                    trace(trace_i) = parent;
                    trace_loc(trace_i) = self.frames(frame_i).locs(abs(parent));
                    trace_Cp(trace_i) = self.frames(frame_i).Cp(abs(parent));
                    trace_Ca(trace_i) = self.frames(frame_i).Ca(abs(parent));
                    trace_Ch(trace_i) = self.frames(frame_i).Ch(abs(parent));
                    trace_i = trace_i - 1;
                    
                    % update parent by reading dp_parent
                    parent = self.dp_parent{self.ZERO_SKIP}{frame_i}(abs(parent));
                    
                end

                % debug print
                %trace
                %trace_loc
                %trace_Cp
                %trace_Ca
                %trace_Ch 

                % remove first element, as it will always be null zero
                Ctp = mean(trace_Cp);
                Cta = mean(trace_Ca);
                Cth = mean(trace_Ch);
   
            else
                % scores not ready
                Ctp = -1;
                Cta = -1;
                Cth = -1;
            end
        end
        
    end
end