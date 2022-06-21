
clear all; 
close all;
clc; 


% This is our filesystem structure
dataset_path = '../competition_data/compdata/';
output_path = '../competition_results/';

% Test Dataset IDs
ID = { 'TEST_S01_T01', 'TEST_S02_T01', 'TEST_S02_T02', 'TEST_S03_T02', ...
       'TEST_S04_T02', 'TEST_S05_T02', 'TEST_S06_T01', 'TEST_S06_T02',...
       'TEST_S07_T02', 'TEST_S08_T01'};        
         
tic()   
for idnb = 1 : 10
     
    % load signal
    sig = load(strcat(dataset_path, ID{idnb}, '.mat'), '-mat');
    sig = sig.sig;
    ppg = sig(1:2, :);
    acc = sig(3:5, :);
    
    % define test setup constants
    srate = 125;                             % 125 Hz
    window   = 8 * srate;                    % window length is 8 seconds
    step     = 2 * srate;                    % step size is 2 seconds
    
    windowNb = (length(sig)-window)/step + 1;  % total number of windows(estimates)
     
    
    % Initialise our system
    analyser = CandidateAnalyser();
    analyser.initialise(srate);
    
    % Compute BPM one analysis block at a time
    BPM = zeros(1, floor(windowNb));        
    for i =   1  :  windowNb
        curSegment = (i-1)*step+1 : (i-1)*step+window;
        BPM(i) = analyser.compute_block(ppg(:, curSegment), acc(:, curSegment));
    end
    
    % Save results
    save(strcat(output_path, 'Result', ID{idnb}(5:end), '.mat'), 'BPM');
    
    
    % clear variables
end
toc()

clear all