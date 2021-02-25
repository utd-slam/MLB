%% OutputData.m
% Saves data into a txt and xlsx file. Run as part of the main experiment 
% script. Author -- Matt H

% CHANGELOG
% 08/08/17  Started keeping changelog. --MH
% 08/08/17  Now allows for mulitple runs. --MH
% 08/10/17  Ready for subject 3. --MH
% 10/10/19  Updating for AM_FST. --MH
% 01/14/20  Forked for YA_FST (new timing scheme). -- MH
% 01/27/20  New key coding means new output data function. -- MH

%% Preallocate all results_mat
% Timing and subject response
reaction = NaN(p.events, p.runsMax); 
response = NaN(p.events, p.runsMax); 

% Xls file
% % headers = {'BLOCK', 'Syntax', 'Babble', 'Rate', 'Sentence', 'Stim', ...
% %     'Jitter key', 'Actual jitter', 'Actual Event duration', ... 
% %     'Subj response', 'RT'}; 
headers = {'BLOCK', ...
    'Jitter key', 'Actual jitter', 'Actual Event duration', ... 
    'Subj response', 'RT'}; 

%% Saving relevant timing information
% Convert to relative time, instead of system
runDur = runEnd - firstPulse; 

real_jitter   = real_stimStart - t1(1:end-1, :) - p.TR; 
real_eventDur = diff(t1); 

% Convert keys from cells to vectors
idx = 1; 
for ii = 1:p.events
    for jj = 1:p.runsMax
        if isempty(real_respKey{ii, jj})
            real_respTime{ii, jj} = NaN; 
        end
        reaction(ii, jj) = real_respTime{ii, jj} - real_stimStart(ii, jj); 
        response(ii, jj) = str2double(real_respKey{ii, jj});
        
        idx = idx + 1; 
    end
end

%% Checks if files already exists to prevent overwrite
cd(dir_results)

while exist(results_xlsx, 'file') == 2
	results_xlsx = [results_xlsx(1:end-5), '_new', results_xlsx(end-4:end)]; 
end

while exist(results_mat, 'file') == 2
	results_mat = [results_mat(1:end-4), '_new', results_mat(end-3:end)]; 
end

%% Begin printing to xlsx file
% Rather than do many sheets, let's just do one long table. 
data = cell(p.events*p.runsMax + 1, length(headers)); 
data(1,:) = headers; 

idx = 1; 
for rr = subj.firstRun:subj.lastRun
%     thisblock = key_events(strcmp(key_events(:, 1), num2str(rr)), :); 
    thisblock = repelem(rr, p.events)'; 
    
    M = horzcat(thisblock, key_jitter(:, rr), real_jitter(:, rr), ... 
        real_eventDur(:, rr), response(:, rr), reaction(:, rr) ...
        ); 
    
    M(ismissing(M)) = "NaN"; 

    for ii = 1:p.events
        for jj = 1:length(headers)
            data{idx+1, jj} = M(ii, jj); 
        end
        
        idx = idx + 1; 
    end
    
end

%% Save data
xlswrite(results_xlsx, data)

clear data ad_all_rms
save(results_mat)
