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
% 05/20/21  Errors coming from Linux, R2021, or ..? had to switch to table
% 06/21/21  Merging with work from Linux laptop
% 06/22/21  Forked into rhythm

%% Preallocate all results_mat
% Timing and stim info
% Xls file
headers = {'BLOCK', 'Stim', 'Actual Event duration', 'Complexity'}; 

%% Saving relevant timing information
% Convert to relative time, instead of system
runDur = runEnd - firstPulse; 

real_jitter   = real_stimStart - t1(1:end-1, :) - p.TR; 
real_eventDur = diff(t1); 

%% Checks if files already exists to prevent overwrite
cd(dir_results)

while exist(results_xlsx, 'file') == 2
	results_xlsx = [results_xlsx(1:end-5), '_new', results_xlsx(end-4:end)]; 
end

while exist(results_mat, 'file') == 2
	results_mat = [results_mat(1:end-4), '_new', results_mat(end-3:end)]; 
end

%% format responses
reaction = nan(p.events, p.runsMax);
response = nan(p.events, p.runsMax);

for ii = 1:p.events
    for jj = subj.firstRun:subj.lastRun
        if isempty(real_respKey{ii, jj})
            real_respKey{ii, jj} = '0'; 
            real_respTime{ii, jj} = NaN; 
        end
        
        if length(real_respKey{ii, jj}) > 1
            temp = str2double(real_respKey{ii, jj});
            response(ii, jj) = 10*temp(1) + temp(2); % eh, it works. 
            reaction(ii, jj) = real_respTime{ii, jj}(1) - abs_stimEnd(ii, jj); 
        else
            response(ii, jj) = str2double(real_respKey{ii, jj});
            reaction(ii, jj) = real_respTime{ii, jj} - abs_stimEnd(ii, jj); 
        end
        
    end
    
end

%% Begin printing to xlsx file
% Rather than do many sheets, let's just do one long table. 
BLOCK      = []; 
Complexity = []; 
Rhythm     = []; 
Response   = []; 
RespTime   = []; 

for rr = subj.firstRun:subj.lastRun
    BLOCK = vertcat(BLOCK, repelem(rr, p.events)'); 
    Rhythm = vertcat(Rhythm, key_stimuli(:, rr)); 
    
    complex = cell(p.events, 1); 
    temp = cellfun((@(x) contains(x, '_i')), key_stimuli(:, rr)); % irregular/complex
    complex(temp) = {'Complex'};
    temp = cellfun((@(x) contains(x, '_r')), key_stimuli(:, rr)); % regular/simple
    complex(temp) = {'Simple'};
    temp = cellfun((@(x) contains(x, 'oddball')), key_stimuli(:, rr)); % oddball
    complex(temp) = {'Oddball'};
    temp = cellfun((@(x) contains(x, 'tone')), key_stimuli(:, rr)); % tone
    complex(temp) = {'Tone'};
    Complexity = vertcat(Complexity, complex); 
    
    Response = vertcat(Response, response(:, rr)); 
    RespTime = vertcat(RespTime, reaction(:, rr)); 
end

T = table(BLOCK, Complexity, Rhythm, Response, RespTime); 

%% Save data
writetable(T, results_xlsx)

clear data ad_all_rms files ad_events fname
save(results_mat)
