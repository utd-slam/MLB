%% MLB_language_v1
% Script to run FST for YA (as part of MLB), ported from my previous 
% isss_multi script. 
% Author - Matt Heard (heardmatthew49@gmail.com)

% MM/DD/YY -- CHANGELOG
% 08/07/17 -- Started changelog. MH
% 07/09/19 -- Cloned for new project. Lots of updates to make. MH
% 10/04/19 -- Updates begin in earnest with a new stimuli set. 
% 10/10/19 -- Debug mode added
% 01/13/20 -- Forked into YA version:
%   + Switched back to traditional ISSS
%   + Implementing fix for timing error
%   + Adjusted output column for easy analysis later
%   + Removed empty column (Real Stim Duration) because timing broke it
%   ~ Still need to counterbalance stimuli presentation!  
%   ~ New jitter which ensures it is no longer than (max allowed by
%     duration) or 1s long
% 01/15/20 -- Pulled to CCBBI scan computer
%   + New scheme for timing based on discussion w/ Xiangrui
%   + Had to add jp_scripts
% 01/27/20 -- New stimuli loading done! Changed to female voice stim. 
% 01/28/20 -- Removed jitter. 
% 01/31/20 -- Brought back to SLAM. New stimuli with SNRs 2 and -2.
%   Modified timing at end. 
% 02/27/20 -- New experiment design. 
% 12/07/20 -- Welcome back after COVID. Surprise, we've moved to UT Dallas.
%   Modified for testing at UT Dallas:
%   + Quicker experiment, 8 events per run of clear OR/SR sentence
%   + Total of 4 tuns, 2x2 factorial (20/32ch, Norm/No-norm)
% 02/23/21 -- And now we're doing something totally different. 
%   - No more SNR condition. Just Noise/OR/SR
%   ~ 24 events: 10 OR, 10 SR, 4 Noise
%   ~ Passive (post-run) task based on memory

clearvars; sca; 
DisableKeysForKbCheck([]); 
try; KbQueueStop; end %#ok<TRYNC,NOSEM>
clc; 

%% Startup
DEBUG = 1; 

try 
    PsychPortAudio('Close'); 
catch
    disp('PsychPortAudio is already closed.')
end

InitializePsychSound

codeStart = GetSecs(); 
Screen('Preference','SkipSyncTests', 1); % may need modification

if DEBUG
    AudioDevice = PsychPortAudio('GetDevices'); 
end

%% Parameters
if DEBUG
    warning('USING DEBUG DEFAULTS!')
    dlg_ans = {'TEST', '3', '4', '0'}; 
else
    prompt = {...
        'Subject number (Naming convention?):', ...
        'First run (1)', ... 
        'Last run (4)', ... 
        'RTBox connected (0/1):', ...
        }; 
    dlg_ans = inputdlg(prompt); 
end

subj.Num  = dlg_ans{1};
subj.firstRun = str2double(dlg_ans{2}); 
subj.lastRun  = str2double(dlg_ans{3}); 
ConnectedToRTBox = str2double(dlg_ans{4}); 

%% Scan paradigm
% Abbreviated Hybrid
% Will test Optimized Hybrid as well
p.TR     = 1.000; 
p.epiNum = 8;  % testing hybrid_isss version!

% Timing
p.runsMax = 4; % New design has more stimuli and longer runs
p.events  = 24; % Events per block. 
p.silent  = 4;  % silent trials
p.sentences  = 20; % how many sentence stimuli per block?
p.structures = 80; % how many sentence structures?

p.presTime = 4.000; % 4 seconds
p.jitter   = 0; 
% No jitter because FIR (Thanks Xiangrui. Sergey, any thoughts?)

p.rxnWindow = 0.000;  % 3 seconds, should we expand this?
% no reaction time window: no response

p.epiTime   = p.TR * p.epiNum;  % 
p.eventTime = p.presTime + p.epiTime;
p.runDuration = p.epiTime + ...   % After first pulse
    p.eventTime * p.events; % Each event
% run duration does not count T1 stabilization

%% Paths
cd ..
dir_exp = pwd; 

dir_stim = fullfile(dir_exp, 'stimuli');
dir_stim_all = fullfile(dir_stim, 'MLB_language'); 

dir_scripts = fullfile(dir_exp, 'scripts');
dir_results = fullfile(dir_exp, 'results');

%% Preallocating timing variables
real_eventStart = NaN(p.events, p.runsMax);
real_stimStart  = NaN(p.events, p.runsMax); 
real_stimEnd    = NaN(p.events, p.runsMax); 
real_eventEnd   = NaN(p.events, p.runsMax); 

real_respTime = cell(p.events, p.runsMax); 
real_respKey  = cell(p.events, p.runsMax); 

abs_eventStart = NaN(p.events, p.runsMax); 
abs_stimStart  = NaN(p.events, p.runsMax); 
abs_stimEnd    = NaN(p.events, p.runsMax); 
abs_rxnEnd     = NaN(p.events, p.runsMax); 
abs_eventEnd   = NaN(p.events, p.runsMax); 

abs_eventStart_delta = NaN(p.events, p.runsMax); 
abs_stimStart_delta  = NaN(p.events, p.runsMax); 
abs_stimEnd_delta    = NaN(p.events, p.runsMax); 
abs_rxnEnd_delta     = NaN(p.events, p.runsMax); 
abs_eventEnd_delta   = NaN(p.events, p.runsMax); 

t1 = NaN(p.events + 1, p.runsMax); 

firstPulse = NaN(1, p.runsMax); 
runEnd     = NaN(1, p.runsMax); 

%% File names
results_xlsx = ['MLB_' subj.Num '_language_v1_hybrid.xlsx']; 
results_mat  = ['MLB_' subj.Num '_language_v1_hybrid.mat']; 

%% Create keys
% Each block consists of 24 sentences. There will be 12 OR and 12 SR
% sentences per block. Also, four silent trials. 
% KEYS TO THE CODE
% key_events: what to play, and at what point
% key_jitter: how much jitter?
% key_eventStart: when does each event start?
% key_stimStart: when does each stimuli start?
% key_stimEnd: when does stimuli end?
% key_eventEnd: when does event end?
% key_stimuli: what stimuli was played per event?

% key_events
% event 321 is silent trial
key_events = nan(p.events, p.runsMax); 
sentences = reshape(1:4:4*p.structures, [p.sentences, p.runsMax]); 
for rr = 1:p.runsMax
    thisrun = sentences(:, rr);
    syntax = Shuffle(repmat(0:3, [1 p.sentences/4])'); 
    thisstim = Shuffle(thisrun) + syntax; 
    key_events(:, rr) = Shuffle([thisstim; repelem(321, p.silent)']); 
end

% CHECK: There should be p.runsMax * p.sentences + 1 unique items
assert(length(unique(key_events)) == p.runsMax * p.sentences + 1)

key_jitter = 0.1*ones(p.events, p.runsMax); % add slight delay
temp = p.epiTime + [0:p.eventTime:((p.events-1)*p.eventTime)]'; %#ok<NBRAK>
key_eventStart = repmat(temp, [1, p.runsMax]); 
key_stimStart  = key_eventStart + key_jitter; 
key_eventEnd   = key_eventStart + p.eventTime;
key_stimuli = cell(p.events, p.runsMax);

%% Load stimuli
% files_all = dir(fullfile(dir_stim_all, '*.wav')); 
ad_events = cell(p.events, p.runsMax); 
fs_events = zeros(p.events, p.runsMax); 
% We are loading audio data in the order of stimuli presentation!

disp('loading all stimuli...')
files = dir(fullfile(dir_stim_all, '*.wav')); fname = {files.name}';
files = fullfile(dir_stim_all, fname)';

for rr = 1:p.runsMax
    for ee = 1:p.events
        key_stimuli{ee, rr} = fname{key_events(ee, rr)}; 
        thisfile = files{key_events(ee, rr)}; 

        [tempAudio, fs_events(ee, rr)] = audioread(thisfile); 
        ad_events{ee, rr} = [tempAudio'; tempAudio']; 
    end
    
end

disp('done!')

% make sure all files have same sampling rate
assert(all(fs_events(1) == fs_events(:)))
fs = fs_events(1); 

% Clean up
dur_all = cellfun((@(x) length(x)/fs), ad_events); 
assert(~any(dur_all(:) > p.presTime))
clear tempAudio

% Few more keys that depend on stim data!
key_stimEnd    = key_stimStart  + dur_all;

%% PTB
if DEBUG
    [wPtr, rect] = Screen('OpenWindow', 1, 185, [0 0 1280 720]);
else
    [wPtr, rect] = Screen('OpenWindow', 0, 185); % may need adjustment
end

DrawFormattedText(wPtr, 'Please wait, preparing experiment...');
Screen('Flip', wPtr);
centerX = rect(3)/2;
centerY = rect(4)/2;
crossCoords = [-30, 30, 0, 0; 0, 0, -30, 30]; 
HideCursor(); 

pahandle = PsychPortAudio('Open', 5, [], [], fs); 
% 5 in lab, 3 at home what is it at CBH? Check audio devices
%%% PTB note:
% MS-Windows you’ll have the choice between up to 5 different audio 
% subsystems:
% WASAPI (on Windows-Vista and later), or WDMKS (on Windows-2000/XP) should
% provide ok latency.
% DirectSound is the next worst choice if you have hardware with 
% DirectSound support.
% If everything else fails, you’ll be left with MME, a completely unusable 
% API for precise or low latency timing.

RTBox('fake', ~ConnectedToRTBox); 
RTBox('UntilTimeout', 1);
Screen('TextSize', wPtr, 42); 

%% Prepare test
try
    for rr = subj.firstRun:subj.lastRun
        % Wait for first pulse
        DrawFormattedText(wPtr, ['Waiting for first pulse. Run ' num2str(rr)], 'center', 'center'); 
        Screen('Flip', wPtr); 
        
        RTBox('Clear'); 
        firstPulse(rr) = RTBox('WaitTR'); 

        % Draw onto screen after recieving first pulse
        Screen('DrawLines', wPtr, crossCoords, 2, 0, [centerX, centerY]);
        Screen('Flip', wPtr); 
        
        WaitTill(firstPulse(rr) + p.epiTime - 1.5*p.TR); % edits Xiangrui
        if ~DEBUG
            t1(1, rr) = RTBox('WaitTR');
        else
            t1(1, rr) = GetSecs(); 
        end
        
        %% Present audio stimuli
        for ev = 1:p.events
            PsychPortAudio('FillBuffer', pahandle, ... 
                ad_events{ev, rr});
            
            abs_stimStart_delta(ev, rr) = t1(ev, rr) + p.TR + key_jitter(ev, rr); 
            WaitTill(abs_stimStart_delta(ev, rr) - 0.1); 
            
            real_stimStart(ev, rr) = PsychPortAudio('Start', pahandle, ... 
                1, abs_stimStart_delta(ev, rr), 1);
            
            % added delay because of scanner error
            WaitTill(t1(ev, rr) + p.eventTime - 0.5*p.TR); % can wait a while!
            if ev ~= p.events
                if ~DEBUG
                    t1(ev + 1, rr) = RTBox('WaitTR'); % does not need 'Clear'--TRs don't buffer
                else
                    t1(ev + 1, rr) = GetSecs(); 
                end
                
            end
            
        end

        runEnd(rr) = GetSecs(); 

        if rr ~= subj.lastRun
            endstring = 'End of run. Let us know when you are ready to continue'; 
            DrawFormattedText(wPtr, endstring, 'center', 'center'); 
            Screen('Flip', wPtr); 
            RTBox('WaitTR'); 
            pause(2)
        end 
        
    end
    
catch err
    sca; 
    runEnd(rr) = GetSecs();  %#ok<NASGU>
    cd(dir_scripts)
    OutputData_MLB_language
    PsychPortAudio('Close'); 
    rethrow(err)
end

%% Closing down
Screen('CloseAll');
PsychPortAudio('Close'); 
DisableKeysForKbCheck([]); 

%% Save data
cd(dir_scripts)
disp('Please wait, saving data...')
OutputData_MLB_language
disp('All done!')
