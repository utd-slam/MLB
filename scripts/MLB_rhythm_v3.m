%% MLB_language_v2
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
% 04/15/21 -- v2 features updates from CBH testing. Now developing for
%   Linux laptop, thank goodness. 
% 05/14/21 -- New tasks
% 05/20/21 -- New language stimuli. Sticking with vocoded stimuli as
%   baseline pending discussion with Yune. New timing again for shorter
%   runs:
%   ~ 20 events: 8 OR, 8 SR, 4 noise
%   ~ active task: typical gender-based OR/SR task
% 06/21/21 -- Merging with work from linux device. Let's just call this one
%   v3. Starting rhythm task
% 06/24/21 -- Picked final stimuli. Changed task. AGAIN. 

clearvars; sca; 
PsychDefaultSetup(2)
DisableKeysForKbCheck([]); 
try; KbQueueStop; end %#ok<TRYNC,NOSEM>
clc; 

%% Startp
DEBUG = 0; 
whichScreen = 1;  % set to 1 at CBH?
whichAudio  = 4; 

try 
    PsychPortAudio('Close'); 
catch
    disp('PsychPortAudio is already closed.')
end

InitializePsychSound

codeStart = GetSecs(); 

if DEBUG
    AudioDevice = PsychPortAudio('GetDevices'); 
    Screen('Preference', 'SkipSyncTests', 1);
else
    Screen('Preference', 'SkipSyncTests', 0); % may need modification
end

%% Parameters
if DEBUG
    warning('USING DEBUG DEFAULTS!')
    dlg_ans = {'TEST', '1', '4'}; 
else
    prompt = {...
        'Subject number:', ...
        'First run (1)', ... 
        'Last run (4)', ... 
        }; 
    dlg_ans = inputdlg(prompt); 
end

subj.Num  = dlg_ans{1};
subj.firstRun = str2double(dlg_ans{2}); 
subj.lastRun  = str2double(dlg_ans{3}); 
ConnectedToRTBox = 0; % ALWAYS OFF

%% Scan paradigm
% Abbreviated Hybrid
% Will test Optimized Hybrid as well
p.TR     = 0.500; 
p.epiNum = 16;  % testing hybrid_isss version!

% Timing
p.runsMax  = 4;  % New design has more stimuli and longer runs
p.events   = 20; % Events per block. 
p.oddball  = 2;  % oddball trials 
p.baseline = 2;  % baseline trials (constant tone)
p.rhythms  = 16; % how many rhythms per block?
p.sounds   = 64; % how many sounds total?

p.presTime = 4.000; % 4 seconds
p.jitter   = 0; 
% No jitter because FIR (Thanks Xiangrui. Sergey, any thoughts?)
% Each trial begins with 0.5s of silence, followed by sentences that are no
% more than 3s long. Each trial ends with at least 0.5s of silence before
% the scanner turns on. 

p.rxnWindow = 3.000;  % 3 seconds after stimuli finishes 
                      % should we expand this?

p.epiTime   = p.TR * p.epiNum;  
p.eventTime = p.presTime + p.epiTime;
p.runDuration = p.epiTime + ...   % After first pulse
    p.eventTime * p.events; % Each event
% run duration does not count T1 stabilization

%% Paths
cd ..
dir_exp = pwd; 

dir_stim        = fullfile(dir_exp, 'stimuli');
dir_stim_rhythm = fullfile(dir_stim, 'MLB_rhythm_48k_rms'); 
dir_scripts     = fullfile(dir_exp, 'scripts');
dir_results     = fullfile(dir_exp, 'results');

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

t1 = NaN(p.events + 1, p.runsMax); 

firstPulse = NaN(1, p.runsMax); 
runEnd     = NaN(1, p.runsMax); 

%% File names
results_xlsx = ['MLB_' subj.Num '_rhythm_v3.xlsx']; 
results_mat  = ['MLB_' subj.Num '_rhythm_v3.mat']; 

%% Create keys
% Each block consists of 16 sentences and 4 noise trials. 8 OR and 8 SR
% KEYS TO THE CODE
% key_events: what to play, and at what point
% key_onset: how much delay before sound onset?
% key_eventStart: when does each event start?
% key_stimStart: when does each stimuli start?
% key_stimEnd: when does stimuli end?
% key_eventEnd: when does event end?
% key_stimuli: what stimuli was played per event?

% key_events
% event 998 is oddball (same) trial
% event 999 is baseline (tone) trial
key_events = nan(p.events, p.runsMax); 
% rhythms = reshape(1:p.sounds, [p.rhythms, p.runsMax]); 
%  1-16 is clarinet (skipping oddball and tone... they are loaded separate)
% 17-32 is piano
% 33-48 is trombone
% 49-64 is vibraphone
% Then, let's choose tone and oddball stimuli. 
% 65 -- clarinet irregular oddball
% 66 -- clarinet regular oddball
% 67 -- clarinet constant tone
% 68 -- piano irreg odd
% 69 -- piano reg odd
% 70 -- piano tone
% 71 -- tbone ""
% 72 -- "" ""
% 73 -- "" ""
% 74 -- vibraphone ""
% 75 -- "" ""
% 76 -- "" "" 
%%% Probably easiest to manually "seed" and shuffle
%%% Let's just go diagnonally. 
rhythms = [1 2 19 20 37 38 55 56  9 10 27 28 45 46 63 64; ...
           3 4 21 22 39 40 49 50 11 12 29 30 47 48 57 58; ...
           5 6 23 24 33 34 51 52 13 14 31 32 41 42 59 60; ...
           7 8 17 18 35 36 53 54 15 16 25 26 43 44 61 62; ...
           ]'; 
%%% We seed these too
oddtones = [65 69 73 76; ... % clar irreg, piano reg, tbone tone, vibe tone
            66 70 73 74; ... % clar reg, piano tone, tbone tone, vibe irreg
            67 70 71 75; ... % clar tone, piano tone, tbone irreg, vibe reg 
            67 68 72 76; ... % clar tone, piano irreg, tbone reg, vibe tone
            ]';       
%%% Shuffle which run and oddtones are presented
runKey = Shuffle(1:4); 
odtKey = Shuffle(1:4); 
for rr = 1:p.runsMax
    thisrun =  rhythms(:, runKey(rr));
    thisodt = oddtones(:, odtKey(rr)); 
    key_events(:, rr) = Shuffle([thisrun; thisodt]); 
end

temp           = p.epiTime + [0:p.eventTime:((p.events-1)*p.eventTime)]'; %#ok<NBRAK>
key_onset      = 0.5*ones(p.events, p.runsMax); % add slight delay
key_eventStart = repmat(temp, [1, p.runsMax]); 
key_stimStart  = key_eventStart + key_onset; 
key_eventEnd   = key_eventStart + p.eventTime;
key_stimuli    = cell(p.events, p.runsMax);

%% Load stimuli
% files_all = dir(fullfile(dir_stim_all, '*.wav')); 
ad_events = cell(p.events, p.runsMax); 
fs_events = zeros(p.events, p.runsMax); 
% We are loading audio data in the order of stimuli presentation!

disp('loading all stimuli...')
files = dir(fullfile(dir_stim_rhythm, '*.wav')); 
fname = {files.name}';

%%% Separate oddballs and tones from rhythms
these_oddtones = contains(fname, 'oddball') | contains(fname, 'tone'); 
foddtones = fname(these_oddtones);
frhythms = fname; frhythms(these_oddtones) = []; 
fname_all = [frhythms; foddtones]; 
fname_all_full = fullfile(dir_stim_rhythm, fname_all);

for rr = 1:p.runsMax
    for ee = 1:p.events
        key_stimuli{ee, rr} = fname_all{key_events(ee, rr)}; 
        thisfile = fname_all_full{key_events(ee, rr)}; 
        
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
key_stimEnd = key_stimStart  + dur_all;
key_rxnEnd  = key_stimEnd + p.rxnWindow; 

%% PTB
if DEBUG
    [wPtr, rect] = Screen('OpenWindow', whichScreen, 185, [0 0 1280 720]);
else
    [wPtr, rect] = Screen('OpenWindow', whichScreen, 185); % may need adjustment
    % 022521 -- set to display 3, and adjust Windows settings
    %%% Will need updating on laptop
    % 052621 -- mirror display in settings, then screen 0
end

DrawFormattedText(wPtr, 'Please wait, preparing experiment...');
Screen('Flip', wPtr);
centerX = rect(3)/2;
centerY = rect(4)/2;
crossCoords = [-30, 30, 0, 0; 0, 0, -30, 30]; 
HideCursor(); 

pahandle = PsychPortAudio('Open', whichAudio, [], [], fs); 
% 5 in lab, 3 at home what is it at CBH? Check audio devices
% 022521 -- CBH uses device 3 (Windows WASAPI) but plays at 48k
%%% Will need updating on laptop
% 052021 -- laptop front speaker is 0, uses ALSA
% 052621 -- using headset output, so device 0 (ALSA based)

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

        abs_rxnEnd(:, rr)  = key_rxnEnd(:, rr)  + firstPulse(rr); 
        
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
            
            abs_stimStart(ev, rr) = t1(ev, rr) + p.TR + key_onset(ev, rr);
            abs_stimEnd(ev, rr)   = abs_stimStart(ev, rr) + dur_all(ev, rr); 
            abs_rxnEnd(ev, rr)    = abs_stimEnd(ev, rr) + p.rxnWindow;
            
            WaitTill(abs_stimStart(ev, rr) - 0.1); 
            
            real_stimStart(ev, rr) = PsychPortAudio('Start', pahandle, ... 
                1, abs_stimStart(ev, rr), 1);
            
            % wait until end of stim
            WaitTill(abs_stimEnd(ev, rr)); 
            RTBox('Clear')
            
            [real_respTime{ev, rr}, real_respKey{ev, rr}] = RTBox(abs_rxnEnd(ev, rr)); 
            
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
            RTBox('Clear')
            RTBox(inf); 
            pause(1)
        end 
        
    end
    
catch err
    sca; 
    runEnd(rr) = GetSecs();  %#ok<NASGU>
    cd(dir_scripts)
    OutputData_MLB_rhythm_v3
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
OutputData_MLB_rhythm_v3
disp('All done!')
