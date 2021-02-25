%% check_cb
% helper function to independently verify stimuli are counterbalanced
% TEST key_stim is balanced across each dimension
tag = {'M_', 'F_', '_O', '_S', '4Hz.wav', p.vocode, '01ch', 'silence', ...
       '_OM_4Hz.wav', '_OF_4Hz.wav', '_SM_4Hz.wav', '_SF_4Hz.wav', ...
       ['_OM_4Hz_' p.vocode], ['_OF_4Hz_' p.vocode], ['_SM_4Hz_' p.vocode], ['_SF_4Hz_' p.vocode]};
num = [8 8 8 8 8 8 4 4 ...
       2 2 2 2 ...
       2 2 2 2];
mask = cell(1, length(tag)); 

for tt = 1:length(tag)
    mask{tt} = cellfun((@(x) contains(x, tag{tt})), stim_all); 
    if any(tt == [7 8])
        for ii = 1:6
            mask{ii}(mask{ii} & mask{tt}) = 0; 
        end
        
    end
    
end

key = cellfun(@find, mask, 'UniformOutput', false); 

for rr = 1:size(key_stim, 2)
    thisrun = key_stim(:, 1); 
    
    for tt = 1:length(mask)
        thiskey = key{tt}; 
        if length(find(ismember(thisrun, thiskey))) ~= num(tt)
            error('check cb')
        end
    end
    
end

disp('cb is good!')