
% % CHECK THE DISTRIBUTION OF TRIAL DURATIONS FROM TRIGGERS % % 
% * Reads from multiple Intan files that have the same base name % 
% * Uses "presentation file"
% 2016_02_01 MJRunfeldt

clear all

% % % Stimulus Parameter File % % % 
dDir = '/Users/mel/Desktop/Rig/scrap/2016_02_03/'; % directory where data is stored
stimBase = 'SAMc_20160203_145320' ; % stimulus base name
base = [dDir,stimBase];  % path to stimulus file

prmFile = [base,'_params.mat']; % stimulus parameter file (necessary)

% % PARAMETERS % % 
interFile = 100 ; % Time between sequential rhd files: set in Intan Software
piLag = 8 ; % (sec) look for intan file initiated within this # of secs after presentation file
maxStepS = 0.05 ; % (sec) no triggers are longer than this value
jitterS = 5e-4 ; % (sec) permitted jitter in trig duration

try params = load(prmFile); % loads stimulus params
catch; disp('Parameter file corresponding to base name not found'); return
end

cd(dDir); playz=dir([base,'*presentation*']); % or look for file with "presentation" in title
if ~isempty(playz) % if there is a file match
    disp(['There are ',num2str(length(playz)),' presentation files'])
    for a = 1:length(playz)
        disp(playz(a).name)
    end
else disp(['There are no presentation files matching "',stimBase,'"']); return
end % END IFF presentation file exists


durs= cell(1,length(playz)); ques = zeros(1,length(playz));
for z = 1:length(playz)
    pres = load([dDir, playz(z).name]);
    presDate = pres.dateTime(1:6); % Date of presentation
    presTime = str2double(pres.dateTime(8:end)); % Time when presentation initiated

% % Locate rhd files that were recorded immediately after the presentation
% file % Time stamps must be EXACTLY "interFile" duration between % %
    cd(dDir); allRhd = dir([stimBase,'*',presDate,'*rhd*']); allRhd= {allRhd.name};
    rhdTimes = cell2mat((cellfun(@(x) str2double(x(29:34)),allRhd,'uniformoutput',0))) ; % time stamps from Intan
    postPres = rhdTimes - presTime;  % Intan files saved after presentation

    keeps = find(postPres == 1); % index of first intan file within "rhds" 
    keeps = find(postPres > 0 & postPres < piLag); % first intan file must be 
    % no longer than "piLag" after presentation file
    delay = postPres(keeps(1)) - 1 ; % adjust for lag/jitter from presentation to intan file
    while length(keeps) < length(postPres(keeps(1):end)) && ...
        (postPres(keeps(end)+1)-delay)/interFile == round( (postPres(keeps(end)+1)-delay)/interFile)
        keeps = [keeps keeps(end)+1]; % look for CONSECUTIVE files
    end
    rhds = allRhd(keeps); % names of all intan files cooresponding to presentation file
% % % % % % % % % FINISHED IDENTIFYING INTAN FILES % % % % % % % % % % % % % % % % % % 

    fs=readIntanParams(strcat(dDir, rhds{1})) ; % Read Intan sampling rate

    % % % Meta-Params % % % 
    maxStep = ceil(maxStepS*fs); jitter = jitterS*fs ;% convert to frames
    %time = [1:length(ttls)]./fs; % time in seconds

    % % % % Pack Trigger parameters into struct % % % % It is important that
    % 1st is recording onset, 2nd is trial onset, and 3rd is trial offset.
    % Changes to this structure necessitate changes in code, beginning at line
    % 103 % % % % %
    trigs.name{1} = 'recON'; trigs.width(1) = ceil(params.recTTL_wS*fs); 
    trigs.reps(1) = params.recTTL_N ;
    trigs.name{2} = 'trialON'; trigs.width(2) = ceil(params.trialTTL_wS*fs); 
    trigs.reps(2) = params.trialTTL_N ;
    trigs.name{3} = 'trialOFF'; trigs.width(3) = ceil(params.endTTL_wS*fs); 
    trigs.reps(3) = params.endTTL_N ;
    % %
    if min(abs(diff(combnk(trigs.width,2)'))) <= jitter; 
        disp('jitter is larger than difference in trigs (see line 26)'); return
    end

    % % % % INITIALIZE % % % %
    triggers = struct; cnt2 = 0 ;
    for i = 1:length(rhds)
        readIntanFile([dDir rhds{i}]) % LOAD Intan File. 

    % % Read time stamps from file. Check amp and board match. % %     
        if sum(t_amplifier==t_board_adc) ~= length(t_amplifier) 
            disp(['Amplifier and ADC Board time misalignment: ',rhd(i).name]); return
        else ttime = t_amplifier; clear t_amplifier t_board_adc t_aux_input t_supply_voltage
        end

        adc = board_adc_data; rawSig = adc(1,:); ttls = adc(2,:); clear adc % recording
        thresh = (max(ttls)+min(ttls))/ 2 ; % trigger threshold (metaParam)

    % % % % % MATCH RECORDED TRIGGERS TO THOSE SAVED IN PARAMETERS % % % % %
    idx = find(ttls>thresh); % all indicies of when signal goes over threshold
    tOn=[idx(1) idx(find(diff(diff(idx))<0) +1)]-1; % step off is when 2nd derivative is negative
    tOff=[idx(diff(idx)>1) idx(end)]; % step off is when diff > 1

    % % % % % IDENTIFY STEPS AND DURATIONS THAT MATCH KNOWN TRIGGERS % % % % % 
    % % % % this is the first trigger validation step: trigger durtion % % % % 
    step = struct ; cnt = 0; % initialize and count
    for a = 1:length(tOn) % for each threshold upswing
        hit = find(tOff> tOn(a)  &  tOff< (tOn(a)+maxStep)) ; % potential off step
        if ~isempty(hit); 
            dur = (tOff(hit(1))-tOn(a)); % duration in frames
            on = tOn(a) ; % onset of step

    % % Match step duration to known trigger durations (+/-jitter) % %  
            match = 0 ; b = 0;
            while match == 0 && b < length(trigs.width) % for each trigger OR until match found
                b = b + 1 ; 
                if dur >= (trigs.width(b)-jitter) && dur <= (trigs.width(b)+jitter) % IF step duration matches known "trig.dur"
                    cnt = cnt+1; % only log if there is a hit
                    step.trig(cnt) = b; step.trigJitter(cnt) = abs(dur-trigs.width(b)); 
                    step.dur(cnt) = dur; % duration in frames
                    step.on(cnt) = ttime(on) ; % onset of step IN SECONDS
                    match = 1; % exit while loop because match was found
                end % END IFF step duration matches trigger;
            end % END WHILE loop for per trigger
        end % END IFF offStep exists
    end % END (a) per step 


    % % % IDENTIFY NUMBER OF CONSECUTIVE STEPS, MATCH TO KNOWN TRIGGERS % % % 
    % % % % this is the second trigger validation step: number of steps % % % %
    c = 1 ; 
    while c <= length(step.trig) % for each step
        stepFor = min(max(trigs.reps), length(step.trig) - c) ; % window to look for steps
        sid = step.trig(c:c+stepFor); % step IDs
        stepDiff = find(cumsum(diff(sid))~=0); % when change in step ID is no longer == 0
        if ~isempty(stepDiff) && trigs.reps(step.trig(c)) == stepDiff(1) % succuessful match
            cnt2 = cnt2 + 1 ; Nreps = stepDiff(1) ;
            triggers.ID(cnt2) = step.trig(c) ; % ID of trigger (corresponding to "trigs" struct)
            triggers.name(cnt2) = trigs.name(step.trig(c)) ; % name of trigger
            triggers.on(cnt2) = step.on(c) ; % time when trigger started 
        elseif isempty(stepDiff) && sum(diff(sid))==0 && trigs.reps(step.trig(c)) == length(sid) % match, last step
            cnt2 = cnt2 + 1 ; Nreps = length(sid) ;
            triggers.ID(cnt2) = step.trig(c) ;
            triggers.name(cnt2) = trigs.name(step.trig(c)) ;
            triggers.on(cnt2) = step.on(c) ; 
        end
        c = c + Nreps  ; % move counter forward
    end % END(c) per single step

    if triggers.ID(1) ~= 1 % error message if recording trigger not the first one
        disp(['Haulted! First trigger does not indicate recording onset; rather "',triggers.name{1},'"'])
        llim = min(length(ttls),triggers.on(1)+(fs*2));
        figure;hold on;plot(ttls(1:llim),'k');plot(tOn,thresh,'r*');plot(tOff,thresh,'b*');xlim([0,llim])
        return
    end
    end % END (ii) per intan file
    % % % % % % % % FINISHED EXTRACTING TRIGGERS FROM INTAN FILE % % % % % % 


    tOns = find(triggers.ID == 2); % index of trial onset
    tOffs = find(triggers.ID== 3); % index of trial offset

    if length(tOns) == length(tOffs)+1
        disp(['Last trial onset has no offset. Last trial removed'])
        tOns(end) = [];
    elseif length(tOns) ~= length(tOffs)
        disp('Trial onsets/offsets mismatch')
    end

    trialDur = triggers.on(tOffs)-triggers.on(tOns) ;

    durs{z} = trialDur; ques(z) = pres.queueDur ; % RESULTS
end % ENd (b) per presentation file

%%
tRange = cell2mat(cellfun(@(x) range(x),durs,'UniformOutput',0));
figure;plot(ques,tRange.*1e3,'.-','markersize',30,'linewidth',3); xlabel('Queue Duration (sec)')
ylabel('Duration of single Trial (ms)'); ylim([0.04 0.06]);xlim([0.5 max(ques)+0.5])
title('Range of Trial Duration')

figure;hist(cell2mat(durs)); xx=xlabel('All Trial Durations (sec)');set(xx,'fontsize',21)