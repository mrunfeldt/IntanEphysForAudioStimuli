%% % % CREATE AND SAVE SAMchord Stimulus and Triggers % % % 
% # * ! % # * ! % # * ! % # * ! % # * ! % # * ! % # * ! % # * ! % # * ! % #
% STEP 1: Set parameters, generate carrier frequencies and trial blocks,
% save parameters. ! % # * ! % # * ! % # * ! % # * ! % # * ! % ! % # * ! % 
% # * ! % # * ! % # * ! % # * ! % # * ! % # * ! % # * ! % # * ! % # * ! % #

clear all

savePath = 'C:\Users\Mel\Desktop\RigConstruction\practice_2\2016_02_03\' ; % path to save file
nameBase = 'SAMc_' ; % base name for saved files
saveType = '.flac' ; % type of file. Supported:   .wav  .flac  .mat  

% % STIMULUS PARAMETERS % % 
maxFs = 192e3; % Max sampling rate
% % % % NOTE: sampling rate is set as 2.* the maximum carrier frequency to
% reduce size of signal and trigger vectors % % % % 
fCenter = 100; % NOTE (Hz)
fmVec = [2.^[1:6], 5, 11]; % modulation frequencies
NperBlock = 2 ; % number of  repeats per block
trialDurS = 1; % (sec) duration of a single trial in seconds
interChordS = 0.5; % time delay between chords (sec)
tBuffS = 1 ; % (sec) delay for initial onset and offset of block
preTrigS = 0.5 ; % add prior to recording trigger (necessary)
Noct = 9 ; % Number of octaves
nPerOct = 2; % N tones per octave
m = 0.99; % modulation depth (percent)
A = 1 ; % modulation amplidute
iPhase = -pi/2; % initial phase offset
rampSec = 5e-3 ; % onset/offsetramp duration in sec

% % TRIGGER PARAMETERS % %
ttlAmp = 1 ; % amp of trigger - it's best for compression to leave at 1 and scale when executing signal
recTTL_wS = 5e-3; % width for recording trigger (sec) - also 1/rate
recTTL_N = 2 ; % number of triggers that mark recording onset
trialTTL_wS = 2e-3;trialTTL_N = 3; % trigs to demarcate beginning of each trial (chord)
endTTL_wS = 3e-3; endTTL_N = 1 ; % number to 

% % METAPARAMS (functions of parameters) % % 
N = (Noct * 2) -1; % Number of notes 
tTrials = length(fmVec) * NperBlock ; % total number of trials

% % Generate Chord Carrier frequencies % % %
[tones,fMax] = generateChordTones(fCenter,nPerOct,N,'noPlot'); %yesPlot noPlot
fMax = fMax + (max(fmVec)*(length(tones)+1)); % account FM into max frequency
fs = round( max(8000,min(fMax.*2,maxFs))  ); disp(['Sampling Rate = ',num2str(fs)])

% % Fs-dependent METAPARAMS (functions of parameters) % % 
frame = 1/fs ; preTrig = ceil(preTrigS*fs);  % convert to frames
recTTL_w = ceil(recTTL_wS*fs); trialTTL_w = ceil(trialTTL_wS*fs); % convert to frames
endTTL_w = ceil(endTTL_wS*fs); %interTime = ceil(interTimeS*fs); % convert to frames
% % Generate TTL/Trigger Waveforms % % 
recTTL= repmat([ones(1,recTTL_w).*ttlAmp, zeros(1,recTTL_w)],1,recTTL_N);
trialTTL = repmat([ones(1,trialTTL_w).*ttlAmp, zeros(1,trialTTL_w)],1,trialTTL_N);
endTTL = repmat([ones(1,endTTL_w).*ttlAmp, zeros(1,endTTL_w)],1,endTTL_N);

% % Generate PseudoRandom Trial Blocks % % 
allTs = repmat(fmVec,NperBlock); allTs = allTs(:)' ; % vector of all Fms
newIdx=randperm(length(allTs)); % pseudorandomize
fMs = allTs(newIdx) ; % vector of all Fms in PseudoRandom chronology
% %
% Generate file name, include date and time stamp and save "params" file
q = datestr(now,'yyyymmdd_HHMMSS') ; % date and time stamp for save name
saveName = [savePath,nameBase,q]; % base for file save name
save([saveName,'_params']) % save parameters

% # * ! % # * ! % # * ! % # * ! % # * ! % # * ! % # * ! % # * ! % # * ! % #
% STEP 2: Generate one entire block of SAM chords and save ! % # * ! % # * 
% # * ! % # * ! % # * ! % # * ! % # * ! % # * ! % # * ! % # * ! % # * ! % #

% % % Determine length of each trial and total block % % %
trT=[0:frame:trialDurS]; % time of single trial
iTrI = [0:frame:interChordS] ; % inter trial time
buffer = [0:frame:tBuffS]; % Block ON Buffer
pre = [0:frame:preTrigS]; % add prior to signal and trigger
durTotal = (length(fMs) * (length(trT) + length(iTrI))) + ...
    (length(buffer)*2) + length(pre);% total # of frames
disp(['Total duration = ',num2str(durTotal/fs),' secs'])

signal = zeros(1,durTotal); % initialize signal vector
trig = signal; % initialize trigger vector

trig(length(pre):length(pre)+length(recTTL)-1) = recTTL; % trigger recording onset
tCnt = length(buffer)+length(pre); % VIP; marks onset of each trial (tracks time within loop)

for j = 1:length(fMs)
    fm = fMs(j) ; % modulation frequency
% % % % CREATE SAM % % % %
    note = zeros(length(tones), length(trT)) ;
    for i = 1:length(tones) % generate chord
        fc = tones(i); % carrier frequency
        note(i,:) = A.*( 1 + (m.*sin( (2*pi.*fm.*trT) + iPhase) ) ).*sin(2*pi*fc.*trT) ;
    end
    chord = sum(note) ./ length(tones); % chord is the mean of notes (sum / # of notes)   

% % % % % APPLY (CONVOLVE) ONSET AND OFFSET RAMPS TO CHORD % % 
    rN = round(rampSec* fs) ; % # of samples for cosine ramp
    rampOn  = sin(0.5*pi*(0:rN)/rN).^2; % square sine ONSET ramp of N samples
    rampOff = fliplr(rampOn); % cosine OFFSET ramp
    head = 1:(rN+1);tail = length(chord)-rN:length(chord); % indicies of ramps
    chord(head) = rampOn.*chord(head) ; chord(tail) = rampOff.*chord(tail) ; % convolve
   
    signal(tCnt:tCnt+length(chord)-1) = chord; 
    trig(tCnt:tCnt+length(trialTTL)-1) = trialTTL; % onset triggers
    trig(tCnt+length(trT):tCnt+length(trT)+length(endTTL)-1) = endTTL; % end of trial triggers
    
    tCnt = tCnt + length(chord) + length(iTrI); % move time counter forward
end


% % % SAVE SIGNAL AND TRIGGER CHANNELS % % 
% VIP: This is where you chose your compression type: I strongly recommend
% that you read MATLAB's documentation on "audiowrite" and pay note to the
% relationship between the type of data (e.g. int32(output)), the file type
% (e.g. .wav), and the 'BitsPerSample" parameter
output(:,1) = signal ./ max(abs(signal));
output(:,2) = trig  ./ max(abs(trig));

if strcmp(saveType,'.wav')
    audiowrite([saveName,saveType],(output),fs,'BitsPerSample',32); 
elseif strcmp(saveType,'.flac')
    audiowrite([saveName,saveType],output,fs,'BitsPerSample',24);
elseif strcmp(saveType,'.mat')
    save(saveName,'output','fs')
else disp('File Format not recognized (line 124)');
end

disp(['Stimulus file name:  ',nameBase,q])

% % # * ! % # * ! % # * ! % # * ! % # * ! % # * ! % # * ! % # * ! % # * ! % 
% STEP 3 (optional): Plot a portion of the signal and power spectra % # * !
% (you don't want to view the entire signal because MATLAB will crash) * ! % 
% # * ! % # * ! % # * ! % # * ! % # * ! % # * ! % # * ! % # * ! % # * ! % #
excise = [1:round(15*fs)]; % segment
data=signal(excise); time = [1:length(excise)].*frame ;

pt = size(data,2); pt = 2^(floor(log2(pt))); range = pt/2;
f = fs*(0:range-1)/pt; Y = fft(data,pt);psd = Y.*conj(Y)/pt;
% % % Restrict range for power spec plot % % %
posPwr = find(psd > mean(psd)) ; 
if ~isempty('posPwr'); 
    posPwr(posPwr>length(f))=[];pIdx = posPwr(1):posPwr(end);
else pIdx = 1:length(f);
end
% % Plot % % 
figure;
subplot(2,1,1);hold on; plot(time,mat2gray(trig(excise)+0.5),'k','linewidth',2);
plot(time,data,'b','linewidth',1.5);xx=xlabel('Time (sec)'); yy=ylabel('Pressure');
leg=legend('Triggers','Signal');set(leg,'fontsize',12,'location','southwest')
%tt=title(['Fm = ',num2str(fm),'Hz. Fcenter = ',num2str(fCenter),'Hz. ',num2str(N),' Chords']);
set(gca,'ytick',[min(data),0,max(data)],'fontsize',14);%set(tt,'fontsize',16)
set(xx,'fontsize',16);set(yy,'fontsize',16); axis tight
subplot(2,1,2);plot(f(pIdx),psd(pIdx),'b','linewidth',2);
xx=xlabel('Frequency(Hz)'); yy=ylabel('Pwr(mV^2)');set(gca,'fontsize',14)
set(xx,'fontsize',16);set(yy,'fontsize',16); axis tight
