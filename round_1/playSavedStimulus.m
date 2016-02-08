

clear all

% % % Set stimulus path and name % % %
stimPath = 'C:\Users\Mel\Desktop\RigConstruction\practice_2\2016_02_03\';
stimName = 'SAMc_20160203_144755'; % base name of stimulus
fileType = '.flac' ; 
stimFile = [stimPath,stimName,fileType]; % full path to stimulus

ai = audioinfo(stimFile);  % audio information
stimDur = ai.Duration;  % duration of audio file (sec)

% % % For each run, save the parameters and a time stamp % % % %
% % % % % documentation that stimulus was played % % % % % 
paramPath = 'C:\Users\Mel\Desktop\RigConstruction\practice_2\2016_02_03\' ;
paramBase = [paramPath,stimName,'_presentation_'] ;

signalScale = 2; % to scale audio signal (only needed when sound is loaded)
ttlScale = 1e3;

% % Set Queue duration manually, or let it equal full duration of stimulus
queueDur = 1 ; % (sec) duration of queue. 
%queueDur = ceil(stimDur); % queue = full stimulus

% if queueDur>1024; queueDur = 1024; 
%     disp('Queue duration was decreased to 1024 (max allowed)')
% end
scBuffer = 8192 ; % sound Card Buffer in Samples
scFs = 192e3; % sound Card Sampling Rate
scBuff = scBuffer / scFs ; % sound card buffer in seconds
disp(['Is your soundcard set to ', num2str(scFs*1e-3), 'kHz?'])    

% % Create Handles % % 
tic
q = audioinfo(stimFile) ; fs= q.SampleRate;
AFR = dsp.AudioFileReader(stimFile,'SamplesPerFrame',ceil(queueDur*fs)); % loads audio info
AP = dsp.AudioPlayer('SampleRate',AFR.SampleRate, ...
			'QueueDuration',queueDur, ...
			'OutputNumUnderrunSamples',true,...
            'BufferSizeSource','property',...
            'BufferSize',scBuffer/2);

loadtime=toc;
dateTime = datestr(now,'yymmdd_HHMMSS') ; % date and time stamp

% % % Save Parameters - will re-save over file name when stimulus is finished
fullPlay = 'no' ;
paramSave = [paramBase,dateTime,'.mat']; % full stimulus presentation save name
save(paramSave,'-regexp', '^(?!(AFR|AP)$).') % save variables except the audioplayer handle
% %
tic
disp('...Playing Audio Stimulus....')

while ~isDone(AFR)
  audio = step(AFR); % 1-step: Load audio file
  audio(:,1)=audio(:,1).*signalScale; % Scale signal
  audio(:,2)=audio(:,2).*ttlScale; %Scale ttl
  nUnderrun = step(AP,audio); % 1-step: Send audio to sound card
  if nUnderrun > 0
    fprintf('Audio player queue underrun by %d samples.\n'...
	     ,nUnderrun);
  end
end
pause(AP.QueueDuration); % wait until audio is played to the end
release(AFR);            % close the input file
release(AP);             % close the audio output device
disp('...Stimulus Completed!!!')

fH=figure; bar(1,1,'facecolor','m');xlim([0.6 1.4])
tt=title('ALL DONE PLAYING STIMULUS'); set(tt,'fontsize',20);
set(gca,'LooseInset',get(gca,'TightInset'),'ytick','','xtick','')
set(gcf,'Visible', 'on'); 
set(fH,'position',[0,0,5e3,5e3]);movegui(fH,'center'); drawnow;shg

% % SAVE % %
fullPlay = 'yes' ;
playtime = toc;
save(paramSave,'-regexp', '^(?!(AFR|AP|audio|fH)$).') % save variables except the audioplayer handle
