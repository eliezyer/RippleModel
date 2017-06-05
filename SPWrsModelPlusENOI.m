function [Wave,times] = SPWrsModelPlusENOI(RegLength, SNR, Fs, SPWrsTimes)

%This algorithm will create a hippocampal sharp-wave ripples (SPWRs) model
%together with Fast Ripples and spike events represented by gaussian
%curves, to use in different detector algorithm, it's possible to vary 
%algorithm main parameters in function input, as described below:
%
%   Output:
%
%   [Wave] = SPWrsModel(...) It returns only the raw signal, without
%   filtering and downsampling.
%                                                                          
%   [Wave,times] = SPWrsModel(...) Returns raw signal and a vector with
%   times, the Sharp Wave Ripples complexes' start indexes are marked with
%   1, others with 0. You'll need this one to know where are those
%   complexes and compare the accuracy of your detector algorithm.
%                                                                          
%   [Wave,times,PPWaves] = SPWrsModel(..) Returns PPWaves array, a
%   filtered and downsampled version of Waves, here is used a fir filter
%   of 400th order.
% 
%   Input:
%
%   SPWrsModel(RegLength) Set the register length with this variable, the
%   acceptable value is in minutes! Default value is 2 minutes.
%
%   SPWrsModel(RegLength, SNR) SNR determine signal-noise ratio (in dB)
%   of the SPWRs complex to basal activity (pink noise). Default value is 
%   -8 dB.
%   
%   SPWrsModel(RegLength, SNR, Fs) Fs is the frequency sample (enter this 
%   in Hz), usually is set to the frequency of your record system. Default 
%   value is 32 kHz
%   
%   SPWrsModel(RegLength, SNR, Fs, SPWrsTimes) SWPrsTimes is related to how
%   many frequency component are in SPWRs complex. Default value is 1.
%
%   SPWrsModel(RegLength, SNR, Fs, SPWrsTimes, Fd) Fd is the desired
%   frequency (in Hz), PPWaves will be downsampled to this frequency.
%   Default value is 1500 Hz. DO NOT SET THIS TO LESS THAN 800 Hz !!!
%
% This m-file is based on the follow paper:
% http://www.ncbi.nlm.nih.gov/pubmed/25570532

%% Parameters

if nargin < 1
    RegLength = 2;%Length of the "record signal", must be in minutes
end

if nargin < 2
    SNR = -8; %Signal-noise ratio, in dB
end

if nargin < 3
    Fs = 30000; %Frequency used to record, in Hertz
end

if nargin < 4
    SPWrsTimes = 1; %this variable determines how many frequency component the SPWrs complex will have.
end

times = zeros(RegLength*60*Fs,1);


%% Creating the SWPrs model

Wave = PinkNoise(RegLength*60*Fs);
Wave(2,:) = zeros(1,length(Wave));
A = 10^(SNR/20)*sqrt(2)*std(Wave(1,:)); % SPWRs Amplitude ||| I was thinking that SNR = 20log10(Asig/Anoise) => 10^(SNR/20) = Asig/Anoise => Asig = Anoise*10^(SNR/20) instead of 10^(SNR/20)*sqrt(2*std(Wave))
% SPWrs occur in a ratio of 0.3 ripples/sec
% (http://learnmem.cshlp.org/content/15/4/222.full), which means of 1
% ripple in 3 seconds
% Generate 1 event (SPWrs, FRs or Spindles) for every 1.5 s, with a
% likelihood of 33% each

for ind = 1:(Fs*1.5):RegLength*60*Fs
    
    event = rand;
    
    if (event <= 1/3)
        
        times(ind,1) = 1;
       
    elseif (event > 1/3 && event <= 2/3)
        
        times(ind,1) = 2;
        
    elseif (event > 2/3)
        
        times(ind,1) = 3;
        
    end
end

aux = find(times > 0);
for ind = 1:length(aux)
    
    if times(aux(ind)) == 1
        
        ComplexWidth = abs(round(100 + 25*randn)); %generates normal random complex width with 100 ms mean
        %and 25 of sd.
        count = 0;
        SPWrs = 0;
        
        while count < SPWrsTimes
            
            Fsharp = (200 + 25*randn); %SWPrs Frequency
            SPWrs = SPWrs + A*sin(2*pi*Fsharp*(0:1000/Fs:ComplexWidth)/1000); %1000/Fs to be in miliseconds
            fm = 1/(2*ComplexWidth/1000); %frequency of sinusoid envelope
            Envelope = sin(2*pi*fm*(0:1000/Fs:ComplexWidth)/1000);
            SPWRs = Envelope .* SPWrs;
            count = count + 1;
            
        end
        
        Wave(1, aux(ind):(aux(ind)+length(SPWRs)-1) ) = Wave(1, aux(ind) :(aux(ind)+length(SPWRs)-1) ) + SPWRs;
        Wave(2, aux(ind):(aux(ind)+length(SPWRs)-1) ) = 1;
        
    elseif times(aux(ind)) == 2 %fast ripples
        
        FrComplexWidth = abs(round(60 + 15*randn)); %generates normal random complex width with 100 ms mean
        %and 25 of sd.
        count = 0;
        Frs = 0;
        
        while count < SPWrsTimes
            
            FrSharp = (425 + 50*randn); %Fr Frequency
            Frs = Frs + 2*A*sin(2*pi*FrSharp*(0:1000/Fs:FrComplexWidth)/1000); %1000/Fs to be in miliseconds
            fm = 1/(2*FrComplexWidth/1000); %frequency of sinusoid envelope
            Envelope = sin(2*pi*fm*(0:1000/Fs:FrComplexWidth)/1000);
            FRs = Envelope .* Frs;
            count = count + 1;
            
        end
        
        Wave(1, aux(ind):(aux(ind)+length(FRs)-1) ) = Wave(1, aux(ind) :(aux(ind)+length(FRs)-1) ) + FRs;
        Wave(2, aux(ind):(aux(ind)+length(FRs)-1) ) = 2;
        
    elseif times(aux(ind)) == 3 %spikeArtifact
        
        x = 0:1/Fs:0.005;
        spindle = 5*A*gaussmf(x,[0.0003 0.0025]);
        
        Wave(1, aux(ind):(aux(ind)+length(spindle)-1) ) = Wave(1, aux(ind) :(aux(ind)+length(spindle)-1) ) + spindle;
        Wave(2, aux(ind):(aux(ind)+length(spindle)-1) ) = 3;
    end
        
end

% %% Plotting raw signal model
% figure;
% plot((1:length(Wave))/Fs,Wave(1,:))
% hold on
% plot((find(Wave(2,:) == 1))/Fs,max(Wave(1,:)),'sr','MarkerSize',3,'MarkerFaceColor','r')
% plot((find(Wave(2,:) == 2))/Fs,max(Wave(1,:)),'sk','MarkerSize',3,'MarkerFaceColor','k')
% plot((find(Wave(2,:) == 3))/Fs,max(Wave(1,:)),'sg','MarkerSize',3,'MarkerFaceColor','g')
% title(['Sharp Wave Ripples signal model with SNR of ' num2str(SNR) ' dB' ]);
% legend('Raw Signal','SPWRs Events','FRs Events','Spindles')
% xlabel('time (s)')
% hold off
% 
% %% Low-pass filtering signal below 260 hz band
% %[y140_220]=eegfilt(Wave,Fs,140,[ ] );
% %[Waves]=eegfilt(y140_220,Fs,[ ],220); %EEG filt
% 
% [b,a] = fir1(400,(260/Fs/2)); %Band pass filter butterworh IIR type, with hamming window
% % figure; freqz(b,1,2048,Fs);
% Waves2 = filter(b,a,Wave); %filtering the signal in both directions
% 
% 
% %% Downsampling to 1500 Hz after filtering!
% 
% % PPWaves = zeros((Fd*RegLength*60)-1);
% 
% for ds = 1:Fd*RegLength*60;
%     
%     PPWaves(ds) = Waves2(round(ds*Fs/Fd)); %significa dividir a frequencia por 2, 2=Fs/Fd, Fs/1000 = 32
%     
% end
% 
% %% Plotting filtered and downsampled signal model
% figure;
% plot((1:length(PPWaves))/Fd,PPWaves)
% hold on
% plot((find(times == 1)-1)/Fs,max(PPWaves),'sr','MarkerSize',3,'MarkerFaceColor','r')
% title(['Sharp Wave Ripples filtered signal model with SNR of ' num2str(SNR) ' dB' ]);
% legend('Filtered Signal','SPWRs Events')
% ylabel('time(s)')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Eliezyer Fermino de Oliveira                                  %
%                                                                         %
%   Version: 1.0.3 $Date: 2015/01/28 $                                    % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
