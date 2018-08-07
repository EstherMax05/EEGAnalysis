%% Data 2 of task-related EEG modulation- Dr. Sunderam
fid = fopen('OpenBCI-RAW-task 2 TN.txt','r');
data1 = textscan(fid,'%d,%f,%f,%f,%f,%f,%f,%f,%f %*[^\n]','HeaderLines',6);
fclose(fid);

Fs = 250;   % Data sample rate (samples/sec, or Hz)
data1 = double(cat(2,data1{:})); % Concatenate cells into an array

%% Example of feature extraction
% Vertical bipolar EOG signal
EOG = data1(:,4);    % EOG for eye movement detection
EOG = detrend(EOG,'constant');   % Remove DC offset
EOGs = smooth(EOG,1*Fs,'lowess');   % Estimate drift
EOGf = EOG-EOGs;     % Remove drift
%EOGf(1,:) = [];
[n,m] = size(EOGf);
t = linspace(0,n/Fs,n); t = t(:);   % Create time vector
figure
subplot(211)
plot(t,EOGf, 'k');
xlim([t(4*Fs) t(end)]);
hold on
EOGfUpperThreshold = 90;
EOGfLowerThreshold = -100;
binaryStateVar = zeros(length(EOGf),1);
count = 0;

%Based on graph EOG threshold can be -80 to -100
for i = 1:length(EOGf)
    if EOGf(i) > EOGfLowerThreshold && EOGf(i) < EOGfUpperThreshold
    binaryStateVar(i) = 0;
    else
        binaryStateVar(i) = 100;
        count = count+1;
        % for every 10 secs if we have >=4 thick ones eyes open else eyes closed
    end
end
binaryStateVa =[];
binaryStateVaStorage =[];
lowbound = 0;
count =1;
filename = 'StairsFunction.xlsx';
sheet = 1;
binaryStateVar = xlsread(filename,sheet, 'A1:A50202');
for i = 1:200
    if mean(binaryStateVar(lowbound+1:lowbound+250))>50
    binaryStateVa = 100;
    binaryStateVaStorage =[binaryStateVaStorage, binaryStateVa];
    else
    binaryStateVa = 0;
    binaryStateVaStorage =[binaryStateVaStorage, binaryStateVa];
    end
    lowbound = lowbound+250;
    count = count+1;
end
stairs(t,binaryStateVar, 'y');
xlim([t(4*Fs) t(end)]);
xlabel('time');
legend('EOG','Binary_Function');
hold off
% Take the difference of the two EEG signals to remove ocular artifact
%% Example of feature extraction
EEG = data1(:,2:3); % Only two channels were recorded
EEG = detrend(EEG,'constant');   % Remove DC offset
[b1,a1] = butter(3,[1 55]/(Fs/2));    % bandpass filter, 1-55 Hz
EEGf = filtfilt(b1,a1,EEG); % Filter data to remove drift and noise
offset = std(EEGf(:,1)); %why?
subplot(212); plot(t,EEGf(:,1),'b',t,EEGf(:,2)+0.1*offset,'r');
axis tight
xlabel('time (s)');
title('EEG and EOG signals (uV) during repeated eye close/open task.');
% Take the difference of the two EEG signals to remove ocular artifact
EEGf = EEGf(:,1)-EEGf(:,2);
t = linspace(0,n/Fs,n); t = t(:);   % Create time vector
% Compute EEG mean-squared power in a moving 1-second window
x1 = smooth(EEGf(:,1).^2,Fs,'moving');
figure
%% 2nd Feature
% Compute EEG FFT in a moving 1-second window
lowBound = 1;
x2= [];
x2_store = [];
while lowBound<(length(x1)-250)
    x2 = fft(x1(lowBound:lowBound+250));
    lowBound = lowBound+250;
    x2_store = [x2_store, x2];
end
N = length(x2_store);
f = (0:N-1)*(Fs/N); 

% for all x2 store
% 8-13hz -> x_2store(9:15, i)
%if amp's less than 100 eyes open if > eyes closed

alpha_store = [];
 for i = 1:200
     alpha = max(abs(x2_store(9:15, i)));
     alpha_store = [alpha_store, alpha];
 end
 stairs(t(1:end-4*Fs+1),4.*binaryStateVar(4*Fs:end), 'y');
 hold on
%xlim([t(6*Fs) t(end)]);
xlabel('time');
legend('EOG','Binary_Function');
 plot(alpha_store(1, 4:end))
 hold off
 alpha_storeStair = zeros(200,1);
 
 threshold = 0:50:900;
 
 DetectedPositives = zeros(1,length(threshold));%eyesclosed
TruePositives = zeros(1,length(threshold));
DetectedNegatives = zeros(1,length(threshold));
TrueNegative = zeros(1,length(threshold));
FalsePositives = zeros(1,length(threshold));
FalseNegatives = zeros(1,length(threshold));
%FalseNegatives = 0;
specificity = zeros(1,length(threshold));
Sensitivity = zeros(1,length(threshold));
for j = 1:length(threshold)
     for i = 1:200
         if alpha_store(i)> threshold(j)
            alpha_storeStair(i) = 100;
         end
     end
DetectedPositives(j) = sum(alpha_storeStair > 0);%eyesclosed
TruePositives(j) = sum(binaryStateVaStorage > 0);
DetectedNegatives(j) = sum(alpha_storeStair == 0);
TrueNegative(j) = sum(binaryStateVaStorage == 0);
FalsePositives(j) = DetectedPositives(j)-TruePositives(j);
if FalsePositives(j)<0
    FalsePositives(j) = 0;
end
FalseNegatives(j) = DetectedNegatives(j)- TrueNegative(j);
if FalseNegatives(j)<0
    FalseNegatives(j) = 0;
end
%FalseNegatives = 0;
Sensitivity(j) = TruePositives(j)./(TruePositives(j)+FalseNegatives(j));
specificity(j) = TrueNegative(j)./(TrueNegative(j)+FalsePositives(j));
 alpha_storeStair = zeros(200,1);
end
figure 
stairs(alpha_storeStair(1:end, 1),'k')
hold on
stairs(binaryStateVaStorage(1:end),'b')
%stairs(t((4*250):end),binaryStateVar((4*250):end), 'y');
xlabel('time');
hold off
figure
plot(1-specificity,Sensitivity, '-*')
xlabel('1-specificity');
%axis tight, xlabel('time (s)'); legend('Raw EEG','Feature');
% Choose a threshold and detect when the feature crosses it
%% Accounting for location
DetectedPositives = zeros(1,length(threshold));%eyesclosed
TruePositives = zeros(1,length(threshold));
DetectedNegatives = zeros(1,length(threshold));
TrueNegative = zeros(1,length(threshold));
FalsePositives = zeros(1,length(threshold));
FalseNegatives = zeros(1,length(threshold));
%FalseNegatives = 0;
specificity = zeros(1,length(threshold));
Sensitivity = zeros(1,length(threshold));
Loc = zeros(1,length(binaryStateVaStorage));
Loc_sum = zeros(1,length(threshold));
for j = 1:length(threshold)
     for i = 1:200
         if alpha_store(i)> threshold(j)
            alpha_storeStair(i) = 100;
         end
     end
     for i = 1:200
         if alpha_storeStair(i) == binaryStateVaStorage(i)
            Loc(i) = 100;
         end
     end
     Loc_sum(j) = sum(Loc==100);
     Loc = zeros(1,length(threshold));
end



