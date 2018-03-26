%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-processing:
% Bandpass, bandstop, notch filter, peak detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load data from text file folder
ecg = load('1002867.txt');

% Signal Variables
L = length(ecg);
fs = 500;
t = (1:L);  %Don't use
t1=(0:L-1)/fs;
Fv = linspace(0, 1, fix(L/2)+1)*(fs/2);
Iv = 1:length(Fv);
Lead=num2cell(ecg,1);

% Filtering variables
order = 6;
fcutlow=0.5;   %low cut frequency in Hz
fcuthigh=100;   %high cut frequency in Hz
[b,a]=butter(order,[fcutlow,fcuthigh]/(fs/2),'bandpass');

%60Hz PLI notch filter
bandstop = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',fs);
          
%60Hz filter magnitude response
fvtool(bandstop,'Fs',fs)
%60Hz filter impulse response
impz(bandstop,50)


bandpass = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',0.5,'HalfPowerFrequency2',100, ...
    'SampleRate',fs);

%Bandpass filter magnitude response
fvtool(bandpass, 'fs', fs);
%60Hz filter impulse response
impz(bandpass,50)

        
i = 2;
filtLead={};

while i < 14
    filtLead{1,i-1} = filtfilt(bandstop,Lead{1,i});
    filtLead{1,i-1} = filtfilt(bandpass,filtLead{1,i-1});
    i = i+1;
end

i = 1;
peak={};
time={};
threshold=[];

while i < 13
    peak{1,i} = [];
    time{1,i} = [];
    threshold(i) = 0.5*max(filtLead{1,i});
    
    for j = 2: L-1
        if((filtLead{1,i}(j,1) > filtLead{1,i}(j+1,1))&&(filtLead{1,i}(j,1) > filtLead{1,i}(j-1,1)) && (filtLead{1,i}(j,1) > threshold(i)))  
            peak{1,i}(j,1) = filtLead{1,i}(j,1);                                   
            time{1,i}(j,1)=(j-1)/fs;     %position stored where peak value met;              
        end
        j=j+1;
    end
  
    peak{1,i}(peak{1,i} == 0) =[];
    time{1,i}(time{1,i}==0)=[];
    
    i = i+1;
end

figure(2);
i = 1;
while i < 13
   subplot(12,1,i);
   plot(t1, filtLead{i});
   hold on;
   plot(time{i}, peak{i}, '*r');
   hold on;
   plot(t1, Lead{1,i+1}, 'LineWidth',3);
   i = i+1;
end

%Tom's
%figure(3)
%plot(Fv,abs(fffleadI(Iv))*2);
%title('Freq Notch + BP Filtered');


% Baseline drift elimination:
% Output is 'D' the detrended signal for 12 leads, after filtering
% C is the cellularized matrix of the 12 leads
filt_lead_mat = cell2mat(filtLead);
C=num2cell(filt_lead_mat,1); %call C{1}...C{12}

%Debugging: Simple plot of original ecg lead 4
F=num2cell(ecg,1);
figure
plot(t1, F{5});

%Initialize detrended lead data matrix 'D' to be zeros
D = zeros(L, 12);

j=1;
figure
while j < 13
    [p,s,mu] = polyfit((1:numel(C{j}))',C{j},6);
    f_y(:,j) = polyval(p,(1:numel(C{j}))',[],mu);
    
    D(:,j) = C{j} - f_y(:,j);   % Detrend data using polyfit curve f_y

    subplot(2,1,1);
    plot(t1, D);
    title('Detrended ECG Signal'), xlabel('Time (sec)'), ylabel('Voltage(uV)')
    axis([0 10 -500 500]);  %Mutual axes
    
    subplot(2,1,2);
    plot(t1, C{j});
    title('Original ECG Signal'), xlabel('Time (sec)'), ylabel('Voltage(uV)')
    axis([0 10 -500 500]);
    
    hold on
    
    j=j+1;
end
j=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Processing
%Misplaced electrode, low and high amplitude, motion artifacts,
%flat line detection
%Input should be D - detended/filtered data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Flags for signal discrepencies
%All flags set to 0 for initial conditions of 12 leads
F_RA_LA = 0;   %RA and LA lead reversal
F_RA_LL = 0;   %RA and LL lead reversal
F_EMG = 0;     %Max amplitude set
F_Flat = 0;    %Flat signal line
F_Min = 0;     %Minimum amplitude on 3 channels
F_Max = 0;     %Saturation on 1 channel

%Indices set
k=1;
j=1;

%Initialize matrices to zeros
max_ampl = zeros(1,12);
min_ampl = zeros(1,12);
max_sum = zeros(1,12);
max_avg = zeros(1,12);
min_sum = zeros(1,12);
min_avg = zeros(1,12);

%Lead Initialization
%filtLead = lead data after all filtering
l_1=D(:,1);    %Lead 1
l_2=D(1,2);    %Lead 2
l_3=D(1,3);    %Lead 3
l_4=D(1,4);    %aVR
l_5=D(1,5);    %avL
l_6=D(1,6);    %avF
l_7=D(1,7);    %V1
l_8=D(1,8);    %V2
l_9=D(1,9);    %V3
l_10=D(1,10);  %V4
l_11=D(1,11);  %V5
l_12=D(1,12);  %V6

% Saturation (high amplitude) condition and flag setting:
% If any lead is considered high amplitude (greater than 2mV for excerpt of continuous 200 msec ),
% F_max is set to 1, and 12-lead ECG is considered unacceptable due to motion artifacts or EMG
% noise.
checkhigh = 0;
x=0;
k=1;

while k < 13
    while x < 50
        count_1 = 1+(x*100);
        count_2 = 100+(x*100);
        %column_x = x+1;
        
        %200msec = 100 samples out of 5000 samples per D value
        %s_mat = D(1:100, 1);    %trial
        submatrix_D = D(count_1:count_2, k);
        max_submatrix_D = max(submatrix_D);
        
        if (max_submatrix_D > 2000)
            F_Max = 1;
        end
        x=x+1;
    end
    x=0;
    k=k+1;
    
end

% Flat line condition and flag setting:
% If any lead is considered to to a constant 0 voltage for at least 1sec
% length, F_Flat is set to 1, and 12-lead ECG is considered unacceptable.
checkflatline = 0;
y=0;
k=1;

while k < 13
    while y < 10
        f_count_1 = 1+(y*500);
        f_count_2 = 500+(y*500);
        %column_y = y+1;
        
        %200msec = 100 samples out of 5000 samples per D value
        %s_mat = D(1:100, 1);    %trial
        f_submatrix_D = D(f_count_1:f_count_2, k);
        f_max_submatrix_D = max(f_submatrix_D);
        
        if (f_max_submatrix_D == 0)
            F_Flat = 1;
        end
        y=y+1;
    end
    y=0;
    k=k+1;
    
end
k=1;

%
%Search for maximum of peak values in peak[] of 12 leads
while k < 13
    max_ampl(k) = max(peak{1,k});
    min_ampl(k) = min(peak{1,k});
    
    % Voltage averaging of peaks for 12 leads
    max_sum(k) = sum(peak{1,k});
    max_avg(k) = max_sum(k)/12;
    
    min_sum(k) = sum(peak{1,k});
    min_avg(k) = min_sum(k)/12;
    
    k=k+1;
end
k=1;

% Low amplitude condition and flag setting:
% If three leads are considered low amplitude (less than 125mV), F_Min is set to 1, and
% 12-lead ECG is considered unacceptable due to possible poor skin-electrode
% contact.
checklow = 0;
while k < 13
    if max_ampl(k) < 125
        checklow = checklow+1;
    end
    if checklow == 3
        F_Min = 1;
    end
    k=k+1;
end


%QRS Detection using thresholding to find peaks of interest
[~,locs_Rwave] = findpeaks(ECG_data,'MinPeakHeight',0.5,...
                                    'MinPeakDistance',200);


% Beats per minute calculation:
% Patients heartbeat, frequency of beats- in gui
t_total = max(t1);  %msec
bpm = length(peak{1,1})*6;  %number of peaks for 6*10sec sample











