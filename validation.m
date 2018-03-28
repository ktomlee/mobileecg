acceptable = load('RECORDS-acceptable.txt');
unacceptable = load('RECORDS-unacceptable.txt');

TotalECG = {};
results = [];
pass = 0;

l = 1;


while l < 226
    
        %TotalECG{1,l} = load([num2str(acceptable(l,1)), '.txt']);
        %if i < 226
        %    TotalECG{2,i} = load([num2str(unacceptable(i,1)), '.txt']);
        %end
    fprintf('Checking record: %i\n', l);

    ecg = load([num2str(unacceptable(l,1)), '.txt']);
    Lead=num2cell(ecg,1);

    % Signal Variables
    L = length(ecg);
    fs = 500;
    t = (1:L);  %Don't use
    t1=(0:L-1)/fs;
    Fv = linspace(0, 1, fix(L/2)+1)*(fs/2);
    Iv = 1:length(Fv);
    Lead=num2cell(ecg,1);

%Flags for signal discrepencies
%All flags set to 0 for initial conditions of 12 leads
F_RA_LA = 0;    %RA and LA lead reversal
F_RA_LL = 0;    %RA and LL lead reversal
F_EMG = 0;      %Max amplitude set
F_Flat = 0;     %Flat signal line
F_Min = 0;      %Minimum amplitude on 3 channels
F_Max = 0;      %Saturation on 1 channel
F_Baseline = 0; %Baseline drift too large on any lead
F_Contact = 0;

F=num2cell(ecg,1);

% Filtering variables
order = 6;
fcutlow=0.5;   %low cut frequency in Hz
fcuthigh=100;   %high cut frequency in Hz
[b,a]=butter(order,[fcutlow,fcuthigh]/(fs/2),'bandpass');

%60Hz PLI notch filter
bandstop = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',fs);

%1.0 - 100Hz bandpass filter
bandpass = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',1.5,'HalfPowerFrequency2',100, ...
    'SampleRate',fs);

i = 2;
filtLead={};
while i < 14
    filtLead{1,i-1} = filtfilt(bandstop,Lead{1,i});
    filtLead{1,i-1} = filtfilt(bandpass,filtLead{1,i-1});
    i = i+1;
end

% Peak finder V2: QRS Detection using thresholding to find peaks of interest&
% Daubechies wavelet transfor of ecg signal beore implementation of findpeaks
% function; 30 = min peak height of R wave; 300 = less than minimum distance 
% between subsequent peaks
% Parameter initialization
db2 = {};
db2rec = {};
ydb2 = {};
peaks = {};
locs = {};
i = 1;

while (i < 13) && (j < 15)
    db2{1,i} = modwt(ecg(:,i+1),'db2',5);
    db2rec{1,i} = zeros(size(db2{1,i}));
    db2rec{1,i}(4:5,:) = db2{1,i}(4:5,:);
    ydb2{1,i} = imodwt(db2rec{1,i},'db2');

    %peak detection
    [peaks{1,i},locs{1,i}] = findpeaks(ydb2{1,i}(1,:),t,'MinPeakHeight',15,'MinPeakDistance',150);
    [speaks{1,i},slocs{1,i}] = findpeaks((-ydb2{1,i}(1,:)),t,'MinPeakHeight',15,'MinPeakDistance',150);
    
    figure(i)
    clf
    subplot(1,2,1);
    plot(t,ydb2{1,i}(1,:))
    hold on
    plot(locs{1,i},peaks{1,i},'ro')
    hold on
    plot(slocs{1,i},-speaks{1,i},'bo')
    grid
    if i==1
        title('Lead 1');
    end
    if i==2
        title('Lead 2');
    end
    if i==3
        title('Lead 3');
    end
    if i==4
        title('Lead 4');
    end
    if i==5
        title('Lead 5');
    end
    if i==6
        title('Lead 6');
    end
    if i==7
        title('Lead 7');
    end
    if i==8
        title('Lead 8');
    end
    if i==9
        title('Lead 9');
    end
    if i==10
        title('Lead 10');
    end
    if i==11
        title('Lead 11');
    end
    if i==12
        title('Lead 12');
    end    
    %legend('ECG Signal','R-waves','S-waves','Location','northeastoutside')
    xlabel('Seconds')

    subplot(1,2,2);
    plot(t, ecg(:,i+1));
     grid
    xlabel('Seconds')
    i = i+1;
    
end 

% Baseline drift elimination:
% Output is 'D' the detrended signal for 12 leads, after filtering
% C is the cellularized matrix of the 12 leads
filt_lead_mat = cell2mat(filtLead);
C=num2cell(filt_lead_mat,1); %call C{1}...C{12}

%Initialize detrended lead data matrix 'D' to be zeros
D = zeros(L, 12);
%f_y = zeros(L, 12);
% Baseline drift greater than 2.5mV in any lead
j=1;

while j < 13
    [p,s,mu] = polyfit((1:numel(C{j}))',C{j},6);
    f_y(:,j) = polyval(p,(1:numel(C{j}))',[],mu);
    
    D(:,j) = C{j} - f_y(:,j);   % Detrend data using polyfit curve f_y
    
    %Baseline drift flag check
    if f_y(:,j) > 2500
        F_Baseline = 1;
    end
    
    j=j+1;
end

%Fast fourier transform to find signal frequency properties pre filtering
i=2;

while i < 14
    fresult = fft(F{i});
    P2 = abs(fresult/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = fs*(0:(L/2))/L;  % frequency domain f
    
    i=i+1;
end
%Fast fourier transform to find signal frequency prpoerties after filtering
i=1;

while i < 13
    fresult_new = fft(filtLead{1,i});
    P2 = abs(fresult_new/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = fs*(0:(L/2))/L;  % frequency domain f
    
    i=i+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Processing
%Misplaced electrode, low and high amplitude, motion artifacts,
%flat line detection
%Input should be D - detended/filtered data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Indices set
i = 1;
j = 1;
k = 1;

%Initialize matrices to zeros
max_ampl = zeros(1,12);
min_ampl = zeros(1,12);
max_sum = zeros(1,12);
max_avg = zeros(1,12);
min_sum = zeros(1,12);
min_avg = zeros(1,12);

%Lead Initialization
%filtLead = lead data after all filtering
l_1=D(:,1);    %Lead I: Small Q wave, medium R, small S, upright T
l_2=D(:,2);    %Lead II: No Q wave, large R, small S, upright T 
l_3=D(:,3);    %Lead III: Small Q wave, variable R, variable S, variable T
l_4=D(:,4);    %aVR: Variable Q wave, small R, large S, inverted T 
l_5=D(:,5);    %aVL: Variable Q wave, variable R wave, none to large S, variable T
l_6=D(:,6);    %aVF: small Q wave, small R, variable S, variable T
l_7=D(:,7);    %V1: QS complex, small R, large S, variable T
l_8=D(:,8);    %V2: No Q wave, larger R than V1, large S, unpright T
l_9=D(:,9);    %V3: No Q wave, variable S, variable R, upright T
l_10=D(:,10);  %V4: No Q wave, larger R than V3, smaller S than V3, upright T
l_11=D(:,11);  %V5: Small Q wave, larger R than V4, smaller S than V4, upright T
l_12=D(:,12);  %V6: small Q wave, smaller R than V5, smaller S than V5, upright T

% Saturation (high amplitude) condition and flag setting:
% If any lead is considered high amplitude (greater than 2mV for excerpt of continuous 200 msec ),
% F_max is set to 1, and 12-lead ECG is considered unacceptable due to motion artifacts or EMG
% noise.
x=0;
while k < 13
    while x < 50
        count_1 = 1+(x*100);
        count_2 = 100+(x*100);
        %column_x = x+1;
        
        %200msec = 100 samples out of 5000 samples per D value
        %s_mat = D(1:100, 1);    %trial
        submatrix_D = D(count_1:count_2, k);
        max_submatrix_D = max(submatrix_D);
        
        if (max_submatrix_D >= 1000) %2000uV = 2mV
            F_Max = 1;
            F_EMG = 1;
        end
        x=x+1;
    end
    x=0;
    k=k+1;
    
end

% Flat line condition and flag setting:
% If any lead is considered to to a constant 0 voltage for at least 1sec
% length, F_Flat is set to 1, and 12-lead ECG is considered unacceptable.
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
i=1;
divisor=12;
%Set peaks and speaks empty cells to 0 element to fix error
while k < 13
    if isempty(peaks{1,k} == 1)
        peaks{1,k}=0;
    end
    if isempty(speaks{1,k} == 1)
        speaks{1,k}=0;
    end
    k=k+1;
end
k=1;
%Set dividing factor to refelct only non-zero cell elements of peaks and
%speaks
while i < 13
   if max(peaks{1,i}) == 0
        divisor = divisor - 1;
   end
   i=i+1;
end        
%Search for maximum of peak values in peak[] of 12 leads
while k < 13    
    max_ampl(k) = max(peaks{1,k});
    min_ampl(k) = min(peaks{1,k});
    max_ampl_s(k) = max(speaks{1,k});
    min_ampl_s(k) = min(speaks{1,k});
    
    % Voltage averaging of peaks for 12 leads    
    max_sum(k) = sum(peaks{1,k});
    max_avg(k) = max_sum(k)/divisor;
    min_sum(k) = sum(peaks{1,k});
    min_avg(k) = min_sum(k)/divisor;
    
    max_sum_s(k) = sum(speaks{1,k});
    max_avg_s(k) = max_sum_s(k)/divisor;
    min_sum_s(k) = sum(speaks{1,k});
    min_avg_s(k) = min_sum_s(k)/divisor;
    
    k=k+1;
end
k=1;

% Low amplitude condition and flag setting:
% If three leads are considered low amplitude (less than 125uV), F_Min is set to 1, and
% 12-lead ECG is considered unacceptable due to possible poor skin-electrode
% contact.
checklow = 0;
while k < 13
    min_length = length(peaks{k});
    if max_ampl(k) <= 15 || min_length < 3  %thresholding check: if not enough peaks found by peak finder or min ampl <= 0.125mV = low saturation
        checklow = checklow+1;
    end
    if checklow == 3
        F_Min = 1;
    end
    k=k+1;
end                                
        
% EMG Noise condition:
% If peak[] length is greater than 100; peaks too frequent to be bpm; set
% EMG noise flag to 1.
% Heart rate between range of 30bpm - 290bpm is reasonably within the human range;
% Resting heartrate of 100bpm and greater can however indicate the patient is tachycardic.
i=1;
while i<13
    fresult_new = fft(filtLead{1,i});
    
    P2 = abs(fresult_new/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    freq_max = max(P1);
    freq_min =min(P1);
    freq_peak = freq_max - freq_min; %WRONG
    peak_num = length(peaks{1,i});
    
    if (peak_num > 50) ||  freq_peak > 175 % 290bpm and greater, or frequency amplitude of lead is above 40
        F_EMG = 1;
    end
    if (length(peaks{1,i}) < 5)
        F_Contact = 1;
    end
    i=i+1;
end

% Reversed RA and LL limb lead check: Inverted P-QRS in Lead II (lead shows
% difference between LL and RA, directed towards LL at 60deg)
% If s-peaks > r-peaks from QRS detection with wavelet, and Leads avR and aVF switch places (I now greater than III), signal is inverted
if (max_avg(2) < max_avg_s(2)) && (max_avg(4) > max_avg(6)) && (max_avg(1)<max_avg_s(1))
    F_RA_LL = 1;
end

% Reversed RA and LA limb lead check: Inverted P-QRS in Lead I (lead shows
% difference between RA and LA, directed towards LA at 0deg), and often
% Lead aVR becomes positive.
% If s-peaks > r-peaks in Lead I AND aVR becomes positive, F_RA_LA set to 1
% May simulate dextrocardia
if (max_avg(1) < max_avg_s(1)) && (max_avg(4) > max_avg_s(4))
    F_RA_LA = 1;
end

% Beats per minute calculation:
% Patients heartbeat, frequency of beats- in gui
t_total = max(t1);  %msec; generally 10sec but can be smaller
sum_bpm = 0;
avg_bpm = 0;
i = 1;
while i < 13
    peaks_i_length = length(peaks{1,i});
    sum_bpm = sum_bpm + peaks_i_length;  %number of peaks for 6*10sec sample = extension to 60sec of data
    i=i+1;
end
avg_bpm = (sum_bpm*6)/12; %averaging of all 12 leads number of peaks, for contingency

    
    if(F_RA_LA || F_RA_LL || F_EMG || F_Flat || F_Min || F_Max || F_Baseline || ~F_Contact == 1)
        fprintf('This ECG data is unacceptable \n');
        F_Acceptable = 0;
        if(F_RA_LA)
            fprintf('F_RA_LA\n');
        end

        if(F_RA_LL)
            fprintf('F_RA_LL\n');
        end

        if(F_EMG)
            fprintf('F_EMG\n');
        end

        if(F_Flat )
            fprintf('F_Flat\n');
        end

        if(F_Min)
            fprintf('F_Min\n');
        end

        if(F_Max)
            fprintf('F_Max\n');
        end

        if(F_Baseline)
            fprintf('F_Baseline\n');
        end

        if(F_Contact)
            fprintf('F_Contact\n');
        end    

    else
        fprintf('This ECG data is acceptable \n');
        F_Acceptable = 1;
        pass = pass + 1;
    end
    
    results(l,:) =  [F_Acceptable F_RA_LA F_RA_LL F_EMG F_Flat F_Max F_Baseline ~F_Contact];

l = l + 1;

end
  
    