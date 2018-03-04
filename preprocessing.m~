%Pre-processing

ecg = load('1002867.txt');

N = (ecg);
L = length(N);
fs = 500;
%t = (1:N);
t = (0:L-1)/fs;

i=2;
j=2;
D = zeros(L, 13);

%row, column
leadI = ecg((1:N),2);
plot(leadI);

C=num2cell(N,1); %call C{1}, C{2} etc... up to 13 (C{1} is time counter)

d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',fs);
figure(1)
while i < 14
    filt_C = filtfilt(d,C{i});
    %plot(t,C{i},t,filt_C);
    plot(t, C{i});
    title 'Leads Plot, Filtered 60Hz', ylabel 'Voltage (uV)', xlabel 'Time (sec)';
    hold on
    i=i+1;
end

figure(2)
while j < 14
    [p,s,mu] = polyfit((1:numel(C{j}))',C{j},6);
    f_y = polyval(p,(1:numel(C{j}))',[],mu);
    
    D = C{j} - f_y;        % Detrend data
    plot(t, D);
    title('Detrended ECG Signal'), xlabel('Time (sec)'), ylabel('Voltage(uV)')
    hold on
    j=j+1;
end

%60Hz filter response
fvtool(d,'Fs',fs)

%<0.5Hz filter and response

%Peak detection - check for max and min amplitudes ; 
%Minimum peak detection -for inverted leads?

%Flat line detection

%QRS detection
%Checks in place for 3 common problems with electrodes switch/placement
%(research another)

% Patients heartbeat, frequency of beats- in gui

