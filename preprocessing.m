%Pre-processing

ecg = load('1002867.txt');

L = length(ecg);
Lead=num2cell(ecg,1);

fs = 500;
t = (1:L);
t1=(0:L-1)/fs;
Fv = linspace(0, 1, fix(L/2)+1)*(fs/2);
Iv = 1:length(Fv);

order = 6;
fcutlow=0.5;   %low cut frequency in Hz
fcuthigh=100;   %high cut frequency in Hz
[b,a]=butter(order,[fcutlow,fcuthigh]/(fs/2),'bandpass');

bandstop = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',fs);

bandpass = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',0.5,'HalfPowerFrequency2',100, ...
    'SampleRate',fs);
        

w = (60/(fs/2));
bw = w;
[num,den]=iirnotch(w,bw); % notch filter implementation 

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
    threshold(i) = 0.6*max(filtLead{1,i});
    
    for j = 2: L-1
        if((filtLead{1,i}(j,1) > filtLead{1,i}(j+1,1))&&(filtLead{1,i}(j,1) > filtLead{1,i}(j-1,1)) && (filtLead{1,i}(j,1) > threshold(i)))  
            peak{1,i}(j,1) = filtLead{1,i}(j,1);                                   
            time{1,i}(j,1)=(j-1)/fs;           %position stored where peak value met;              
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

%plot(Fv,abs(fffleadI(Iv))*2);
%title('Freq Notch + BP Filtered');


