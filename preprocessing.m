%Pre-processing

ecg = load('1002867.txt');

L = length(ecg);
Lead=num2cell(ecg,1);

fs = 500;
t = (1:L);
t1=(0:L-1)/fs;
Fv = linspace(0, 1, fix(L/2)+1)*(fs/2);
Iv = 1:length(Fv);

%row, column
leadI = ecg((1:L),2);
leadII = ecg((1:L),3);
leadIII = ecg((1:L),4);

order = 4;
fcutlow=0.5;   %low cut frequency in Hz
fcuthigh=100;   %high cut frequency in Hz
[b,a]=butter(order,[fcutlow,fcuthigh]/(fs/2),'bandpass');

ftleadI = filter(b,a,leadI);  %BP filter time domain
ftleadII = filter(b,a,leadII);  %BP filter time domain

i = 2;

filtLead={};

while i < 14
    filtLead{i-1} = filter(b,a,C{i});
    i = i+1;
end

fleadI = fft(leadI)/L;
ffleadI = filter(b,a,fleadI);  %BP filter freq domain

w = (60/(fs/2));
bw = w;
[num,den]=iirnotch(w,bw); % notch filter implementation 

fftleadI = filter(num, den, ftleadI); % notch filter time
fftleadII = filter(num, den, ftleadII); % notch filter time

fffleadI = filter(num, den, ffleadI); % notch filter freq

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
        %if((filtLead{1,i}(j,1) > filtLead{1,i}(j+1,1)) && (filtLead{1,i}(j,1) > filtLead{1,i}(j-1,1)) && ) 
            peak{1,i}(j,1) = filtLead{1,i}(j,1);                                   
            time{1,i}(j,1)=(j-1)/fs;           %position stored where peak value met;              
        end
        j=j+1;
    end
  
    peak{1,i}(peak{1,i} == 0) =[];
    time{1,i}(time{1,i}==0)=[];
    
    i = i+1;
end

hhI = fftleadI;
hhII = fftleadII;
peakI=[];           %loop initialing, having all the value zero in the array
peakII=[];
timeI=0;            %loop initialing, having all the value zero in the array
timeII=0;
thI=0.6*max(fftleadI);  %thresold setting at 60 percent of maximum value
thII=0.6*max(fftleadII);



for i=2:L-1   
    if((hhI(i)>hhI(i+1))&&(hhI(i)>hhI(i-1))&&(hhI(i)>thI))  
        peakI(i)=hhI(i);                                   
        timeI(i)=[i-1]/fs;           %position stored where peak value met;              
    end
    
    if((hhII(i)>hhII(i+1))&&(hhII(i)>hhII(i-1))&&(hhII(i)>thII))  
        peakII(i)=hhII(i);                                   
        timeII(i)=[i-1]/fs;           %position stored where peak value met;              
    end
end
peakI(peakI==0)=[];           % neglect all zeros from array;
peakII(peakII==0)=[];

timeI(timeI==0)=[];     % neglect all zeros from array;
timeII(timeII==0)=[]; 


figure(1);
i = 1;
while i < 13
   subplot(12,1,i);
   plot(t1, filtLead{i});
   hold on;
   plot(time{i}, peak{i}, '*r');
   i = i+1;
end

%plot(Fv,abs(fffleadI(Iv))*2);
%title('Freq Notch + BP Filtered');


