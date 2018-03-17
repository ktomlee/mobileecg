%Pre-processing

ecg = load('1002867.txt');

L = length(ecg);
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

fleadI = fft(leadI)/L;
ffleadI = filter(b,a,fleadI);  %BP filter freq domain

w = (60/(fs/2));
bw = w;
[num,den]=iirnotch(w,bw); % notch filter implementation 

fftleadI = filter(num, den, ftleadI); % notch filter time
fftleadII = filter(num, den, ftleadII); % notch filter time

fffleadI = filter(num, den, ffleadI); % notch filter freq

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
subplot(2,1,1);
plot(t1,fftleadI);
hold on;                 % hold the plot and wait for next instruction;
plot(timeI,peakI,'*r'); 
title('Peak Detection on Filtered LeadI');
xlabel('Time (sec)');

subplot(2,1,2);
plot(t1,fftleadII);
hold on;                 % hold the plot and wait for next instruction;
plot(timeII,peakII,'*r'); 
title('Peak Detection on Filtered LeadII');
xlabel('Time (sec)');


%figure(2);
%subplot(3,2,1);
%plot(t,leadI);
%title('Time Not Filtered');

%subplot(3,2,3);
%plot(t, ftleadI);
%title('Time BP Filtered');

%subplot(3,2,5);
%plot(t, fftleadI);
%title('Notch + BP Filtered');
%subplot(3,2,2);
%plot(Fv,abs(fleadI(Iv))*2);
%title('Freq Not Filtered');
%subplot(3,2,4);
%plot(Fv,abs(ffleadI(Iv))*2);
%title('Freq BP Filtered');
%subplot(3,2,6);
%plot(Fv,abs(fffleadI(Iv))*2);
%title('Freq Notch + BP Filtered');


