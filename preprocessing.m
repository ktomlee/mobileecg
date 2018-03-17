%Pre-processing

ecg = load('1002867.txt');

L = length(ecg);
fs = 500;
fn = 250;
t = (1:L);
t1=(0:L-1)/fs;
Fv = linspace(0, 1, fix(L/2)+1)*fn;
Iv = 1:length(Fv);

%row, column
leadI = ecg((1:L),2);
leadII = ecg((1:L),3);
leadIII = ecg((1:L),4);

order = 4;
fcutlow=0.5;   %low cut frequency in Hz
fcuthigh=100;   %high cut frequency in Hz
[b,a]=butter(order,[fcutlow,fcuthigh]/(fs/2),'bandpass');

ftleadI = filter(b,a,leadI);
fleadI = fft(leadI)/L;
ffleadI = filter(b,a,fleadI);

w = (60/(fs/2));
bw = w;
[num,den]=iirnotch(w,bw); % notch filter implementation 

fftleadI = filter(num, den, ftleadI); % notch filter time
fffleadI = filter(num, den, ffleadI); % notch filter freq

hh = fftleadI;
j=[];           %loop initialing, having all the value zero in the array
time=0;          %loop initialing, having all the value zero in the array
th=0.6*max(fftleadI);  %thresold setting at 45 percent of maximum value

for i=2:L-1 % length selected for comparison  
    % dropping first ie i=1:L-1  point because hh(1-1) 
   % in the next line  will be zero which is not appreciable in matlab 
    if((hh(i)>hh(i+1))&&(hh(i)>hh(i-1))&&(hh(i)>th))  
% condition, i should be> then previous(i-1),next(i+1),thrsold point;
        j(i)=hh(i);                                   
%if condition satisfy store hh(i)in place of j(i)value whichis initially 0;
       
        time(i)=[i-1]/fs;           %position stored where peak value met;              
      
    end
end
 j(j==0)=[];               % neglect all zeros from array;
 time(time==0)=[];     % neglect all zeros from array;
m=(time)';               % converting rows in column;

figure(1);
plot(t1,fftleadI);
hold on;                 % hold the plot and wait for next instruction;
plot(time,j,'*r'); 
title('Peak Detection on Filtered Data');
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


