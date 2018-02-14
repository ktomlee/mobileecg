%fileID=fopen('1002867.txt');
%B=fread(fileID);
ecg=load('1002867.txt');
ecg=ecg(:,2);
Fs=1000;
N=length(ecg);
t=(0:N-1)/Fs; %time period(total sample/Fs )

subplot(3,1,1)
plot(t, VarName2)
ylabel('Amplitude')
title('Raw Single Lead ECG Data')

w=50/(Fs/2);
bw=w;
[num,den]=iirnotch(w,bw); % notch filter implementation 
ecg_notch=filter(num,den,ecg);
[e,f]=wavedec(ecg_notch,10,'db6');% Wavelet implementation
g=wrcoef('a',e,f,'db6',8); 
ecg_wave=ecg_notch-g; % subtracting 10th level aproximation signal
                       %from original signal                  
ecg_smooth=smooth(ecg_wave); % using average filter to remove glitches
                             %to increase the performance of peak detection 
N1=length(ecg_smooth);
t1=(0:N1-1)/Fs;

% Peak detection algorithm 
% For more detailsor detailed explanation on this look into 
% Matlab for beginers 
hh=ecg_smooth;
 j=[];           %loop initialing, having all the value zero in the array
time=0;          %loop initialing, having all the value zero in the array
th=0.6*max(hh);  %thresold setting at 45 percent of maximum value
 
for i=2:N1-1 % length selected for comparison  
    % deopping first ie i=1:N-1  point because hh(1-1) 
   % in the next line  will be zero which is not appreciable in matlab 
    if((hh(i)>hh(i+1))&&(hh(i)>hh(i-1))&&(hh(i)>th))  
% condition, i should be> then previous(i-1),next(i+1),thrsold point;
        j(i)=hh(i);                                   
%if condition satisfy store hh(i)in place of j(i)value whichis initially 0;
       
        time(i)=[i-1]/Fs;           %position stored where peak value met;              
      
    end
end
 j(j==0)=[];               % neglect all zeros from array;
 time(time==0)=[];     % neglect all zeros from array;
m=(time)';               % converting rows in column;
k=length(m);

z = gausswin(50);
y = filter(z,1,ecg);
subplot(3,1,2);
plot(t,y); title('Gaussian Smoothed Data');

 % Plot frequency spectrum
subplot(3,1,3);
plot(t1,ecg_smooth)
hold on;                 % hold the plot and wait for next instruction;
plot(time,j,'*r'); title('Peak Detection on Averaged Data')    
xlabel('Time (sec)')

