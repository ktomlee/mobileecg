%Pre-processing

ecg = load('1002867.txt');

L = length(ecg);
Lead=num2cell(ecg,1);

fs = 500;
t = (1:L);
t1=(0:(L/2)-1)/fs;
Fv = linspace(0, 1, fix(L/2)+1)*(fs/2);
Iv = 1:length(Fv);



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
db2Lead={};

while i < 14
    
    [Lo_D,Hi_D] = wfilters('haar','d'); 
    [ca1,cd1] = dwt(Lead{1,i},Lo_D,Hi_D);
    db2Lead{1,i-1} = ca1;
    i = i+1;
end

i = 1;
peak={};
time={};
threshold=[];
while i < 13
    peak{1,i} = [];
    time{1,i} = [];
    threshold(i) = 0.5*max(db2Lead{1,i});
    
    for j = 2: (L/2)-1
        if((db2Lead{1,i}(j,1) > db2Lead{1,i}(j+1,1))&&(db2Lead{1,i}(j,1) > db2Lead{1,i}(j-1,1)) && (db2Lead{1,i}(j,1) > threshold(i)))  
            peak{1,i}(j,1) = db2Lead{1,i}(j,1);                                   
            time{1,i}(j,1)=(j-1)/(fs);           %position stored where peak value met;              
        end
        j=j+1;
    end
  
    peak{1,i}(peak{1,i} == 0) =[];
    time{1,i}(time{1,i} == 0)=[];
    
    i = i+1;
end

figure(1);
i = 1;
while i < 2
   subplot(12,1,i);
   plot(t1, db2Lead{i});
   hold on;
   plot(time{i}, peak{i}, '*r');
   i = i+1;
end

%plot(Fv,abs(fffleadI(Iv))*2);
%title('Freq Notch + BP Filtered');


