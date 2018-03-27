ecg = load('1063073.txt');

L = length(ecg);
fs = 500;
t=(0:L-1)/(fs);

db2 = {};
db2rec = {};
ydb2 = {};
peaks = {};
locs = {};
i = 1;

while i < 13

    db2{1,i} = modwt(ecg(:,i+1),'db2',5);
    db2rec{1,i} = zeros(size(db2{1,i}));
    db2rec{1,i}(4:5,:) = db2{1,i}(4:5,:);
    ydb2{1,i} = imodwt(db2rec{1,i},'db2');

    [peaks{1,i},locs{1,i}] = findpeaks(ydb2{1,i}(1,:),t,'MinPeakHeight',30,'MinPeakDistance',0.300);
    [speaks{1,i},slocs{1,i}] = findpeaks((-ydb2{1,i}(1,:)),t,'MinPeakHeight',30,'MinPeakDistance',0.300);
    
    figure(i)
    clf
    subplot(1,2,1);
    plot(t,ydb2{1,i}(1,:))
    hold on
    plot(locs{1,i},peaks{1,i},'ro')
    hold on
    plot(slocs{1,i},-speaks{1,i},'bo')
    grid
    %legend('ECG Signal','R-waves','S-waves','Location','northeastoutside')
    xlabel('Seconds')

    subplot(1,2,2);
    plot(t, ecg(:,i+1));
     grid
    xlabel('Seconds')
    
    if(size(peaks{1,i}) < 3)
        fprintf('Not enough peaks were detected at %i\n', i);
    end
    
    if(max(peaks{1,i}) > 500)
        fprintf('High amplitude detected at %i\n', i);
    end
    
    i = i + 1;
    
end
