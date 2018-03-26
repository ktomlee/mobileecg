ecg = load('1002867.txt');

L = length(ecg);
fs = 500;
t=(0:L-1)/(fs);

figure(1);
subplot(3,1,1);
plot(t, ecg(:,2));
subplot(3,1,2);
plot(t, ecg(:,3));
subplot(3,1,3);
plot(t, ecg(:,4));

[CA,CD] = dwt(ecg(:,2),'db1');
ecgrec = zeros(size(ecg(:,2)));

db2 = modwt(ecg(:,3),'db2',5);
db2rec = zeros(size(db2));
db2rec(4:5,:) = db2(4:5,:);
y = imodwt(db2rec,'db2');

%y = abs(y).^2;
[qrspeaks,locs] = findpeaks(y,t,'MinPeakHeight',50,'MinPeakDistance',0.300);
figure(2)
subplot(1,2,1);
plot(t,y)
hold on
plot(locs,qrspeaks,'ro')
xlabel('Seconds')

subplot(1,2,2);
plot(t, ecg(:,3));
xlabel('Seconds')


%plot(t,CD);