%ecg = load('1002603.txt'); %Normal
%ecg = load('1002867.txt'); %Normal
%ecg = load('1003574.txt'); %Random signals
%ecg = load('1004502.txt');
%ecg = load('1005639.txt'); %good, maybe low amplitude
%ecg = load('1006637.txt');
%ecg = load('1007823.txt');
%ecg = load('1009404.txt');

% Signal Variables
L = length(ecg);
fs = 500;
t = (1:L);  %Don't use
t1=(0:L-1)/fs;

l_1=ecg(:,2);    %Lead I: Small Q wave, medium R, small S, upright T
l_2=ecg(:,3);    %Lead II: No Q wave, large R, small S, upright T 
l_3=ecg(:,4);    %Lead III: Small Q wave, variable R, variable S, variable T
l_4=ecg(:,5);    %aVR: Variable Q wave, small R, large S, inverted T 
l_5=ecg(:,6);    %aVL: Variable Q wave, variable R wave, none to large S, variable T
l_6=ecg(:,7);    %aVF: small Q wave, small R, variable S, variable T
l_7=ecg(:,8);    %V1: QS complex, small R, large S, variable T
l_8=ecg(:,9);    %V2: No Q wave, larger R than V1, large S, unpright T
l_9=ecg(:,10);    %V3: No Q wave, variable S, variable R, upright T
l_10=ecg(:,11);  %V4: No Q wave, larger R than V3, smaller S than V3, upright T
l_11=ecg(:,12);  %V5: Small Q wave, larger R than V4, smaller S than V4, upright T
l_12=ecg(:,13);  %V6: small Q wave, smaller R than V5, smaller S than V5, upright T

%12 Lead plot
figure
subplot(6,2,1);
plot(t1, l_1);
title('12 Lead ECG Plot')
subplot(6,2,2);
plot(t1, l_2);
subplot(6,2,3);
plot(t1, l_3);
ylabel('Voltage (mV)')
subplot(6,2,4);
plot(t1, l_4);
subplot(6,2,5);
plot(t1, l_5);
subplot(6,2,6);
plot(t1, l_6);
subplot(6,2,7);
plot(t1, l_7);
subplot(6,2,8);
plot(t1, l_8);
subplot(6,2,9);
plot(t1, l_9);
subplot(6,2,10);
plot(t1, l_10);
subplot(6,2,11);
plot(t1, l_11);
subplot(6,2,12);
plot(t1, l_12);
xlabel('Time (sec)')
