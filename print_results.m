ecg = load('2846714.txt');

% Signal Variables
L = length(ecg);
fs = 500;
t1=(0:L-1)/fs;
Fv = linspace(0, 1, fix(L/2)+1)*(fs/2);
Iv = 1:length(Fv);
%Lead=num2cell(ecg,1);

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

figure(1)
subplot(6,2,1);
plot(t1, l_1);
title('Lead I','FontSize', 8)
axis([0 10 -500 500]);
ylabel('Voltage (mV)')
grid on
grid minor
subplot(6,2,2);
plot(t1, l_2);
title('Lead II','FontSize', 8)
axis([0 10 -500 500]);
grid minor
subplot(6,2,3);
plot(t1, l_3);
title('Lead III','FontSize', 8)
axis([0 10 -500 500]);
grid minor
subplot(6,2,4);
plot(t1, l_4);
title('Lead aVR','FontSize', 8)
axis([0 10 -500 500]);
grid minor
subplot(6,2,5);
plot(t1, l_5);
title('Lead aVL','FontSize', 8)
axis([0 10 -500 500]);
grid minor
subplot(6,2,6);
plot(t1, l_6);
title('Lead aVF','FontSize', 8)
axis([0 10 -500 500]);
grid minor
subplot(6,2,7);
plot(t1, l_7);
title('Lead V1','FontSize', 8)
axis([0 10 -500 500]);
grid minor
subplot(6,2,8);
plot(t1, l_8);
title('Lead V2','FontSize', 8)
axis([0 10 -500 500]);
grid minor
subplot(6,2,9);
plot(t1, l_9);
title('Lead V3','FontSize', 8)
axis([0 10 -500 500]);
grid minor
subplot(6,2,10);
plot(t1, l_10);
title('Lead V4','FontSize', 8)
axis([0 10 -500 500]);
grid minor
subplot(6,2,11);
plot(t1, l_11);
title('Lead V5','FontSize', 8)
axis([0 10 -500 500]);
xlabel('Time (sec)')
grid minor
subplot(6,2,12);
plot(t1, l_12);
title('Lead V6','FontSize', 8)
axis([0 10 -500 500]);
grid minor
xlabel('Time (sec)')
    

