q = load('1002867.txt');

EKG1 = q(:,2);
EKG2 = q(:,3);
Fs   = 500;
t    = [0:length(EKG1)-1]/Fs;

[R1,TR1]  = findpeaks( EKG1, t, 'MinPeakHeight',150);
[Q1,TQ1]  = findpeaks((-EKG1), t, 'MinPeakHeight',100);              % NOTE: No ?S? Waves In EKG1
[R2,TR2]   = findpeaks( EKG2, t, 'MinPeakHeight', 100);
[QS2,TQS2] = findpeaks((-EKG2), t, 'MinPeakHeight', 75);
figure(1)
subplot(2,1,1)
plot(t, EKG1)
hold on
plot(TR1, R1, '^r')
plot(TQ1, -Q1, 'vg')
hold off
grid
axis([0  10    ylim])
legend('EKG', 'R', 'Q')

subplot(2,1,2)
plot(t, EKG2)
hold on
plot(TR2, R2, '^r')
plot(TQS2, -QS2, 'vg')
plot(TQS2, -QS2, 'vb')
hold off
grid
axis([0  10    ylim])
legend('EKG', 'R', 'Q', 'S')