%Pre-processing

ecg = load('1002867.txt');

L = length(ecg);
fs = 500;
fn = 250;
t = (1:L);
Fv = linspace(0, 1, fix(L/2)+1)*fn;
Iv = 1:length(Fv);


%row, column
leadI = ecg((1:L),2);
leadII = ecg((1:L),3);
leadIII = ecg((1:L),4);

fleadI = fft(leadI)/L;

%plot(t,leadI);
plot(Fv, abs(fleadI(Iv))*2)

