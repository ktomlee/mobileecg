%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QRS Detection using thresholding to find peaks of interest
% 100 = min peak height of R wave; 300 = less than minimum distance between
% subsequent peaks

ecg = load('1002867.txt');
fs = 500;
L= length(ecg);
time = (0:L-1)/fs;
Lead=num2cell(ecg,1); %lead 1 = Lead{2}

%INitialize leads from cells of 'Lead'
lead1 = Lead{2};
lead2 = Lead{3};
lead3 = Lead{4};
lead4 = Lead{5};
lead5 = Lead{6};
lead6 = Lead{7};
lead7 = Lead{8};
lead8 = Lead{9};
lead9 = Lead{10};
lead10 = Lead{11};
lead11 = Lead{12};
lead12 = Lead{13};

%i = 1;
%while i < 2
    [~,locs_Rwave] = findpeaks(lead1,'MinPeakHeight',100,...
                                        'MinPeakDistance',300);
    
    f_inverted = -lead1;
    [~,locs_Swave] = findpeaks(f_inverted,'MinPeakHeight',30,...
                                        'MinPeakDistance',300);
                                    
    figure
    hold on 
    plot(lead1)
    plot(locs_Rwave,lead1(locs_Rwave),'rv','MarkerFaceColor','r')
    plot(locs_Swave,lead1(locs_Swave),'rs','MarkerFaceColor','b')
    grid on
    legend('ECG Signal','R-waves','S-waves')
    xlabel('Samples')
    ylabel('Voltage(mV)')
    title('R-wave and S-wave in Noisy ECG Signal')
    
    
%end
