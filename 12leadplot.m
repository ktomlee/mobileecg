%plot(VarName14, VarName15)

figure 
subplot(12,1,1)
plot(VarName14, VarName15)
title('12 Lead ECG')
ylabel('Lead 1 Voltage (mV)', 'fontsize', 8)
yh = get(gca,'YLabel'); % Handle of the y label
set(yh, 'Units', 'Normalized')
pos = get(yh, 'Position');
set(yh, 'Position',pos.*[1,1,1],'Rotation',0)
set(gca,'xtick',[],'ytick',[])
set(gca,'ytickMode', 'auto')
subplot(12,1,2)
plot(VarName14, VarName16)
ylabel('Lead 2', 'fontsize', 8)
yh = get(gca,'YLabel'); % Handle of the y label
set(yh, 'Units', 'Normalized')
pos = get(yh, 'Position');
set(yh, 'Position',pos.*[1,1,1],'Rotation',0)
set(gca,'xtick',[],'ytick',[])
subplot(12,1,3)
plot(VarName14, VarName17)
ylabel('Lead 3', 'fontsize', 8)
yh = get(gca,'YLabel'); % Handle of the y label
set(yh, 'Units', 'Normalized')
pos = get(yh, 'Position');
set(yh, 'Position',pos.*[1,1,1],'Rotation',0)
set(gca,'xtick',[],'ytick',[])
subplot(12,1,4)
plot(VarName14, VarName18)
ylabel('Lead 4', 'fontsize', 8)
yh = get(gca,'YLabel'); % Handle of the y label
set(yh, 'Units', 'Normalized')
pos = get(yh, 'Position');
set(yh, 'Position',pos.*[1,1,1],'Rotation',0)
set(gca,'xtick',[],'ytick',[])
subplot(12,1,5)
plot(VarName14, VarName19)
ylabel('Lead 5', 'fontsize', 8)
yh = get(gca,'YLabel'); % Handle of the y label
set(yh, 'Units', 'Normalized')
pos = get(yh, 'Position');
set(yh, 'Position',pos.*[1,1,1],'Rotation',0)
set(gca,'xtick',[],'ytick',[])
subplot(12,1,6)
plot(VarName14, VarName20)
ylabel('Lead 6', 'fontsize', 8)
yh = get(gca,'YLabel'); % Handle of the y label
set(yh, 'Units', 'Normalized')
pos = get(yh, 'Position');
set(yh, 'Position',pos.*[1,1,1],'Rotation',0)
set(gca,'xtick',[],'ytick',[])
subplot(12,1,7)
plot(VarName14, VarName21)
ylabel('Lead 7', 'fontsize', 8)
yh = get(gca,'YLabel'); % Handle of the y label
set(yh, 'Units', 'Normalized')
pos = get(yh, 'Position');
set(yh, 'Position',pos.*[1,1,1],'Rotation',0)
set(gca,'xtick',[],'ytick',[])
subplot(12,1,8)
plot(VarName14, VarName22)
ylabel('Lead 8', 'fontsize', 8)
yh = get(gca,'YLabel'); % Handle of the y label
set(yh, 'Units', 'Normalized')
pos = get(yh, 'Position');
set(yh, 'Position',pos.*[1,1,1],'Rotation',0)
set(gca,'xtick',[],'ytick',[])
subplot(12,1,9)
plot(VarName14, VarName23)
ylabel('Lead 9', 'fontsize', 8)
yh = get(gca,'YLabel'); % Handle of the y label
set(yh, 'Units', 'Normalized')
pos = get(yh, 'Position');
set(yh, 'Position',pos.*[1,1,1],'Rotation',0)
set(gca,'xtick',[],'ytick',[])
subplot(12,1,10)
plot(VarName14, VarName24)
ylabel('Lead 10', 'fontsize', 8)
yh = get(gca,'YLabel'); % Handle of the y label
set(yh, 'Units', 'Normalized')
pos = get(yh, 'Position');
set(yh, 'Position',pos.*[1,1,1],'Rotation',0)
set(gca,'xtick',[],'ytick',[])
subplot(12,1,11)
plot(VarName14, VarName25)
ylabel('Lead 11', 'fontsize', 8)
yh = get(gca,'YLabel'); % Handle of the y label
set(yh, 'Units', 'Normalized')
pos = get(yh, 'Position');
set(yh, 'Position',pos.*[1,1,1],'Rotation',0)
set(gca,'xtick',[],'ytick',[])
subplot(12,1,12)
plot(VarName14, VarName26)
ylabel('Lead 12', 'fontsize', 8)
yh = get(gca,'YLabel'); % Handle of the y label
set(yh, 'Units', 'Normalized')
pos = get(yh, 'Position');
set(yh, 'Position',pos.*[1,1,1],'Rotation',0)
set(gca,'xtick',[],'ytick',[])
set(gca,'xtickMode', 'auto')
xlabel('Time (msec)')
