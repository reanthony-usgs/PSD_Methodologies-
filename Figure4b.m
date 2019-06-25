% Creates the difference mean and median PSD figure for earthquake removed
% PSD datasets using both 1 and 3 hour windows

clear all

% load in the data

load EQ_Differences/1Hr_Mean_Diff.mat
Hr1_Mean = Mean_Diff;
load EQ_Differences/1Hr_Median_Diff.mat
Hr1_Median = Median_Diff;
load EQ_Differences/Periods.mat
load EQ_Differences/3Hr_Mean_Diff.mat
Hr3_Mean = Mean_Diff;
load EQ_Differences/3Hr_Median_Diff.mat
Hr3_Median = Median_Diff;



figure(15); clf
semilogx(Periods,Hr1_Mean,'b--','LineWidth',2)
hold on
semilogx(Periods,Hr1_Median,'b','LineWidth',2)
semilogx(Periods,Hr3_Mean,'r--','LineWidth',2)
semilogx(Periods,Hr3_Median,'r-','LineWidth',2)
xlabel('Period (s)')
ylabel('dB Difference')
legend('Mean: 1 Hour Windows','Median: 1 Hour Windows','Mean: 3 Hour Windows','Median: 3 Hour Windows')


set(gca,'FontSize',20) 
xlim([2.5 500])