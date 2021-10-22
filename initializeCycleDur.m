function cycleDur = initializeCycleDur(Sm,Ssd)
%INITIALIZECYCLEDUR Returns Gaussian distributed duration of a cell cycle
%Also used to initialize death counter for immune cells
% Sm = mean duration, Ssd = standard deviation
meanCycle = Sm;
stdCycle  = Ssd;

meanCycleDur = meanCycle;
stdCycleDur = stdCycle;

cycleDur = -1;
while cycleDur<0
    cycleDur = floor(normrnd(meanCycleDur,stdCycleDur));
end
end