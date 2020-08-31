function spinRate = getSpinRate(launchSpinRate, time)

%function spinRate = getSpinRate(launchSpinRate_RPM, time_sec)

load('ballParameters','T');
spinRate=launchSpinRate*exp(-time/T);
end