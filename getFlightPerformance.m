function [range, endDeviation, flightTime, maxHeight, landingAngle] = getFlightPerformance(flightData,bounce_beg,roll_beg)

%function [range, endDeviation, flightTime, maxHeight, landingAngle] = getFlightPerformance(flightData_SI)

load('ballParameters','dt');
range=flightData(bounce_beg,1);
endDeviation=flightData(end, 3);
% flightTime=dt*(length(flightData(1:roll_beg,1)) -1 + 5*(length(flightData((roll_beg+1):end,1))));
flightTime=dt*(length(flightData(1:bounce_beg,1)) -1);
maxHeight=max(flightData(1:bounce_beg,5));
landingAngle = -atan(flightData(bounce_beg-1,6)/...
    norm([flightData(bounce_beg-1,2) flightData(bounce_beg-1,4)]))*180/pi;
end