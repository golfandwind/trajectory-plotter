function wind = getWind(windInfo, height)

% function wind = getWind(windInfo, height)
load('ballParameters','z0');
windSpeed = windInfo(1); windHeading = windInfo(2)*pi/180; windElevation = windInfo(3)*pi/180;
wind = [windSpeed*cos(windElevation)*cos(windHeading), windSpeed*cos(windElevation)*sin(windHeading), windSpeed*sin(windElevation)];
if height<z0
    height = z0;
end

if windInfo(4)
    wind = wind*log(height/z0)/log(10/z0);
end
end