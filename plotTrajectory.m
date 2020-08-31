function [landingPoints, flightPerformance, flightData, flag] = plotTrajectory(launchSpeed, launchSpinRate, launchHeading, loft, spinAxis, windInfo)
%function [landingPoints=[first_x, first_y, last_x, last_y], flightPerformance=[range, endDeviation, flightTime, maxHeight, landingAngle], flightData=[X, X_dot, Y, Y_dot, Z, Z_dot]] = plotTrajectory(launchSpeed_m/s, launchSpinRate_RPM, launchHeading_degrees, loft_degrees, spinAxis_degrees, [windSpeed_m/s, windHeading_degrees, windElevation_degrees, windModel_const(0)/log(1))

%% Initialization
%tic
load('ballParameters','g','rho','mass','dt','radius','CDmodel','CLmodel','mu','mu_roll')
plotMode = 1;   % 0 (none), 1 (plot)
flag = 0;   % 0 (normal), 1 (killed) [IGNORE]

loft = loft*pi/180;
launchHeading = launchHeading*pi/180;
spinAxis = spinAxis*pi/180;
spin_cap = [cos(spinAxis)*sin(launchHeading), -cos(spinAxis)*cos(launchHeading), sin(spinAxis)];

%% Computing the trajectory
i = 1;beg = 1;j = 1;beg_g = 1;
Y = [0, launchSpeed*cos(loft)*cos(launchHeading), 0, launchSpeed*cos(loft)*sin(launchHeading), 0, launchSpeed*sin(loft)];
localHeight = 0; ballHeight = 0;

while 1
    % flight
    while ((i<=beg+1)||((Y(i,5)>=localHeight)))
        wind = getWind(windInfo, Y(i,5));
        airSpeed = [Y(i,2)-wind(1), Y(i,4)-wind(2), Y(i,6)-wind(3)];
        L = 0.5*rho*norm(cross(spin_cap,airSpeed))*cross(spin_cap,airSpeed)*pi*(radius^2)*CLmodel([norm(airSpeed),getSpinRate(launchSpinRate,(i-1)*dt)]);
        D = 0.5*rho*norm(airSpeed)*(-airSpeed)*pi*(radius^2)*CDmodel([norm(airSpeed),getSpinRate(launchSpinRate,(i-1)*dt)]);
        Y_dot = [Y(i,2), (L(1)+D(1))/mass, Y(i,4), (L(2)+D(2))/mass, Y(i,6), (L(3)+D(3)-mass*g)/mass];
        Y(i+1,:) = Y(i,:) + dt*Y_dot;
        i = i + 1;
        ballHeight(i) = Y(i,5) - localHeight;
    end

    if Y(i,5)<=localHeight
        Y(i,5) = localHeight;
        h_max = max(ballHeight(beg:end));
        beg = i;
        if j==1
            beg_g = beg;
        end
        if h_max<0.005
            break
        end
        % bounce
        % break
        v = Y(i,[2,4,6]);
        y_cap = [0,0,1]; y_cap = y_cap/norm(y_cap);
        z_cap = cross(v,y_cap); z_cap = z_cap/norm(z_cap);
        x_cap = cross(y_cap,z_cap);
        v_ix = dot(v,x_cap); v_iy = dot(v,y_cap);
        w_ix = (2*pi/60)*getSpinRate(launchSpinRate,(i-1)*dt)*dot(spin_cap,x_cap);
        w_y = (2*pi/60)*getSpinRate(launchSpinRate,(i-1)*dt)*dot(spin_cap,y_cap);
        w_iz = (2*pi/60)*getSpinRate(launchSpinRate,(i-1)*dt)*dot(spin_cap,z_cap);
        theta_c = deg2rad(15.4)*(norm(v)/18.6)*(atan2(v_ix,-v_iy)/deg2rad(44.4));
        v_ixp = v_ix*cos(theta_c) + v_iy*sin(theta_c);
        v_iyp = v_iy*cos(theta_c) - v_ix*sin(theta_c);
        e_z = (v_iyp>=-20)*(0.51 + 0.0375*v_iyp + 0.000903*(v_iyp^2)) ...
            + (v_iyp<-20)*0.12;
        e_x = (v_iy>=-20)*(0.51 + 0.0375*v_iy + 0.000903*(v_iy^2)) ...
            + (v_iy<-20)*0.12;
        mu_cz = -2*(v_ixp + radius*w_iz)/(7*v_iyp*(1 + e_z));
        mu_cx = 2*radius*w_ix/(7*v_iy*(1 + e_x));
        
        if mu < mu_cz
            v_rxp = v_ixp + mu*v_iyp*(1 + e_z);
            v_ryp = -e_z*v_iyp;
            w_rz = w_iz + 5*mu*v_iyp*(1 + e_z)/(2*radius);
        else
            v_rxp = 5*v_ixp/7 - 2*radius*w_iz/7;
            v_ryp = -e_x*v_iyp;
            w_rz = -v_rxp/radius;
        end
        if mu < mu_cx
            v_rz = mu*v_iy*(1 + e_x);
            w_rx = w_ix - 5*mu*v_iy*(1 + e_x)/(2*radius);
        else
            v_rz = 2*radius*w_ix/7;
            w_rx = v_rz/radius;
        end
        
        v_rx = v_rxp*cos(theta_c) - v_ryp*sin(theta_c);
        v_ry = v_rxp*sin(theta_c) + v_ryp*cos(theta_c);
        Y(i,[2,4,6]) = v_rx*x_cap + v_ry*y_cap + v_rz*z_cap;
        w = w_rx*x_cap + w_y*y_cap + w_rz*z_cap;
        launchSpinRate = (30/pi)*norm(w); spin_cap = w/norm(w);
    else
        flag = 1;
        break
    end
    j = j + 1;
end

% rolling
if ~flag
    v = Y(i,[2,4,6]); v_dot = 0;
    while (i==beg) || (norm(v)>0.5) || (norm(v_dot)>1)
        % break
        gx = 0; gy = 0;
        if norm(v)
            y_cap = [-gx,-gy,1]; y_cap = y_cap/norm(y_cap);
            z_cap = cross(v,y_cap); z_cap = z_cap/norm(z_cap);
            x_cap = cross(y_cap,z_cap);
        end
        t_cap = [gx,gy,(gx^2 + gy^2)];
        if norm(t_cap)
            t_cap = t_cap/norm(t_cap);
        end
        if norm(dot(v,x_cap))<0.5
            break;
        end
        if i==beg
            v = dot(v,x_cap)*x_cap;
        else
            v = norm(v)*x_cap;
        end
        
        g_t = dot([0,0,-g],t_cap)*t_cap;
        f = -(5/7)*mu_roll*abs(dot([0,0,-g],y_cap))*x_cap*(norm(v)~=0);
        v_dot = g_t + f;
        
        Y_dot = [v(1),v_dot(1),v(2),v_dot(2),v(3),v_dot(3)];
        Y(i+1,:) = Y(i,:) + 5*dt*Y_dot;
        i = i + 1;
        Y(i,5) = 0;
        v = Y(i,[2,4,6]);
    end
end

if beg==1
    beg = i;beg_g = i;
end

flightData = Y;
[range, endDeviation, flightTime, maxHeight, landingAngle] = getFlightPerformance(flightData,beg_g,beg);
flightPerformance = [range, endDeviation, flightTime, maxHeight, landingAngle];
landingPoints = [Y(beg_g,1), Y(beg_g,3), Y(end,1), Y(end,3)];
%toc

%% Plotting the trajectory
if plotMode
    if plotMode==1
        if flag==1  % killed [IGNORE]
            plot3(Y(:,1),Y(:,3),Y(:,5),'r','LineWidth',2)
        else
            % normal
            plot3(Y(:,1),Y(:,3),Y(:,5),'k','LineWidth',2)
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
        end
        daspect([1 1 1]);grid on;hold on
    end
end

end
