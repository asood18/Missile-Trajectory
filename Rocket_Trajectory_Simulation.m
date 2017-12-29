% Description:
% Predicts trajectories based on the given 
% rocket mass, engine thrust, launch parameters, and drag coefficient.

% Variable List:
% Delta = Time step (s)
% t = Time (s)
% Thrust = Thrust (kN)
% Mass = Mass (kton)
% Mass_Rocket_With_Motor = Mass with motor (kton)
% Mass_Rocket_Without_Motor = Mass without motor (kton)
% Theta = Angle (deg)
% C = Drag coefficient
% Rho = Air density (kg/m^3)
% A = Rocket projected area (m^2)
% Gravity = Gravity (m/s^2)
% Launch_Rod_Length = Length of launch rod (m)
% n = Counter
% Fn = Normal force (kN)
% Drag = Drag force (kN)
% Fx = Sum of forces in the horizontal direction (N)
% Fy = Sum of forces in the vertical direction (N)
% Vx = Velocity in the horizontal direction (m/s)
% Vy = Velocity in the vertical direction (m/s)
% Ax = Acceleration in the horizontal direction (m/s^2)
% Ay = Acceleration in the vertical direction (m/s^2)
% x = Horizontal position (m)
% y = Vertical position (m)
% Distance_x = Horizontal distance travelled (m)
% Distance_y = Vertical travelled (m)
% Distance = Total distance travelled (m)
% Memory_Allocation = Maximum number of time steps expected


%----------------------------------------------------Parameters----------------------------------------------------------------------%

clear all;       % Clear command window and workspace
close all;

% Parameters
Delta = 0.001;                  % Time step 
Memory_Allocation = 30000;      % Maximum number of time steps expected

% Preallocate memory for arrays
t = zeros(1, Memory_Allocation);
Thrust = zeros(1, Memory_Allocation);
Mass = zeros(1, Memory_Allocation);
Theta = zeros(1, Memory_Allocation);
Fn = zeros(1, Memory_Allocation);
Drag = zeros(1, Memory_Allocation);
Fx = zeros(1, Memory_Allocation);
Fy = zeros(1, Memory_Allocation);
Ax = zeros(1, Memory_Allocation);
Ay = zeros(1, Memory_Allocation);
Vx = zeros(1, Memory_Allocation);
Vy = zeros(1, Memory_Allocation);
x = zeros(1, Memory_Allocation);
y = zeros(1, Memory_Allocation);

Distance_x = zeros(1, Memory_Allocation);
Distance_y = zeros(1, Memory_Allocation);
Distance = zeros(1, Memory_Allocation);
%----------------------------------------------------------Constants--------------------------------------------------------------------%

C = input('Enter the value of Drag Coefficient : ');   % Drag coefficient=0.4
Rho = 1.2;                              % Air density (kg/m^3)
A = 4.9*10^-4;                          % Rocket projected area (Km^2)
Gravity = 9.81;                         % Gravity (m/s^2)
Launch_Rod_Length = 1;                  % Length of launch rod (Km)
Mass_Rocket_With_Motor =input('Enter the mass of rocket with motor : ');       % Mass with motor= 0.01546 k ton
Mass_Rocket_Without_Motor = input('Enter the mass of rocket without motor : ');     % Mass without motor=0.0117 k ton

Theta(1) =input('Enter the angle of Projection : ');% Initial angle = 60(deg)
Vx(1) = 0;                      % Initial horizontal speed (km/min)
Vy(1) = 0;                      % Initial vertical speed (km/min)
x(1) = 0;                       % Initial horizontal position (km)
y(1) = 0.1;                     % Initial vertical position (km)
Distance_x(1) = 0;              % Initial horizontal distance travelled (km)
Distance_y(1) = 0;              % Initial vertical distance travelled (km)
Distance(1) = 0;                % Initial  distance travelled (km)
Mass(1) = Mass_Rocket_With_Motor;       % Initial rocket mass (k ton)


%-------------------------------------------------------------------------------------------------------------------------------%
n = 1;                          % Initial time step            

while y(n) > 0                  % Run until rocket hits the ground
    n = n+1;                    % Increment time step
 
    t(n)= (n-1)*Delta;          % Elapsed time                     

   
    if t(n) <= 0.1                              % Launch phase 1
        Thrust(n) = 56*t(n);  
        Mass(n) = Mass_Rocket_With_Motor;
     elseif t(n) > 0.1 && t(n) < 0.5            % Launch phase 2
        Thrust(n) = 5.6;                          
        Mass(n) = Mass_Rocket_With_Motor;
    elseif t(n) >= 0.5 && t(n) < 3.5            % Launch phase 3
        Thrust(n) = 0;
        Mass(n) = Mass_Rocket_With_Motor;
    elseif t(n) >= 3.5                          % Launch phase 4                        
        Thrust(n) = 0;                                         
        Mass(n) = Mass_Rocket_Without_Motor;    % Rocket motor ejects
    end

    
    % Normal force calculations  
    if Distance(n-1) <= Launch_Rod_Length       % Launch rod normal force
        Fn(n) = Mass(n)*Gravity*cosd(Theta(1));
    else
        Fn(n) = 0;                              % No longer on launch rod
    end
    
    % Drag force calculation
    Drag(n)= 0.5*C*Rho*A*(Vx(n-1)^2+Vy(n-1)^2); % Calculate drag force
    
    % forces calculations Horiontal and Vertical 
    Fx(n)= Thrust(n)*cosd(Theta(n-1))-Drag(n)*cosd(Theta(n-1))-Fn(n)*sind(Theta(n-1));
                                   
    Fy(n)= Thrust(n)*sind(Theta(n-1))-(Mass(n)*Gravity)-Drag(n)*sind(Theta(n-1))+Fn(n)*cosd(Theta(n-1));
            
        
    % Acceleration calculations
    Ax(n)= Fx(n)/Mass(n);                       % Net acceleration in x direction 
    Ay(n)= Fy(n)/Mass(n);                       % Net acceleration in y direction
	
    % Velocity calculations
    Vx(n)= Vx(n-1)+Ax(n)*Delta;                 % Velocity in x direction
    Vy(n)= Vy(n-1)+Ay(n)*Delta;                 % Velocity in y direction
	
    % Position calculations
    x(n)= x(n-1)+Vx(n)*Delta;                   % Position in x direction
    y(n)= y(n-1)+Vy(n)*Delta;                   % Position in y direction
    
    % Distance calculations    
    Distance_x(n) = Distance_x(n-1)+abs(Vx(n)*Delta);      % Distance in x 
    Distance_y(n) = Distance_y(n-1)+abs(Vy(n)*Delta);      % Distance in y 
    Distance(n) = (Distance_x(n)^2+Distance_y(n)^2)^(1/2); % Total distance

    % Rocket angle calculation
    Theta(n)= atand(Vy(n)/Vx(n));      % Angle defined by velocity vector
end


%-------------------------------------------------------Plots----------------------------------------------------------------------%
figure('units','normalized','outerposition',[0 0 1 1]) 

% Figure 1
subplot(3,3,1)
plot(x(1:n),y(1:n)); 
xlim([0 inf]);
ylim([0 inf]);
xlabel({'Range (km)'});
ylabel({'Altitude (km)'});
title({'Trajectory'});

% Figure 2
subplot(3,3,2)
plot(t(1:n),Vx(1:n));
xlabel({'Time (min)'});
ylabel({'Vx (km/min)'});
title({'Horizontal Velocity'});

% Figure 3
subplot(3,3,3)
plot(t(1:n),Vy(1:n));
xlabel({'Time (min)'});
ylabel({'Vy (km/min)'});
title({'Vertical Velocity'});

% Figure 4
subplot(3,3,4)
plot(t(1:n),Theta(1:n));
xlabel({'Time (min)'});
ylabel({'Theta (Deg)'});
title({'Theta'});

% Figure 5
subplot(3,3,5)
plot(Distance(1:n),Theta(1:n));
xlim([0 2]);
ylim([59 61]);
xlabel({'Distance (km)'});
ylabel({'Theta (Deg)'});
title({'Theta at Launch'});

% Figure 6
subplot(3,3,6)
plot(t(1:n),Mass(1:n));
ylim([.0017 .02546]);
xlabel({'Time (min)'});
ylabel({'Mass (kton)'});
title({'Rocket Mass'});

% Figure 7
subplot(3,3,7)
plot(t(1:n),Thrust(1:n));
xlim([0 0.8]);
xlabel({'Time (min)'});
ylabel({'Thrust (kN)'});
title({'Thrust'});

% Figure 8
subplot(3,3,8)
plot(t(1:n),Drag(1:n));
xlabel({'Time (min)'});
ylabel({'Drag (kN)'});
title({'Drag Force'});

% Figure 9
subplot(3,3,9)
plot(Distance(1:n),Fn(1:n));
xlim([0 2]);
xlabel({'Distance (m)'});
ylabel({'Normal Force (kN)'});
title({'Normal Force'});
