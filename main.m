% Person properties
m_person = 77; %mass of person in chair (kg)

% Chair Properties
m_chair = 10; %mass of chair (kg)

% Environment properties
Temp_room = 300; %temperature of room (K)
Pressure_room = 101.325e3; %ambient pressure (Pa)
mu_s = 0.3; %static coefficient of friction
mu_k = 0.25; %kinetic coefficient of friction
g = 9.801; %gravitational constant (m/s/s)

%Gas properties
R = 287.058; %Specific gas constant for air (J/(kg*K))
T = Temp_room; %Temperature inside can

% Can properties
n = 350; %number of cans
m0 = 0.6; %initial mass of can (kg)
p0 = 1e6; %initial gage pressure (Pa)
V = pi*(0.028)^2*0.225; %volume of can (m^3)
A = pi*(1e-3)^2; %area of nozzle (m^2)

% Pack up parameters
params.g = g;
params.n = n;
params.R = R;
params.T = T;
params.V = V;
params.A = A;
params.p0 = p0;
params.m0 = m0;
params.mu_s = mu_s;
params.mu_k = mu_k;
params.m_person = m_person;
params.m_chair = m_chair;

% Integrate
tspan = [0,10];
x0 = [0;0];
[tout,xout] = ode45(@(t,x)odefun(t,x,params),tspan,x0);

figure
plot(tout,xout(:,1))
xlabel('time (sec)')
ylabel('distance (m)')


% N = 1e4;
% ns = linspace(1,1e5,N);
% maxdist = zeros(N,1);
% 
% for i = 1 : N
% params.n = ns(i);
% % Integrate
% tspan = [0,10];
% x0 = [0;0];
% [tout,xout] = ode45(@(t,x)odefun(t,x,params),tspan,x0);
% maxdist(i) = max(xout(:,1));
% end
% 
% close all
% figure
% subplot(1,2,1)
% plot(ns,maxdist)
% xlim([1,500])
% 
% subplot(1,2,2)
% plot(ns,maxdist)



function xdot = odefun(t,x,params)
%Unpack struct
n = params.n;
A = params.A;
m_person = params.m_person;
m_chair = params.m_chair;
g = params.g;

%Get force due to can thrust
p_can = get_can_pressure(t,params);
F_cans = n*A*p_can;


%Get total mass
m_cans = n * get_can_mass(p_can,params);
m_tot = m_person + m_chair + m_cans;

%Normal force
Fn = g * m_tot;

%Static friction force
Fs = Fn * params.mu_s;



%Kinetic friction force
Fk = Fn * params.mu_k;

%Acceleration
tol = 1e-10;
if F_cans - Fs <= 0 && x(2) < tol
    xdot = [0;0];
else
    a = 1/m_tot * (F_cans-Fk); 
    xdot = [ x(2) ; a ];
end
end

function p = get_can_pressure(t,params)
R = params.R;
V = params.V;
T = params.T;
A = params.A;
p0 = params.p0;

%Get current pressure
p = p0*exp( -A/V*sqrt(R*T) * t );
end

function m = get_can_mass(p,params)
R = params.R;
T = params.T;
m0 = params.m0;
p0 = params.p0;
V = params.V;

%Mass Loss
ml = (p0-p)*V^2/(R*T);

%Get current mass
m = m0 - ml;
end



