function [Time_new, R, Diagnostics] = Traj_3DAGCT_forward_deepen_time(Particle, Fields, Time, Numerics)
    %  Calculates particle's trajectory at times specified by Time[]
    %  array.
    % 
    % Returns:
    %    2D array R - positions of particle at requested times. 
    %      1st string R(1,:) = x, [RE] 
    %      2nd string R(2,:) = y, [RE]
    %      3rd string R(3,:) = z, [RE]
    %                 
    %    structure Diagnostics - all additional diagnostics can be added here
    %          vector W [MeV] - particle's kinetic energy at requested times
    %          vector p_par [m0*c] - particles's parallel momentum
    % 
    % Accepts: 
    %    Particle, Fields, Numerics stuctures and Time[] array, all w/ units:
    % 
    % structure Particle	describes the particle and its initial state
    %      Charge  [e]                - particle charge
    %      Mass     [el-n mass]   - particle mass
    %      W_ini     [MeV]           - initial particle energy
    %      Pitch_ini [deg]            - initial pitch angle
    %      column vector R         - initial particle position
    %      x   [RE]                      - x-coord
    %      y   [RE]                      - y-coord
    %      z   [RE]                      - z-coord
    % 
    % structure Fields	- describes e/m fields
    % structure B_field - B-field parameters	
    %      Field_Model  = 1; % 0 - analytical dipole, 1 - Tsyganenko
    %      year = 2000; day = 10; hour = 10; min = 10; sec = 0; dipole momentum vector parameters
    %      autotilt = 1; angle = 0; % 1 - the universal time is used to calculate dipole tilt, 0 - tilt angle is specified explicitly through 'angle'
    %      intB = 1; % 0 - pure dipole, 1 - IGRF
    %      extB = 1; % 0 - none, 1 - Tsyganenko T01_01
    %      Sol_P = 3.0; % [nPa] solar wind ram pressure
    %      DST = -20.0; % DST-index
    %      IMF_By = 3.0; % [nT] IMF By
    %      IMF_Bz = -5.0; % [nT] IMF Bz
    %      structure E_field - dummy (for future development)
    % 
    % vector Time [sec]	- a vector of times to return the trajectory at, 
    % 
    % structure Numerics
    % X_AbsTol  [R.E.]	- absolute particle position tolerance (in units of R.E.)
    % Integrator ['RK45']     - integrator to use (dummy)
    % 	
    %      Variables are scaled inside 'Trajectory' for actual calculations,
    %      but from the outside the fuction deals with non-normalized variables.
    %      Time[] can be either an array of times to return the trajectory at,
    %      or (if scalar) the number of periods of E-field to run the calculation for. 
    %      In the latter case, the trajectory points are returned at every period.
    %======================================================================
    % created by yzf 20200315
    

    % caluclate normalization coefficients
    Scaling = Calc_Scaling(Particle); 
    
    % normalize input parameters
    [Particle, Fields, Time] = Normalize(Particle, Fields, Time, Scaling);
    
    % calculate initial magnetic momentum (invariant)
    W0 = Particle.W_ini;
    %keyboard
    B0 = norm(Calc_B_deepen(Particle.R, Fields.B_field, Time(1), Scaling));
    Pitch0 = Particle.Pitch_ini;
    M0 = W0./B0.*(W0/2 + 1).*(sin(Pitch0)).^2 ; %magnetic momentum
    p_par0 = sqrt(W0.*(W0+2)).*cos(Pitch0); % initial parallel momentum

    X_ini = [Particle.R; p_par0; W0]; % inital condition 4-vector [x;y;z;p_par]
     
    % set particle position tolerances
     AbsTol(1) = Numerics.X_AbsTol; % abs.tol. in X
     AbsTol(2) = Numerics.X_AbsTol; % abs.tol. in Y
     AbsTol(3) = Numerics.X_AbsTol; % abs.tol. in Z
     AbsTol(4) = max([1e-12,abs(min(p_par0)).*Numerics.X_AbsTol]); % abs.tol. in p_par 
     AbsTol(5) = min(W0).*Numerics.X_AbsTol; % abs.tol. in W0

     options = odeset('RelTol',1e-22,'AbsTol',AbsTol,'OutputFcn',@odetpbar,'Events',@events);
    % calculate trajectory
    Charge = Particle.Charge;

     [t, X, te, Xe, ie] = ode45(@Calc_dX_dt, Time, X_ini, options, M0, Charge, Fields, Scaling); 
    
    % scale variables back to those with units
    X = X'; %transpose solution array
    % calculate kinetic energy array
    Diagnostics.W = Calc_W(X, M0, Fields.B_field, Time, Scaling)*Scaling.Energy; % calculate K.E. array
    % scale parallel momentum back to w/units 
    Diagnostics.p_par = X(4,:)*Scaling.p_par;
    Diagnostics.p_perp = Calc_Pperp(X, M0, Fields.B_field, Time, Scaling)*Scaling.p_par;
    % electric field along the trajectory
    Diagnostics.E = Calc_E_deepen(X, Fields.E_field, Time, Scaling)*Scaling.E;
    % magnetic field along the trajectory
    Diagnostics.B = Calc_B_deepen(X, Fields.B_field, Time, Scaling)*Scaling.B;
    % calculate the pitch angle of particles at each step
    Diagnostics.PA = Calc_PA(X, M0, Fields.B_field, Time, Scaling);

    % record the terminal contion
    Diagnostics.isloss = te;
    % scale position vectors back to w/units
    R(1,:) = X(1,:)*Scaling.r;
    R(2,:) = X(2,:)*Scaling.r;
    R(3,:) = X(3,:)*Scaling.r;
    
    Diagnostics.W1 = Calc_W1(X, M0, Fields.B_field, Scaling)*Scaling.Energy; % calculate K.E. array    

% check the boundary here.
    if ~isempty(te)
%         sprintf('At t = %1.2f seconds the concentration of A is %1.2f mol/L.', TE,VE);
        disp(['The boundary has been hit at time: ', num2str(te*Scaling.Time), ' s.']);
        disp(['The actual position (X, Y, Z) is: (', num2str(R(1,end)), ' ', num2str(R(2,end)), ' ', num2str(R(3,end)), ') RE.']);
%         disp(ie);
    end

    Time_new = t'*Scaling.Time;

function dX_dt = Calc_dX_dt(t, X, M, Charge, Fields, Scaling)
    % this is the function used by the ODE solver
    % calculates r.h.s. of the g.c. drift equations
    % X = [x; y; z; p_par], dX_dt = [dx/dt; dy/dt; dz/dt; dp_par/dt]
    % operates on normalized quantities
    R = X(1:3,1); % extract coordinates from 4-vector X
    p_par = X(4,1); % extract parallel momentum from 4-vector X
    % q = Charge;
    q = sign(Charge);
    dr = 1e-6; % half-range over which spatial central differences are calculated

    E = Calc_E_deepen(R, Fields.E_field, t, Scaling); % E-field vector
    B = Calc_B_deepen(R, Fields.B_field, t, Scaling); % B-field vector
    norm_B = norm(B); % magnitude of B
    % gamma = X(5,1)+1;
    gamma = sqrt(2*M*norm_B + p_par^2 + 1); % gamma
    dB_dr = Calc_dB_dr(R, Fields.B_field, Scaling, dr, t); % grad-B vector
    b = B/norm_B; % unit vector along field line
    B_2 = norm_B^2; % magnitude of B squared
    [db_ds, dB_ds] = Calc_dB_ds(R, Fields.B_field, Scaling, dr, t); % grad-B and vector grad-b along field line
 
    % calculate r.h.s. of the drift equations for dR/dt and dp_par/dt    
    dR_dt = cross(E,B)/B_2 + M/q/gamma*cross(B,dB_dr)/B_2 + p_par^2/q/gamma*cross(B,db_ds)/B_2 + b*p_par/gamma;
%{
Note: here we would like to add a term that is lacked in the parallel momentum equation. See the notes for detail.
Frank Wang.
%}
%     dp_dt = dot(E,b)*q - M/gamma*dB_ds;
    u_e=cross(E,B)/B_2;
%     dp_dt = dot(E,b)*q - M/gamma*dB_ds + p_par*(dot(u_e,db_ds));
    dp_dt = 0 - M/gamma*dB_ds + p_par*(dot(u_e,db_ds));
%     dgamma_dt
%     dW_dt
%     They are actually the same here.
%     dW_dt = q*sum(dR_dt.*E);
    % dW_dt = q*sum((dR_dt-b*p_par/gamma).*E);
%{
 dW_dt should include dB_dt related term. 
 Zefan Yin
%}
    dB_dt = Calc_dB_dt(R, Fields.B_field, Scaling, dr, t);
    dW_dt = q*sum(dR_dt.*E) + M/gamma*dB_dt;
%     disp(sign(sum(dR_dt.*E)));
%     disp(sum(dR_dt.*E));
    % build return vector dX_dt
%     if dB_dt ==0
%         disp(1);
%     end

    dX_dt = [dR_dt; dp_dt; dW_dt];
    
function dB_dr = Calc_dB_dr(R, B_field, Scaling, dr, t)
    % calculates gradients of B-field magnitude
    % operates on normalized quantities
    % using central differences

    dx = [dr;0;0];
    dB_dx = (norm(Calc_B_deepen(R + dx, B_field, t, Scaling)) - norm(Calc_B_deepen(R - dx, B_field, t, Scaling)))/2/dr;
    
    dy = [0;dr;0];
    dB_dy = (norm(Calc_B_deepen(R + dy, B_field, t, Scaling)) - norm(Calc_B_deepen(R - dy, B_field, t, Scaling)))/2/dr;

    dz = [0;0;dr];
    dB_dz = (norm(Calc_B_deepen(R + dz, B_field, t, Scaling)) - norm(Calc_B_deepen(R - dz, B_field, t, Scaling)))/2/dr;
    
    dB_dr = [dB_dx; dB_dy; dB_dz];
    
function dB_dt = Calc_dB_dt(R, B_field, Scaling, dr, t)
    % calculates dBdt
    % operates on normalized quantities
    % using central differences

    dt = 1e-4;
    dB_dt = (norm(Calc_B_deepen(R, B_field, t+dt, Scaling)) - norm(Calc_B_deepen(R, B_field, t-dt, Scaling)))/2/dt;
    
function [db_ds, dB_ds] = Calc_dB_ds(R, B_field, Scaling, dr, t)    
    % db_ds: calculates vector gradient of b along the field line
    % dB_ds: calculates gradient of B-magnitude along the field line
    % using central difference along the field line
    B = Calc_B_deepen(R, B_field, t, Scaling);
    b = B/norm(B);
    
    B_plus = Calc_B_deepen(R + dr*b, B_field, t, Scaling);
    norm_B_plus = norm(B_plus);
    b_plus = B_plus/norm_B_plus;
    
    B_minus = Calc_B_deepen(R - dr*b, B_field, t, Scaling);
    norm_B_minus = norm(B_minus);
    b_minus = B_minus/norm_B_minus;
        
    db_ds = (b_plus - b_minus)/2/dr;
    dB_ds = (norm_B_plus - norm_B_minus)/2/dr;  
    
function Scaling = Calc_Scaling(Particle)
    % calculates scaling coefficients for the variables
%     B0 = .275; % [G]  Earth's B-field at the surface (equator)
%     B0 = 0.3106; % [G]  Earth's B-field at the surface (equator)
    B0 = 0.3106; % [G]  Earth's B-field at the surface (equator)
    RE = 6.37e8; % [cm] - Earth's radius
%     RE = 6.37814e8; % [cm] - Earth's radius
    m0 = 9.1094e-28; % [g] - electron mass
    e  = 4.8032e-10; % [statcoul] - elementary charge
    c  = 2.99792e10; % [cm/sec] - speed of light
%     MeV = 1.6022e-6; % [erg] - energy [erg] equivalent to 1 MeV
    MeV = 1.6022e-6; % [erg] - energy [erg] equivalent to 1 MeV
    mV_m = 3.33564e-8; % [statvolt/cm] - E-field [statV/cm] equavalent to 1 mV/m
    Mass = Particle.Mass;
    Charge = abs(Particle.Charge);
    
    Scaling.B = Mass*m0*c^2/Charge/e/RE/B0; % [B.E.]
    Scaling.r = 1.0;    % [R.E.]
    Scaling.Energy = Mass*m0*c^2/MeV; % [MeV]
    Scaling.E = Mass*m0*c^2/Charge/e/RE/mV_m;   % [mV/m]
    Scaling.Frequency = c/RE/0.001; % [mHz]
    Scaling.Time = RE/c;    % [sec]
    Scaling.Angle = 180/pi; % [deg/rad]
    Scaling.p_par = 1; %[m0*c]

function [Particle_out, Fields_out, Time_out] = Normalize(Particle, Fields, Time, Scaling)
    % normalizes the quantities by dividing them by Scaling coefficients
    % Scaling - structure, containing scaling coefficients, generetaed by the 'Scaling fucntion'

    Particle.W_ini = Particle.W_ini/Scaling.Energy;
    Particle.Pitch_ini = Particle.Pitch_ini/Scaling.Angle;
    Particle.R = Particle.R/Scaling.r;

    Time = Time/Scaling.Time;
    
    Particle_out = Particle;
    Fields_out = Fields;
    Time_out = Time;
     
function W = Calc_W(X, M, B_field, Time, Scaling)
    % calculates kinetic energy 
    % operates on normalized quantities
    p_par = X(4,:);
    R = X(1:3,:);
    B = Calc_B_deepen(R, B_field, Time, Scaling);
    norm_B = sqrt(B(1,:).^2 + B(2,:).^2 + B(3,:).^2);
    W = sqrt(1 + 2.0*M*norm_B + p_par.^2) - 1.0;
    
    
function  Pperp = Calc_Pperp(X, M, B_field, Time, Scaling)
    B = Calc_B_deepen(X, B_field, Time, Scaling);
    norm_B = sqrt(B(1,:).^2 + B(2,:).^2 + B(3,:).^2);
    Pperp = sqrt(2.0*M*norm_B);
    
function PA = Calc_PA(X, M, B_field, Time, Scaling)
    % calculates the pitch angle of particles
    p_par = X(4,:);
    R = X(1:3,:);
    B = Calc_B_deepen(R, B_field, Time, Scaling);
    norm_B = sqrt(B(1,:).^2 + B(2,:).^2 + B(3,:).^2);
    PA = atan(sqrt(2.0*M*norm_B)./p_par);
    PA(find(PA<0)) = pi+PA(find(PA<0));
        
    
    
        
function W = Calc_W1(X, M, B_field, Scaling)
    % calculates kinetic energy 
    % operates on normalized quantities
    W = (X(5,:));
    
function [value,isterminal,direction] = events(t, y, M, Charge, Fields, Scaling)
% M, Charge, Fields, Scaling)
% Locate the boundary the particle crosses and stop integration.  
% R_bound = 10.5; %RE
% R_bound = 1+100e5/6.37e8; %RE
% R_bound = 5.5;
value = y(1)+5.5;
% value = y(2)+3;
isterminal = 1;   % stop the integration
direction = 0;   % Detect all zero crossings
% direction = -1;   % negative direction
% direction = 1;   % positive direction