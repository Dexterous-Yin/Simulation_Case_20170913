%{
This code is used to reproduce the particle observation from RBSP during 20170913 event.
Zefan
2021-03-18
%}

% load the satellite's location
filename = 'data\rbspb_def_MagEphem_TS04D_20170912_v1.0.0.h5';
fileinfo = h5info(filename);
rgsm = h5read(filename,'/Rsm'); rgsm = rgsm(:,1:end-1);
utc = h5read(filename,'/UTC'); utc = utc(1:end-1)'; % hours from 0:00 of the current day
Lm = h5read(filename,'/L'); Lm = Lm(1,1:end-1); % Mcllwain L-shell value.
Lstar = h5read(filename,'/Lstar'); Lstar = Lstar(1,1:end-1); % L*
xgsm = rgsm(1,:);
ygsm = rgsm(2,:);
zgsm = rgsm(3,:);

filename_add = 'data\rbspb_def_MagEphem_TS04D_20170913_v1.0.0.h5';
fileinfo_add = h5info(filename_add);
rgsm_add = h5read(filename_add,'/Rsm'); rgsm_add = rgsm_add(:,1:end-1);
utc_add = h5read(filename_add,'/UTC'); utc_add = utc_add(1:end-1)'; % hours from 0:00 of the current day
Lm_add = h5read(filename_add,'/L'); Lm_add = Lm_add(1,1:end-1); % Mcllwain L-shell value.
Lstar_add = h5read(filename_add,'/Lstar'); Lstar_add = Lstar_add(1,1:end-1); % L*
xgsm_add = rgsm_add(1,:);
ygsm_add = rgsm_add(2,:);
zgsm_add = rgsm_add(3,:);

utc=[utc,24+utc_add];
xgsm = [xgsm,xgsm_add];
ygsm = [ygsm,ygsm_add];
zgsm = [zgsm,zgsm_add];
Lm = [Lm,Lm_add];
Lstar = [Lstar,Lstar_add];

utc_sm = utc(1):1/60:utc(end);
meettime = 24+40/60; % 07-13/00:40:00
meeti = find(abs(utc_sm-meettime)<0.0001);

r = sqrt(xgsm.^2+ygsm.^2+zgsm.^2);
Theta = pi/2-atan(zgsm./(sqrt(xgsm.^2+ygsm.^2)));
Phi_sat = atan2(ygsm,xgsm);
L_sat = r./sin(Theta).^2;
MLAT_sat = atan(zgsm./(sqrt(xgsm.^2+ygsm.^2)))/pi*180; %[degrees]
MLT_sat = (Phi_sat+pi)./(2*pi)*24;

L_sat_sm = interp1(utc,Lm,utc_sm,'pchip'); L_sat_sm = smoothdata(L_sat_sm,'movmean',10);
MLT_sat_sm = interp1(utc,MLT_sat,utc_sm);

temp = find(utc_sm >=23 & utc_sm <=24+4);
L_sat_focus = L_sat_sm(temp);
MLT_sat_focus = MLT_sat_sm(temp);
time_sat_focus = utc_sm(temp); %[hr]

% The time when satellite observes the deepest magnetic dip
meetid = find(abs(time_sat_focus-meettime)<1e-5);
meetMLT = MLT_sat_focus(meetid); % 12.7419
meetL = L_sat_focus(meetid); % Ldipole-5.75 Lm-5.4695 L*-4.4496

MLT_ob = MLT_sat_focus; % satellite MLT
L_ob = L_sat_focus; % satellite L
time_ob = time_sat_focus;

E_ob = [20.6866,24.1037,28.0853,32.7245,38.1301,44.4287,51.7677]; %[keV]

flux_ob = zeros(length(E_ob),length(time_ob));

t_rec = cell(length(E_ob),length(L_ob));
x_rec = cell(length(E_ob),length(L_ob));
y_rec = cell(length(E_ob),length(L_ob));
MLT_rec = cell(length(E_ob),length(L_ob));
L_rec = cell(length(E_ob),length(L_ob));
W_rec = cell(length(E_ob),length(L_ob));
PA_rec = cell(length(E_ob),length(L_ob));
B_rec = cell(length(E_ob),length(L_ob));
time_set = zeros(length(E_ob),length(L_ob));
% set the distributation factors of magnetic dip
% basic settings
BE = 0.3106*10^5; % [nT] equatorial magnetic field at L=1
RE = 6.37*1e6; 

%%
% model settings
global L0_s sigmaL_s MLT0_s deltaB0_s sigMLT_s sigmaT t0
global L0_l sigmaL_l MLT0_l delMLT_l deltaB0_l factor_l deltaB0_l_c delMLT_l_c factor_l_c sigmal_l_c
global velocity
global dip_num MLT_interval
dip_num=18;
MLT_interval=0.5;

% set the time zero at midnight and reaches the deepest situation at 12MLT
% NOTE the sigmaT should corresponds with the deltaB0_s
sigmaT = 2*60*60; %[s]
maxMLT=12;
v_MLT_hr = 2; %[MLT/hr]
bound_MLT = 0;
Field_A = 0.9; % equatorial electric field strength % [mV/m]

t0 = (24-maxMLT)/v_MLT_hr*60*60; %[s]
velocity = v_MLT_hr/24*360/60/60; %[degree/s]

% if want to reproduce the energy spectrum of proton fluxes, use E_id=[3,4,5,6] and T_id=31:300
for E_id = [4,6] %[3,4,5,6] 
    for T_id = 102 %31:300 

        % set the parameters for the smaller hole
        L0=5.5;
        L0_s = L0;
        sigmaL_s = 0.2;
        deltaB0_s = -50*1e-9; % [T] % set this as zeros for bouncing protons
        MLT0_s = meetMLT-(time_ob(T_id)-meettime)*v_MLT_hr; % the location of smaller dip at the observation time.
        sigMLT_s = 0.3;

        % set the parameters for the larger hold
        L0_l = L0;
        sigmaL_l = 0.4;
        sigmal_l_c = 0.4;
        deltaB0_l = -30*1e-9;  % [T] % set this as zeros for bouncing protons
        deltaB0_l_c = -50*1e-9; % [T] % set this as zeros for bouncing protons
        MLT0_l = MLT0_s+0.5;
        delMLT_l = 0.2;
        factor_l = 3;
        delMLT_l_c = 0.2;
        factor_l_c = 3;
        
        no_particles = 1;
        
        simTime = 25000; %(24-MLT0_s)/2*60*60; %(T_test-15)/2*60*60; %10000; %[s]
        MLT_launch = MLT_ob(T_id); 
        phi_launch = MLT_launch/24*2*pi-pi;
        lambda_launch = 0; % launch at equator
        
        L = L_ob(T_id);        
        x_start=L*cos(lambda_launch)^3*cos(phi_launch);
        y_start=L*cos(lambda_launch)^3*sin(phi_launch);
        z_start=L*cos(lambda_launch)^2*sin(lambda_launch);
        
        for i = 1:no_particles
            Particle(i).Charge = 1.0; %[elementary charge]
            Particle(i).Mass = 1.6726E-27/9.1094E-31; %[el-n mass], here for proton.
            Particle(i).W_ini = E_ob(E_id)*10^-3; %[MeV]
            Particle(i).Pitch_ini = 90;  %initial pitch angle [deg] % for perpendicular-moving protons
            % Particle(i).Pitch_ini = 27; %initial pitch angle [deg] % for bouncing protons
            Particle(i).R(1,1) = x_start; %x - [R.E.]
            Particle(i).R(2,1) = y_start; %y - [R.E.]
            Particle(i).R(3,1) = z_start; %y - [R.E.]
        end
        
        Fields.B_field.Field_Model  = 0; % 0 - analytical dipole
        
        % The background E field we shall first use simple convection plus
        % corotation electric field.
        % corotation E.
        Fields.E_field.V_cor = 14.4427; % unit: mV/m.
        % convection E.
        Fields.E_field.A = Field_A; % -->[mV/m]        
        
        % set absolute particle position accuracy as X_AbsTol
        % in units of R.E.
        Numerics.X_AbsTol = 1.e-6; %[R.E.]
        
        Time = 0:-0.5:-simTime; % output particle position at specified moments of time
        
        % solve trajectory
        for i = 1:no_particles
            [Time_new, R, Diagnostics] = Traj_3DAGCT_forward_deepen_time(Particle(i), Fields, Time, Numerics);
            RR(i,:,:) = R;
        end
        
        if isempty(Diagnostics.isloss) 
            clear RR;
            continue;
        end
        
        Time_old = Time;
        Time = Time_new;
        
        X = squeeze(RR(:,1,:));
        Y = squeeze(RR(:,2,:));
        Z = squeeze(RR(:,3,:));
        
        Lambda = atan2(Z,sqrt(X.^2+Y.^2));
        Phi = atan2(Y,X);
        Theta = pi/2-atan(Z./(sqrt(X.^2+Y.^2)));
        sinTheta=sin(Theta);
        cosTheta=cos(Theta);
        sinPhi=sin(Phi);
        cosPhi=cos(Phi);
        MLT = (Phi+pi)./(2*pi).*24;
        
        tempMLT = find(MLT>6 & MLT<12);
        if ~isempty(tempMLT)
            clear RR;
            continue;
        end

        if Y(end)<0
            clear RR;
            continue;
        end
        
        R = sqrt(X.^2+Y.^2+Z.^2);
        R0 = R./cos(Lambda).^2;
        Height = (R-1)*6.37e3;
        
        %%%
        % Diagnose the pitch angle and loss cone
        PitchA = Diagnostics.PA; % annotate this if for 90 degree particles
        KEnergy = Diagnostics.W;
        
        t_rec{E_id,T_id} = Time;
        x_rec{E_id,T_id} = X;
        y_rec{E_id,T_id} = Y;
        L_rec{E_id,T_id} = R0;
        MLT_rec{E_id,T_id} = MLT;
        W_rec{E_id,T_id} = KEnergy.*1e3;
        PA_rec{E_id,T_id} = PitchA./pi*180;
        B_rec{E_id,T_id} = (Diagnostics.B(1,:)'.^2+Diagnostics.B(2,:)'.^2+Diagnostics.B(3,:)'.^2).^0.5.*BE; %[nT]
        time_set(E_id,T_id) = Time(end);
        
        flux_ob(E_id,T_id) = Time(end)+time_ob(T_id)*60*60; % Time from 09-12/00:00
        
        disp([num2str(E_id),'  ',num2str(T_id),'  ',datestr(datenum('2017-09-12/00:00:00')+flux_ob(E_id,T_id)/60/60/24)]);
        
        clear RR;
    end
end
flux_ob_rec = flux_ob;

%%
% set rbsp location
formatIn = 'yyyy-mm-dd/HH:MM:SS';
start_time = datenum('2017-09-12/22:00:00',formatIn);
end_time = datenum('2017-09-13/04:00:00',formatIn);
time_b = utc_sm./24+datenum('2017-09-12/00:00:00',formatIn);
pick_b = find(time_b>=start_time & time_b<=end_time);
time_b_plot = time_b(pick_b);
Phi_plot = (MLT_sat_sm(pick_b)-12)./12.*pi;
L_sat_plot = L_sat_sm(pick_b);
x_b_eq = L_sat_plot.*cos(Phi_plot);
y_b_eq = L_sat_plot.*sin(Phi_plot);

figure;
cmap = colormap('jet');
% plot the magnetic dip and the particle trajectories
% the magnetic dip structure
xmin = -8; xmax = 8;
ymin = -8; ymax = 8;
x_plot = [xmin:0.05:xmax]';
y_plot = [ymin:0.05:ymax];
L_plot = sqrt(x_plot.^2+y_plot.^2);
Phi_plot = atan2(y_plot,x_plot);
MLT_plot_s_rec = (Phi_plot+pi)/(2*pi)*24;
MLT_plot_l_rec = (Phi_plot+pi)/(2*pi)*24;

v_MLT = velocity/360*24;

t_plot = [0,-2000,-4000];

T_id = 102;
MLT0_s_p = meetMLT-(time_ob(T_id)-meettime)*v_MLT_hr;
for t= t_plot
    MLT_plot_s = MLT_plot_s_rec;
    MLT_plot_l = MLT_plot_s_rec;
    
    MLTC_s = MLT0_s_p-v_MLT*t;
    MLT_plot_s(abs(MLT_plot_s-MLTC_s)>12) = MLT_plot_s(abs(MLT_plot_s-MLTC_s)>12)+24;
    
    t_dip=(24-MLTC_s)/v_MLT;
    dip_B_s = Calc_deltaB_deepen(L_plot,MLT_plot_s,deltaB0_s,L0_s,sigmaL_s,MLTC_s,sigMLT_s,t_dip,t0,sigmaT).*1e9; %[nT]
    dip_B = dip_B_s;
    
    for j=0:dip_num-1
        if j == 0
            MLTC_l = MLTC_s+0.5+j*MLT_interval;
            MLT_plot_l(abs(MLT_plot_l-MLTC_l)>12) = MLT_plot_l(abs(MLT_plot_l-MLTC_l)>12)+24;
            dip_B_l = Calc_deltaB(L_plot,MLT_plot_l,deltaB0_l_c,L0_l,sigmal_l_c,MLTC_l,delMLT_l_c,factor_l_c).*1e9;
        else
            MLTC_l = MLTC_s+0.5+j*MLT_interval;
            MLT_plot_l(abs(MLT_plot_l-MLTC_l)>12) = MLT_plot_l(abs(MLT_plot_l-MLTC_l)>12)+24;
            dip_B_l = Calc_deltaB(L_plot,MLT_plot_l,deltaB0_l,L0_l,sigmaL_l,MLTC_l,delMLT_l,factor_l).*1e9;
        end
        dip_B = dip_B+dip_B_l;
    end
    
    figure('Color',[1 1 1]);
    pco = pcolor(y_plot',x_plot,dip_B);
    set(pco,'LineStyle','none');
    colormap bone;
    c = colorbar('southoutside');
    c.Label.String = 'dB (nT)';
    set(gca,'XDir','reverse');
    set(gca,'XAxisLocation','origin'); set(gca,'YAxisLocation','origin');
    set(gca,'Layer','top');
    caxis([-60,20]);
    xlabel({'Y';'Dawn'}); ylabel({'X';'Noon'});
    hold on;
    
    % plot the rbsp-trajectories
    plottime_b = time_ob(T_id)./24+datenum('2017-09-12/00:00:00',formatIn)+t/60/60/24;
    marker_b_x = interp1(time_b_plot,x_b_eq,plottime_b);
    marker_b_y = interp1(time_b_plot,y_b_eq,plottime_b);
    plot(y_b_eq,x_b_eq,'--','Color','c','LineWidth',1); hold on;
    plot(marker_b_y,marker_b_x,'x','MarkerEdgeColor','c','MarkerSize',8,'LineWidth',1); hold on;
    
    % plot the boundary
    bound_x = [-5.5,-5.5,-6,-6];
    bound_y = [0,4,4,0];
    patch(bound_y,bound_x,[0.9290 0.6940 0.1250],'EdgeColor','none'); hold on; % source region

    % plot the particle trajectories
    for E_id=[4,6]
        for T_id= 102
            t_stop = find(abs(t_rec{E_id,T_id})>abs(t));
            if ~isempty(t_stop)
                t_stop = t_stop(1);
            end
            if E_id==4
                colorE='b';
            else
                colorE='r';
            end
            % plot(y_rec{E_id,T_id}(1:t_stop),x_rec{E_id,T_id}(1:t_stop),'Color',colorE,'LineWidth',1.5); hold on;
            plot(y_rec{E_id,T_id}(t_stop:end),x_rec{E_id,T_id}(t_stop:end),'Color',colorE,'LineWidth',1.5); hold on;
            % plot(y_rec{E_id,T_id}(950:2150),x_rec{E_id,T_id}(950:2150),'Color',colorE,'LineWidth',1.5); hold on;
            plot(y_rec{E_id,T_id}(end),x_rec{E_id,T_id}(end),'x','MarkerEdgeColor',colorE,'MarkerSize',8,'LineWidth',1); hold on;
        end
    end
        %
    axis equal;
    ylim([xmin,xmax]); xlim([ymin,ymax]);
    % title(['E=',num2str(E_ob(E_id)),'keV',' T=',num2str(t),'s',' MLT=',num2str(MLTC_s,'%2.2f')]);
    title(['T=',num2str(t),'s']);
    tempp = get(gca,'position');
    % set(gca,'position',[tempp(1),tempp(2),1,0.2]);
    set(gca,'Box','on');
    set(gca,'XTick',[-8:2:8]);
    % print(['simulation_result\trajectories_',num2str(find(t_plot==t))],'-dpng','-r1000'); % demo output
    % close;         
end


%%%%%%%%%%%%%%%%%%%%%%functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the distribution of magnetic dip
function res = Calc_deltaB(L,MLT,deltaB0,L0,sigmaL,MLT0,delMLT,factor)
    res = deltaB0.*exp(-(L-L0).^2./(2.*sigmaL^2)).*((tanh(factor*(MLT-MLT0+delMLT))+tanh(factor*(-MLT+MLT0+delMLT)))/(2.*tanh(factor*delMLT)));
end

% the distribution of magnetic dip
function res = Calc_deltaB_deepen(L,MLT,deltaB0,L0,sigmaL,MLT0,sigMLT,t,t0,sigmaT)
    res = deltaB0.*exp(-(t-t0).^2./(2.*sigmaT.^2)).*exp(-(L-L0).^2./(2.*sigmaL^2)).*(1-((MLT-MLT0).^2)./sigMLT.^2).*exp(-((MLT-MLT0).^2)./(2.*sigMLT.^2));
end