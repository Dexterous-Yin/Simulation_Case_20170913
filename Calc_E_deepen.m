function E = Calc_E_deepen(R, E_field, t, Scaling)

BE = 0.3106E5; % 0.275e5; % [nT] equatorial Earth B-field revised.

size_R=size(R,2);
E_background=zeros(3,size_R);

V_cor = E_field.V_cor;
V_con_A = E_field.A;

global L0_s sigmaL_s MLT0_s deltaB0_s sigMLT_s sigmaT t0
global L0_l sigmaL_l MLT0_l delMLT_l deltaB0_l factor_l deltaB0_l_c delMLT_l_c factor_l_c sigmal_l_c
global velocity
global dip_num MLT_interval

v_MLT = velocity/360*24;

for i_E = 1:size_R
    X = R(1,i_E);
    Y = R(2,i_E);
    Z = R(3,i_E);
    time=t(i_E)*Scaling.Time;
    
    r = sqrt(X.^2+Y.^2+Z.^2);
    Theta = pi/2-atan(Z./(sqrt(X.^2+Y.^2)));
    Phi = atan2(Y,X);
    s = cos(Theta);
    L = r./sin(Theta).^2;
    sinTheta=sin(Theta);
    cosTheta=cos(Theta);
    sinPhi=sin(Phi);
    cosPhi=cos(Phi);
    MLAT = atan(Z./(sqrt(X.^2+Y.^2)))/pi*180; %[degrees]
    MLT_s = (Phi+pi)/(2*pi)*24;
    MLT_l = (Phi+pi)/(2*pi)*24;
    MLTC_s = MLT0_s-v_MLT*time;
    MLTC_l = MLT0_l-v_MLT*time;

    % E_r.
    E_r = -V_con_A/sinTheta^2*sin(Phi+pi) - V_cor./r.^2*sinTheta.^2; %[mV/m]
    % E_Theta.
    E_Theta = 2*V_con_A*cosTheta/sinTheta^3*sin(Phi+pi)+V_cor./r.^2*2.*sinTheta.*cosTheta; %[mV/m]
    % E_Phi.
    E_Phi = -V_con_A/sinTheta^3*cos(Phi+pi);

    
    % transformation matrix.
    matrix_conv_sphere=[sinTheta*cosPhi, cosTheta*cosPhi, -1*sinPhi;...
        sinTheta*sinPhi, cosTheta*sinPhi, cosPhi;...
        cosTheta, -1*sinTheta, 0];

    E_tmp_VMSC = matrix_conv_sphere*[E_r;E_Theta;E_Phi];  
    
    t_dip = (24-MLTC_s)/v_MLT; %[s]   
    
    if abs(MLT_s-MLTC_s)>12
        if MLTC_s>=12
            MLT_s = MLT_s+24;
        else
            MLT_s = -(24-MLT_s);
        end
    end
    dip_E_s = Calc_deltaEr_deepen(s,L,MLT_s,deltaB0_s,L0_s,sigmaL_s,MLTC_s,sigMLT_s,v_MLT,t_dip,t0,sigmaT)*1e3;
    E_2 = dip_E_s;
    
    for j=0:dip_num-1
        MLT_l_j = MLT_l;
        if abs(MLT_l_j-MLTC_l-j*MLT_interval)>12
            if MLTC_l+j*MLT_interval>=12
                MLT_l_j = MLT_l_j+24;
            else
                MLT_l_j = -(24-MLT_l_j);
            end
        end
        if j==0
            dip_E_l = Calc_deltaEr(s,L,MLT_l_j,deltaB0_l_c,L0_l,sigmal_l_c,MLTC_l+j*MLT_interval,delMLT_l_c,v_MLT,factor_l_c)*1e3;
        else
            dip_E_l = Calc_deltaEr(s,L,MLT_l_j,deltaB0_l,L0_l,sigmaL_l,MLTC_l+j*MLT_interval,delMLT_l,v_MLT,factor_l)*1e3;
        end 
        E_2=E_2+dip_E_l; % mV/m
    end

    E_11 = 0; % mV/m
    E_3 = 0;
    
    matrix_conv = [-3*sinTheta*cosTheta*cosPhi/sqrt(1+3*cosTheta^2), (sinTheta^2-2*cosTheta^2)*cosPhi/sqrt(1+3*cosTheta^2), -sinPhi; ...
        -3*sinTheta*cosTheta*sinPhi/sqrt(1+3*cosTheta^2), (sinTheta^2-2*cosTheta^2)*sinPhi/sqrt(1+3*cosTheta^2), cosPhi; ...
        (sinTheta^2-2*cosTheta^2)/sqrt(1+3*cosTheta^2), 3*sinTheta*cosTheta/sqrt(1+3*cosTheta^2), 0];
    
    E_tmp = matrix_conv*[E_11;E_2;E_3];
    E_background(:,i_E) = E_tmp_VMSC+E_tmp;
end

E=E_background;


% Normalize E.
E = E/Scaling.E;
end


% h_phi
function res = h3(s,L)
    RE = 6.37e6;
    res = (1-s.^2).^(3/2).*L.*RE;
end

% the distribution of electric field
function res = Calc_deltaEr(s,L,MLT,deltaB0,L0,sigmaL,MLT0,delMLT,v_MLT,factor)
    RE = 6.37e6;
    res = pi/12.*deltaB0.*h3(s,L).*exp(-(L-L0).^2./(2.*sigmaL^2)).*v_MLT.*(tanh(factor*(MLT-MLT0+delMLT))+tanh(factor*(-MLT+MLT0+delMLT)))/(2.*tanh(factor*(delMLT)));
end

% the distribution of electric field
function res = Calc_deltaEr_deepen(s,L,MLT,deltaB0,L0,sigmaL,MLT0,sigMLT,v_MLT,t,t0,sigmaT)
    RE = 6.37e6;
    res_1 = pi/12.*deltaB0.*h3(s,L).*exp(-(L-L0).^2./(2.*sigmaL^2)).*(-(t-t0)/sigmaT^2).*exp(-(t-t0).^2./(2.*sigmaT.^2)).*(MLT-MLT0).*exp(-((MLT-MLT0).^2)./(2.*sigMLT.^2));
    res_2 = pi/12.*deltaB0.*h3(s,L).*exp(-(L-L0).^2./(2.*sigmaL^2)).*exp(-(t-t0).^2./(2.*sigmaT.^2)).*v_MLT.*(1-((MLT-MLT0).^2)./sigMLT.^2).*exp(-((MLT-MLT0).^2)./(2.*sigMLT.^2));   
    res=res_1+res_2;
end
