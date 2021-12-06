function B = Calc_B_deepen(R, B_field, t, Scaling)
    % calculates B-field at positions specified by array R
    % B = [Bx; By; Bz] is a 2D array where first row is Bx-s, second row is
    % By-s, third row is Bz-s.
    % X = [x; y; z] is a 2D array where first row is x-s, second row is
    % y-s, third row is z-s.
    % This function operates on normalized
    % quantities, i.e., x,y,z are in [R.E.], B in [m0*c^2/q/RE], t is in [R.E./c]
    % created from Calc_B_2.m by yzf 20200315 
    Field_Model  = B_field.Field_Model;
%     YEAR =         B_field.year;
%     DAY =          B_field.day;
%     HOUR =         B_field.hour;
%     MIN =          B_field.min;
%     SEC =          B_field.sec;
%     autotilt =     B_field.autotilt;
%     angle =        B_field.angle;
%     intB =         B_field.intB;
%     extB =         B_field.extB;
%     Sol_P =        B_field.Sol_P;
%     DST =          B_field.DST;
%     IMF_By =       B_field.IMF_By;
%     IMF_Bz =       B_field.IMF_Bz;

    % BE = 0.30115e5; % 0.275e5; % [nT] equatorial Earth B-field
    BE = 0.3106E5; % 0.275e5; % [nT] equatorial Earth B-field revised.

    x = R(1,:);
    y = R(2,:);
    z = R(3,:);

    %     R_size = size(x,2);

    if Field_Model == 0 % calculate B as analytical dipole
        r = sqrt(x.^2 + y.^2 + z.^2);
        r_3 = r.^3;
        r_5 = r.^5;

        % calculate field components in [B.E.]
        B_x = - 3.*x.*z./r_5;
        B_y = - 3.*y.*z./r_5;
        B_z = 1./r_3 - 3*z.*z./r_5;

        % build B-vector and normalize it
        %         B = [B_x; B_y; B_z]/Scaling.B;
        B = [B_x; B_y; B_z];

        % else % calculate B using FORTRAN compiled routines (Tsyganenko)
    else % calculate B with geopack_Matlab. Or even we can use IRBEM-4.4.0 if needed.

        R_size = size(x,2);  %number of (x,y,z) points in R-array - calculate B for each of them
        PARMOD = [Sol_P, DST, IMF_By, IMF_Bz]; % parameters for external field

        %         Choose internal magnetic field model.
        switch intB
            case 0
                INNAME='GEOPACK_DIP';
            case 1
                INNAME='GEOPACK_IGRF_GSM';
            otherwise
                disp('Undefined choice of intB!');
                keyboard;
        end

        %         Choose external magnetic field model.
        switch extB
            case 0
                EXNAME='T89';
            case 1
                EXNAME='T96';
            case 2
                EXNAME='T01';
            otherwise
                disp('Undefined choice of extB!');
                keyboard;
        end

        for i = 1:R_size
            X_gsm = x(i); Y_gsm = y(i); Z_gsm = z(i);
            [B_gsm_x,B_gsm_y,B_gsm_z]=tsycalc(X_gsm,Y_gsm,Z_gsm,INNAME,EXNAME,PARMOD,YEAR,DAY,HOUR,MIN,SEC);
            %             [B_gsm_x, B_gsm_y, B_gsm_z] = B_GSM(X_gsm, Y_gsm, Z_gsm, intB, extB, parmod, geopack1, geopack2);
            B_x(i) = B_gsm_x; B_y(i) = B_gsm_y; B_z(i) = B_gsm_z;
        end

        % build B-vector and normalize it
        B = [B_x; B_y; B_z]/BE; % B-field in units of [BE]
    end

    % set the magnetic dip here
    
    size_R=size(R,2);
    dip_B = zeros(3,size_R);
    global L0_s sigmaL_s MLT0_s deltaB0_s sigMLT_s sigmaT t0
    global L0_l sigmaL_l MLT0_l delMLT_l deltaB0_l factor_l deltaB0_l_c delMLT_l_c factor_l_c sigmal_l_c
    global velocity  
    global dip_num MLT_interval
 
    v_MLT = velocity/360*24;

    for i=1:size_R
        X = R(1,i);
        Y = R(2,i);
        Z = R(3,i);
        time=t(i)*Scaling.Time; % [s]
        if size(X,1)>1
            error('Length is more than one.');
        end
        r = sqrt(X.^2+Y.^2+Z.^2);
        Theta = pi/2-atan(Z./(sqrt(X.^2+Y.^2)));
        Phi = atan2(Y,X);
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
        
        t_dip = (24-MLTC_s)/v_MLT; %[s] 
               
        if abs(MLT_s-MLTC_s)>12 
            if MLTC_s>=12
                MLT_s = MLT_s+24;
            else
                MLT_s = -(24-MLT_s);
            end      
        end
        
        dip_B_s = Calc_deltaB_deepen(L,MLT_s,deltaB0_s,L0_s,sigmaL_s,MLTC_s,sigMLT_s,t_dip,t0,sigmaT).*1e9; %[nT]
        B11 = dip_B_s; % [nT]
        
        for j=0:dip_num-1
            MLT_l_j = MLT_l; % particle MLT
            if abs(MLT_l_j-MLTC_l-j*MLT_interval)>12
                if MLTC_l+j*MLT_interval>=12
                    MLT_l_j = MLT_l_j+24;
                else
                    MLT_l_j = -(24-MLT_l_j);
                end
            end
            if j==0
                dip_B_l = Calc_deltaB(L,MLT_l_j,deltaB0_l_c,L0_l,sigmal_l_c,MLTC_l+j*MLT_interval,delMLT_l_c,factor_l_c).*1e9;
            else
                dip_B_l = Calc_deltaB(L,MLT_l_j,deltaB0_l,L0_l,sigmaL_l,MLTC_l+j*MLT_interval,delMLT_l,factor_l).*1e9;
            end
            B11 = B11+dip_B_l; %[nT]
        end
         
        B2 = 0;
        B3 = 0;
        
        matrix_conv = [-3*sinTheta*cosTheta*cosPhi/sqrt(1+3*cosTheta^2), (sinTheta^2-2*cosTheta^2)*cosPhi/sqrt(1+3*cosTheta^2), -sinPhi; ...
            -3*sinTheta*cosTheta*sinPhi/sqrt(1+3*cosTheta^2), (sinTheta^2-2*cosTheta^2)*sinPhi/sqrt(1+3*cosTheta^2), cosPhi; ...
            (sinTheta^2-2*cosTheta^2)/sqrt(1+3*cosTheta^2), 3*sinTheta*cosTheta/sqrt(1+3*cosTheta^2), 0];        
        
        dip_B(:, i) = matrix_conv*[B11;B2;B3]./BE;
    end

    % Sum the background B and dip B.
    B = B + dip_B;

    % Normalize B.
    B = B/Scaling.B;
end

% the distribution of magnetic dip
function res = Calc_deltaB(L,MLT,deltaB0,L0,sigmaL,MLT0,delMLT,factor)
    res = deltaB0.*exp(-(L-L0).^2./(2.*sigmaL^2)).*((tanh(factor*(MLT-MLT0+delMLT))+tanh(factor*(-MLT+MLT0+delMLT)))/(2.*tanh(factor*delMLT)));
end

% the distribution of magnetic dip
function res = Calc_deltaB_deepen(L,MLT,deltaB0,L0,sigmaL,MLT0,sigMLT,t,t0,sigmaT)
    res = deltaB0.*exp(-(t-t0).^2./(2.*sigmaT.^2)).*exp(-(L-L0).^2./(2.*sigmaL^2)).*(1-((MLT-MLT0).^2)./sigMLT.^2).*exp(-((MLT-MLT0).^2)./(2.*sigMLT.^2));
end