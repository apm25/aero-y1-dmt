clear
clc

R = 0.243;
B = 2;
Uinf = 10;
nElements = 24;
nu = 1.48e-5;
rho = 1.225;
max_iters = 20000;

alpha_cyl  = deg2rad(-180:1:180)';
Cl_cyl_all = zeros(size(alpha_cyl));
Cd_cyl_all = ones(size(alpha_cyl)) * 1.2;

raw_sg   = readmatrix("sg6040_polars.csv");
raw_s834 = readmatrix("s834_polars_new.csv");

active_condition = 'blended';
switch active_condition
    case 'clean',    idx = 1:91;    disp('Running: FREE TRANSITION');
    case 'blended',  idx = 92:182;  disp('Running: BLENDED');
    case 'rough',    idx = 183:273; disp('Running: FULLY TURBULENT');
end

sg_20  = raw_sg(raw_sg(:,1) == 20000,  :);
sg_50  = raw_sg(raw_sg(:,1) == 50000,  :);
sg_100 = raw_sg(raw_sg(:,1) == 100000, :);

s834_20  = raw_s834(raw_s834(:,1) == 20000,  :);
s834_50  = raw_s834(raw_s834(:,1) == 50000,  :);
s834_100 = raw_s834(raw_s834(:,1) == 100000, :);

mean_chord_approx = 0.1;
AR_blade          = R / mean_chord_approx;

[alpha_sg20,  Cl_sg20_ext,  Cd_sg20_ext]  = viterna_extend(deg2rad(sg_20(idx,2)),    sg_20(idx,3),    sg_20(idx,4),    AR_blade);
[alpha_sg50,  Cl_sg50_ext,  Cd_sg50_ext]  = viterna_extend(deg2rad(sg_50(idx,2)),    sg_50(idx,3),    sg_50(idx,4),    AR_blade);
[alpha_sg100, Cl_sg100_ext, Cd_sg100_ext] = viterna_extend(deg2rad(sg_100(idx,2)),   sg_100(idx,3),   sg_100(idx,4),   AR_blade);

[alpha_s834_20,  Cl_s834_20_ext,  Cd_s834_20_ext]  = viterna_extend(deg2rad(s834_20(idx,2)),  s834_20(idx,3),  s834_20(idx,4),  AR_blade);
[alpha_s834_50,  Cl_s834_50_ext,  Cd_s834_50_ext]  = viterna_extend(deg2rad(s834_50(idx,2)),  s834_50(idx,3),  s834_50(idx,4),  AR_blade);
[alpha_s834_100, Cl_s834_100_ext, Cd_s834_100_ext] = viterna_extend(deg2rad(s834_100(idx,2)), s834_100(idx,3), s834_100(idx,4), AR_blade);

% 1 = Cylinder, 2 = SG6040, 3 = S834
alpha20_dict  = {alpha_cyl, alpha_sg20,  alpha_s834_20};
alpha50_dict  = {alpha_cyl, alpha_sg50,  alpha_s834_50};
alpha100_dict = {alpha_cyl, alpha_sg100, alpha_s834_100};

Cl20_dict  = {Cl_cyl_all, Cl_sg20_ext,  Cl_s834_20_ext};
Cd20_dict  = {Cd_cyl_all, Cd_sg20_ext,  Cd_s834_20_ext};

Cl50_dict  = {Cl_cyl_all, Cl_sg50_ext,  Cl_s834_50_ext};
Cd50_dict  = {Cd_cyl_all, Cd_sg50_ext,  Cd_s834_50_ext};

Cl100_dict = {Cl_cyl_all, Cl_sg100_ext, Cl_s834_100_ext};
Cd100_dict = {Cd_cyl_all, Cd_sg100_ext, Cd_s834_100_ext};

rVal         = linspace(0.0163, R, nElements);
mu           = rVal / R;
elementWidth = rVal(2) - rVal(1);

r_hub = 0.025;

% 60% SG6040 + 40% S834 over lifting span
r_split = r_hub + 0.60*(R - r_hub);
r_tip_start = r_split;

airfoil_map = zeros(1, nElements);
airfoil_map(rVal <= r_hub)                  = 1; % cylinder
airfoil_map(rVal > r_hub & rVal < r_split)  = 2; % SG6040
airfoil_map(rVal >= r_split)                = 3; % S834 tip

fprintf('r_hub=%.4f, r_split=%.4f, r_tip_start=%.4f, R=%.4f\n', r_hub, r_split, r_tip_start, R);
fprintf('Map counts: cyl=%d, SG6040=%d, S834=%d\n', ...
    nnz(airfoil_map==1), nnz(airfoil_map==2), nnz(airfoil_map==3));

design_TSRs_to_test  = [3.5 4 4.5 5];
operating_TSRs       = 0.5:0.1:7;
nOp  = length(operating_TSRs);
nDes = length(design_TSRs_to_test);

all_Cp_results              = zeros(nDes, nOp);
all_Torque_results          = zeros(nDes, nOp);
all_Thrust_results          = zeros(nDes, nOp);
all_TangentialForce_results = zeros(nDes, nOp);

spanwise_Thrust_TSR35    = zeros(nElements, nOp);
spanwise_ForceTang_TSR35 = zeros(nElements, nOp);
spanwise_Torque_TSR35    = zeros(nElements, nOp);
spanwise_Thrust_TSR40    = zeros(nElements, nOp);
spanwise_ForceTang_TSR40 = zeros(nElements, nOp);
spanwise_Torque_TSR40    = zeros(nElements, nOp);
spanwise_Thrust_TSR45    = zeros(nElements, nOp);
spanwise_ForceTang_TSR45 = zeros(nElements, nOp);
spanwise_Torque_TSR45    = zeros(nElements, nOp);
spanwise_Thrust_TSR50    = zeros(nElements, nOp);
spanwise_ForceTang_TSR50 = zeros(nElements, nOp);
spanwise_Torque_TSR50    = zeros(nElements, nOp);

for d = 1:nDes

    Design_TSR = design_TSRs_to_test(d);
    lambda_r   = Design_TSR .* mu;

    % Target Optimum Tip Performance (S834)
    LD_s834 = s834_100(idx,3) ./ s834_100(idx,4);
    LD_s834(s834_100(idx,3) < 0.3) = 0;
    [~, opt_idx] = max(LD_s834);
    ClDes    = s834_100(idx(opt_idx),3);
    alphaDes = deg2rad(s834_100(idx(opt_idx),2));

    % Target Zero-Lift Angle (SG6040 root)
    Cl_col    = sg_50(idx, 3);
    alpha_col = deg2rad(sg_50(idx, 2));
    zc = find(diff(sign(Cl_col)), 1, 'first');
    if ~isempty(zc)
        alpha0 = alpha_col(zc) + (0 - Cl_col(zc)) * ...
                 (alpha_col(zc+1) - alpha_col(zc)) / (Cl_col(zc+1) - Cl_col(zc));
    else
        alpha0 = 0;
    end

    r_blade_start      = 0.0163;   
    chord_at_r_start   = 0.0102;    
    chord_at_hub_outer = 0.01655;   

    r_aero_start = 0.06;
    idx_tip      = rVal > r_tip_start;
    idx_root     = rVal < r_aero_start;
    aero_mask    = ~idx_root & ~idx_tip;

    % GLAUERT CHORD
    chordDisttheory = (16.*pi.*R^2) ./ (9.*B.*Design_TSR^2.*rVal.*ClDes);
    
    % Manufacturing cap
    c_max_aero  = 0.09;
    chordCapped = min(chordDisttheory, c_max_aero);
    c_R = interp1(rVal, chordCapped, r_aero_start, 'linear', 'extrap');
    c_R = min(c_R, c_max_aero);   
    c = chordCapped;   

    % TIP TAPER
    c(idx_tip) = chordCapped(idx_tip) .* ...
                 sqrt(1 - ((rVal(idx_tip) - r_tip_start) ./ (R - r_tip_start)).^2);
    c(idx_tip) = max(c(idx_tip), 0.006);

    % ROOT BLENDING
    idx_hub_body = (rVal <= r_hub);
    for i = find(idx_hub_body)
        t_i = (rVal(i) - r_blade_start) / (r_hub - r_blade_start);
        t_i = max(0, min(1, t_i));
        c(i) = chord_at_r_start + t_i * (chord_at_hub_outer - chord_at_r_start);
    end

    idx_root_blend = (rVal > r_hub) & (rVal < r_aero_start);
    blend_fraction = (rVal(idx_root_blend) - r_hub) ./ (r_aero_start - r_hub);
    easing         = 0.5 .* (1 - cos(pi .* blend_fraction));
    c(idx_root_blend) = chord_at_hub_outer + easing .* (c_R - chord_at_hub_outer);
    c(idx_root_blend) = min(c(idx_root_blend), c_max_aero);   

    % ENFORCE MONOTONIC CHORD (Aero Section)
    aero_idx = find(aero_mask);
    for k = 2:length(aero_idx)
        if c(aero_idx(k)) > c(aero_idx(k-1))
            c(aero_idx(k)) = c(aero_idx(k-1));
        end
    end

    % TWIST
    beta_ideal = atan(2 ./ (3 .* Design_TSR .* mu)) - alphaDes;
    max_twist  = deg2rad(35);
    betaArray  = max(beta_ideal, -max_twist);
    betaArray  = min(betaArray,   max_twist);
    betaArray(airfoil_map == 1) = beta_ideal(airfoil_map == 1);
    alpha0_span = alpha0 .* ones(1, nElements);

    localSolidity = (B .* c) ./ (2 .* pi .* rVal);

    fprintf('\n--- Design TSR = %.1f ---\n', Design_TSR);
    fprintf('ClDes=%.3f  alphaDes=%.2f deg  alpha0=%.2f deg\n', ...
        ClDes, rad2deg(alphaDes), rad2deg(alpha0));
    fprintf('chord:  r=16.55mm: %.2fmm  r=25mm: %.2fmm  aero_start: %.1fmm  mid: %.1fmm  tip: %.1fmm\n', ...
        c(1)*1000, c(find(rVal>=r_hub,1))*1000, c(find(~idx_root,1))*1000, c(round(nElements/2))*1000, c(end)*1000);
    fprintf('beta:   root=%.1f  mid=%.1f  tip=%.1f deg\n', ...
        rad2deg(betaArray(1)), rad2deg(betaArray(round(nElements/2))), rad2deg(betaArray(end)));
    
    figure;
    plot(rVal*1000, c*1000, 'b-', 'LineWidth', 2); hold on;
    plot(rVal*1000, chordCapped*1000, 'r--', 'LineWidth', 1);
    xline(r_blade_start*1000, 'k:', 'Blade start (16.55mm)', 'LineWidth', 1);
    xline(r_hub*1000,         'k--','Hub outer (25mm)',      'LineWidth', 1);
    xline(r_aero_start*1000,  'k-.','Aero start (60mm)',     'LineWidth', 1);
    yline(chord_at_r_start*1000,   'g--', '10.2mm at r=16.55mm', 'LineWidth', 1);
    yline(chord_at_hub_outer*1000, 'm--', '16.55mm at r=25mm',   'LineWidth', 1);
    xlabel('r (mm)'); ylabel('chord (mm)');
    title(['Chord - Design TSR = ', num2str(Design_TSR)]);
    legend('Physical chord','Capped Glauert (80mm max)', ...
           'Blade start','Hub outer','Aero start', ...
           'Chord=10.2mm','Chord=16.55mm','Location','northeast');
    ylim([0 90]); grid on;

    % INNER BEM LOOP
    for op = 1:nOp
        Operating_TSR   = operating_TSRs(op);
        omega           = (Operating_TSR * Uinf) / R;
        localSpeedRatio = Operating_TSR .* mu;

        a           = (1/3) .* ones(1, nElements);
        a_dash      = zeros(1, nElements);
        phi         = atan2((1 - a), (localSpeedRatio .* (1 + a_dash)));
        phi         = min(phi, deg2rad(89));

        iterating   = true;
        itercount   = 0;
        relax_a     = 0.25;
        relax_phi   = 0.25;
        res_a_prev  = 1.0;
        res_a       = 1.0;
        res_adash   = 1.0;
        res_phi     = 1.0;
        a_prev2     = a;
        a_prev1     = a;
        high_a_prev = false(1, nElements);

        while iterating
            alpha = phi - betaArray;

            W_local  = sqrt((Uinf.*(1-a)).^2 + (omega.*rVal.*(1+a_dash)).^2);
            Re_local = (W_local .* c) ./ nu;

            w20  = max(0, min(1, (50000  - Re_local) ./ (50000  - 20000)));
            w100 = max(0, min(1, (Re_local - 50000)  ./ (100000 - 50000)));
            w50  = max(0, 1 - w20 - w100);

            Cl_20_eval  = zeros(1,nElements); Cl_50_eval  = zeros(1,nElements); Cl_100_eval = zeros(1,nElements);
            Cd_20_eval  = zeros(1,nElements); Cd_50_eval  = zeros(1,nElements); Cd_100_eval = zeros(1,nElements);

            for foil_idx = 1:3
                mask = (airfoil_map == foil_idx);
                if any(mask)
                    sa    = alpha(mask);
                    lo20  = alpha20_dict{foil_idx}(1);  hi20  = alpha20_dict{foil_idx}(end);
                    lo50  = alpha50_dict{foil_idx}(1);  hi50  = alpha50_dict{foil_idx}(end);
                    lo100 = alpha100_dict{foil_idx}(1); hi100 = alpha100_dict{foil_idx}(end);

                    sa_20  = max(lo20,  min(hi20,  sa));
                    sa_50  = max(lo50,  min(hi50,  sa));
                    sa_100 = max(lo100, min(hi100, sa));

                    Cl_20_eval(mask)  = interp1(alpha20_dict{foil_idx},  Cl20_dict{foil_idx},  sa_20,  'linear');
                    Cd_20_eval(mask)  = interp1(alpha20_dict{foil_idx},  Cd20_dict{foil_idx},  sa_20,  'linear');
                    Cl_50_eval(mask)  = interp1(alpha50_dict{foil_idx},  Cl50_dict{foil_idx},  sa_50,  'linear');
                    Cd_50_eval(mask)  = interp1(alpha50_dict{foil_idx},  Cd50_dict{foil_idx},  sa_50,  'linear');
                    Cl_100_eval(mask) = interp1(alpha100_dict{foil_idx}, Cl100_dict{foil_idx}, sa_100, 'linear');
                    Cd_100_eval(mask) = interp1(alpha100_dict{foil_idx}, Cd100_dict{foil_idx}, sa_100, 'linear');
                end
            end

            Cl_current = w20.*Cl_20_eval + w50.*Cl_50_eval + w100.*Cl_100_eval;
            Cd_current = w20.*Cd_20_eval + w50.*Cd_50_eval + w100.*Cd_100_eval;

            c_over_r_safe = min(c ./ rVal, 0.35);
            Cl_attached   = 2.*pi.*(alpha - alpha0_span);
            delta_Cl      = max(0, Cl_attached - Cl_current);
            Cl_current    = Cl_current + 3.*(c_over_r_safe).^2.*delta_Cl.*(1-mu);

            Cl_current(airfoil_map == 1) = 0;
            Cd_current(airfoil_map == 1) = 1.2;

            sin_phi_safe = max(abs(sin(phi)), 1e-4);
            f_tip  = (B.*(R - rVal))     ./ (2.*rVal.*sin_phi_safe);
            f_root = (B.*(rVal - r_hub)) ./ (2.*rVal.*sin_phi_safe);
            f_root = max(f_root, 0);
            F_tip  = (2/pi).*acos(exp(-max(f_tip,  0)));
            F_root = (2/pi).*acos(exp(-max(f_root, 0)));
            F      = max(F_tip.*F_root, 1e-4);

            Ct = Cl_current.*sin(phi) - Cd_current.*cos(phi);
            Cn = Cl_current.*cos(phi) + Cd_current.*sin(phi);

            CT_BEM   = (localSolidity.*Cn.*(1-a).^2) ./ max(sin(phi).^2, 1e-8);
            CT_BEM_p = max(CT_BEM, 0);
            CT_cross = 0.96.*F;

            high_a = (CT_BEM_p > 1.02.*CT_cross) | ...
                     (high_a_prev & (CT_BEM_p > 0.98.*CT_cross));
            low_a  = ~high_a;
            high_a_prev = high_a;

            K     = (localSolidity.*Cn) ./ (4.*F.*max(sin(phi).^2, 1e-8));
            a_new = zeros(1, nElements);
            a_new(low_a) = K(low_a) ./ (1 + K(low_a));

            F_h      = F(high_a);
            CT_h     = CT_BEM_p(high_a);
            radicand = max(0, CT_h.*(50-36.*F_h) + 12.*F_h.*(3.*F_h-4));
            a_new(high_a) = (18.*F_h - 20 - 3.*sqrt(radicand)) ./ (36.*F_h - 50);

            a_new  = max(-0.5, min(a_new, 0.95));
            a_prev = a;
            a      = relax_a.*a_new + (1-relax_a).*a_prev;
            a      = max(-0.5, min(a, 0.95));

            if itercount > 3 && mod(itercount, 2) == 0
                delta1       = a_prev1 - a_prev2;
                delta2       = a_prev  - a_prev1;
                denom_aitken = delta2  - delta1;
                aitken_mask  = abs(denom_aitken) > 1e-12;
                if any(aitken_mask)
                    a_aitken = a_prev;
                    a_aitken(aitken_mask) = a_prev(aitken_mask) - ...
                        delta1(aitken_mask).^2 ./ denom_aitken(aitken_mask);
                    a_aitken = max(-0.5, min(a_aitken, 0.95));
                    if max(abs(a_aitken - a_prev)) < max(abs(a - a_prev))
                        a = a_aitken;
                    end
                end
            end
            a_prev2 = a_prev1;
            a_prev1 = a_new;

            if itercount > 10
                if res_a > res_a_prev * 0.99
                    relax_a   = max(relax_a   * 0.85, 0.05);
                    relax_phi = max(relax_phi * 0.85, 0.05);
                else
                    relax_a   = min(relax_a   * 1.02, 0.25);
                    relax_phi = min(relax_phi * 1.02, 0.25);
                end
            end
            res_a_prev = res_a;

            denom_dash = 4.*F.*sin(phi).*cos(phi);
            denom_dash = sign(denom_dash).*max(abs(denom_dash), 1e-6);
            K_dash     = (localSolidity.*Ct) ./ denom_dash;
            a_dash_new = K_dash ./ (1 - K_dash);
            a_dash_new = max(-0.5, min(a_dash_new, 0.5));

            a_dash_prev = a_dash;
            a_dash      = relax_a.*a_dash_new + (1-relax_a).*a_dash_prev;

            phi_prev = phi;
            phi_new  = atan2((1-a), (localSpeedRatio.*(1+a_dash)));
            phi_new  = min(phi_new, deg2rad(89));
            phi      = relax_phi.*phi_new + (1-relax_phi).*phi_prev;

            res_a     = max(abs(a_new      - a_prev));
            res_adash = max(abs(a_dash_new - a_dash_prev));
            res_phi   = max(abs(phi_new    - phi_prev));

            if res_a < 1e-4 && res_adash < 1e-4 && res_phi < 1e-4
                iterating = false;
            end

            itercount = itercount + 1;
            if itercount >= max_iters
                iterating = false;
            end
        end

        W                    = sqrt((omega.*rVal.*(1+a_dash)).^2 + (Uinf.*(1-a)).^2);
        deltaLift            = 0.5.*rho.*W.^2.*Cl_current.*c.*elementWidth;
        deltaDrag            = 0.5.*rho.*W.^2.*Cd_current.*c.*elementWidth;
        axialThrust          = B .* (deltaLift.*cos(phi) + deltaDrag.*sin(phi));
        tangentialForce      = B .* (deltaLift.*sin(phi) - deltaDrag.*cos(phi));
        tangentialTorque     = tangentialForce .* rVal;

        totalThrust          = sum(axialThrust);
        totalTangentialForce = sum(tangentialForce);
        totalTorque          = sum(tangentialTorque);
        turbinePower         = totalTorque * omega;
        availablePower       = 0.5 * rho * pi * R^2 * Uinf^3;
        powerCoefficient     = turbinePower / availablePower;

        all_Cp_results(d, op)              = powerCoefficient;
        all_Torque_results(d, op)          = totalTorque;
        all_Thrust_results(d, op)          = totalThrust;
        all_TangentialForce_results(d, op) = totalTangentialForce;

        if Design_TSR == 3.5
            spanwise_Thrust_TSR35(:,op)    = axialThrust'     ./ B;
            spanwise_ForceTang_TSR35(:,op) = tangentialForce' ./ B;
            spanwise_Torque_TSR35(:,op)    = tangentialTorque'./ B;
        elseif Design_TSR == 4.0
            spanwise_Thrust_TSR40(:,op)    = axialThrust'     ./ B;
            spanwise_ForceTang_TSR40(:,op) = tangentialForce' ./ B;
            spanwise_Torque_TSR40(:,op)    = tangentialTorque'./ B;
        elseif Design_TSR == 4.5
            spanwise_Thrust_TSR45(:,op)    = axialThrust'     ./ B;
            spanwise_ForceTang_TSR45(:,op) = tangentialForce' ./ B;
            spanwise_Torque_TSR45(:,op)    = tangentialTorque'./ B;
        elseif Design_TSR == 5.0
            spanwise_Thrust_TSR50(:,op)    = axialThrust'     ./ B;
            spanwise_ForceTang_TSR50(:,op) = tangentialForce' ./ B;
            spanwise_Torque_TSR50(:,op)    = tangentialTorque'./ B;
        end
    end

    if Design_TSR == 3.5
        spanwise_Chord_TSR35 = c;
        spanwise_Twist_TSR35 = betaArray;
    elseif Design_TSR == 4.0
        spanwise_Chord_TSR40 = c;
        spanwise_Twist_TSR40 = betaArray;
    elseif Design_TSR == 4.5
        spanwise_Chord_TSR45 = c;
        spanwise_Twist_TSR45 = betaArray;
    elseif Design_TSR == 5.0
        spanwise_Chord_TSR50 = c;
        spanwise_Twist_TSR50 = betaArray;
    end
    disp(['Finished Design TSR: ', num2str(Design_TSR)]);
end

% FLAT PLATE SIMULATION
flat_TSRs       = [4.0, 4.5];
flat_plate_chord = 0.090;        
flat_plate_alpha = deg2rad(8);   

all_Cp_flat    = zeros(length(flat_TSRs), nOp);
all_Torque_flat= zeros(length(flat_TSRs), nOp);
all_Thrust_flat= zeros(length(flat_TSRs), nOp);

for fd = 1:length(flat_TSRs)
    FP_TSR   = flat_TSRs(fd);
    lambda_r_fp = FP_TSR .* mu;

    c_fp = flat_plate_chord .* ones(1, nElements);
    idx_hub_fp = (rVal <= r_hub);
    for i = find(idx_hub_fp)
        t_i = (rVal(i) - r_blade_start) / (r_hub - r_blade_start);
        t_i = max(0, min(1, t_i));
        c_fp(i) = chord_at_r_start + t_i*(chord_at_hub_outer - chord_at_r_start);
    end

    beta_fp = atan(2 ./ (3 .* FP_TSR .* mu)) - flat_plate_alpha;
    beta_fp = max(beta_fp, -deg2rad(35));
    beta_fp = min(beta_fp,  deg2rad(35));
    beta_fp(airfoil_map == 1) = atan(2 ./ (3 .* FP_TSR .* mu(airfoil_map==1))) - flat_plate_alpha;
    
    localSolidity_fp = (B .* c_fp) ./ (2 .* pi .* rVal);

    fprintf('\n--- FLAT PLATE TSR = %.1f ---\n', FP_TSR);
    
    for op = 1:nOp
        Operating_TSR   = operating_TSRs(op);
        omega           = (Operating_TSR * Uinf) / R;
        localSpeedRatio = Operating_TSR .* mu;

        a      = (1/3) .* ones(1, nElements);
        a_dash = zeros(1, nElements);
        phi    = atan2((1 - a), (localSpeedRatio .* (1 + a_dash)));
        phi    = min(phi, deg2rad(89));

        iterating   = true;
        itercount   = 0;
        relax_a     = 0.25;
        relax_phi   = 0.25;
        res_a_prev  = 1.0;
        res_a       = 1.0;
        res_adash   = 1.0;
        res_phi     = 1.0;
        a_prev2     = a;
        a_prev1     = a;
        high_a_prev = false(1, nElements);

        while iterating
            alpha_fp_local = phi - beta_fp;

            Cl_fp = 2.*pi.*sin(alpha_fp_local).*cos(alpha_fp_local);
            Cd_fp = 2.*pi.*sin(alpha_fp_local).^2 + 0.02;  
            
            Cl_fp(airfoil_map == 1) = 0;
            Cd_fp(airfoil_map == 1) = 1.2;

            sin_phi_safe = max(abs(sin(phi)), 1e-4);
            f_tip  = (B.*(R - rVal))     ./ (2.*rVal.*sin_phi_safe);
            f_root = (B.*(rVal - r_hub)) ./ (2.*rVal.*sin_phi_safe);
            f_root = max(f_root, 0);
            F_tip  = (2/pi).*acos(exp(-max(f_tip,  0)));
            F_root = (2/pi).*acos(exp(-max(f_root, 0)));
            F      = max(F_tip.*F_root, 1e-4);

            Ct_fp = Cl_fp.*sin(phi) - Cd_fp.*cos(phi);
            Cn_fp = Cl_fp.*cos(phi) + Cd_fp.*sin(phi);

            CT_BEM   = (localSolidity_fp.*Cn_fp.*(1-a).^2) ./ max(sin(phi).^2, 1e-8);
            CT_BEM_p = max(CT_BEM, 0);
            CT_cross = 0.96.*F;

            high_a = (CT_BEM_p > 1.02.*CT_cross) | (high_a_prev & (CT_BEM_p > 0.98.*CT_cross));
            low_a  = ~high_a;
            high_a_prev = high_a;

            K     = (localSolidity_fp.*Cn_fp) ./ (4.*F.*max(sin(phi).^2, 1e-8));
            a_new = zeros(1, nElements);
            a_new(low_a) = K(low_a) ./ (1 + K(low_a));

            F_h      = F(high_a);
            CT_h     = CT_BEM_p(high_a);
            radicand = max(0, CT_h.*(50-36.*F_h) + 12.*F_h.*(3.*F_h-4));
            a_new(high_a) = (18.*F_h - 20 - 3.*sqrt(radicand)) ./ (36.*F_h - 50);

            a_new  = max(-0.5, min(a_new, 0.95));
            a_prev = a;
            a      = relax_a.*a_new + (1-relax_a).*a_prev;
            a      = max(-0.5, min(a, 0.95));

            if itercount > 3 && mod(itercount, 2) == 0
                delta1       = a_prev1 - a_prev2;
                delta2       = a_prev  - a_prev1;
                denom_aitken = delta2  - delta1;
                aitken_mask  = abs(denom_aitken) > 1e-12;
                if any(aitken_mask)
                    a_aitken = a_prev;
                    a_aitken(aitken_mask) = a_prev(aitken_mask) - ...
                        delta1(aitken_mask).^2 ./ denom_aitken(aitken_mask);
                    a_aitken = max(-0.5, min(a_aitken, 0.95));
                    if max(abs(a_aitken - a_prev)) < max(abs(a - a_prev))
                        a = a_aitken;
                    end
                end
            end
            a_prev2 = a_prev1;
            a_prev1 = a_new;

            if itercount > 10
                if res_a > res_a_prev * 0.99
                    relax_a   = max(relax_a   * 0.85, 0.05);
                    relax_phi = max(relax_phi * 0.85, 0.05);
                else
                    relax_a   = min(relax_a   * 1.02, 0.25);
                    relax_phi = min(relax_phi * 1.02, 0.25);
                end
            end
            res_a_prev = res_a;

            denom_dash = 4.*F.*sin(phi).*cos(phi);
            denom_dash = sign(denom_dash).*max(abs(denom_dash), 1e-6);
            K_dash     = (localSolidity_fp.*Ct_fp) ./ denom_dash;
            a_dash_new = K_dash ./ (1 - K_dash);
            a_dash_new = max(-0.5, min(a_dash_new, 0.5));

            a_dash_prev = a_dash;
            a_dash      = relax_a.*a_dash_new + (1-relax_a).*a_dash_prev;

            phi_prev = phi;
            phi_new  = atan2((1-a), (localSpeedRatio.*(1+a_dash)));
            phi_new  = min(phi_new, deg2rad(89));
            phi      = relax_phi.*phi_new + (1-relax_phi).*phi_prev;

            res_a     = max(abs(a_new      - a_prev));
            res_adash = max(abs(a_dash_new - a_dash_prev));
            res_phi   = max(abs(phi_new    - phi_prev));

            if res_a < 1e-4 && res_adash < 1e-4 && res_phi < 1e-4
                iterating = false;
            end
            itercount = itercount + 1;
            if itercount >= max_iters
                iterating = false;
            end
        end

        W_fp             = sqrt((omega.*rVal.*(1+a_dash)).^2 + (Uinf.*(1-a)).^2);
        deltaLift_fp     = 0.5.*rho.*W_fp.^2.*Cl_fp.*c_fp.*elementWidth;
        deltaDrag_fp     = 0.5.*rho.*W_fp.^2.*Cd_fp.*c_fp.*elementWidth;
        axialThrust_fp   = B .* (deltaLift_fp.*cos(phi) + deltaDrag_fp.*sin(phi));
        tangForce_fp     = B .* (deltaLift_fp.*sin(phi) - deltaDrag_fp.*cos(phi));
        tangTorque_fp    = tangForce_fp .* rVal;

        all_Cp_flat(fd, op)     = sum(tangTorque_fp)*omega / availablePower;
        all_Torque_flat(fd, op) = sum(tangTorque_fp);
        all_Thrust_flat(fd, op) = sum(axialThrust_fp);
    end
    fprintf('Finished flat plate TSR=%.1f\n', FP_TSR);
end

wind_speeds   = 3:0.5:15;         
nWind         = length(wind_speeds);
all_Power_vs_Wind   = zeros(nDes, nWind);  
flat_Power_vs_Wind  = zeros(length(flat_TSRs), nWind);   

for d = 1:nDes
    Design_TSR = design_TSRs_to_test(d);
    [Cp_peak, peak_op_idx] = max(all_Cp_results(d,:));
    peak_op_TSR = operating_TSRs(peak_op_idx);
    for w = 1:nWind
        Uwind = wind_speeds(w);
        P_avail = 0.5 * rho * pi * R^2 * Uwind^3;
        all_Power_vs_Wind(d, w) = Cp_peak * P_avail;
    end
end

for fd = 1:length(flat_TSRs)
    [Cp_peak_fp, ~] = max(all_Cp_flat(fd,:));
    for w = 1:nWind
        Uwind = wind_speeds(w);
        P_avail = 0.5 * rho * pi * R^2 * Uwind^3;
        flat_Power_vs_Wind(fd, w) = Cp_peak_fp * P_avail;
    end
end

P_available_wind = 0.5 * rho * pi * R^2 .* wind_speeds.^3;
P_betz           = 16/27 .* P_available_wind;   

% CAD BLUEPRINT EXPORT
profile_names = strings(nElements, 1);
profile_names(airfoil_map == 1) = "Cylinder (10.2mm)";
profile_names(airfoil_map == 2) = "SG6040";
profile_names(airfoil_map == 3) = "S834";

beta_deg_35 = rad2deg(spanwise_Twist_TSR35)';
CAD_Table_35 = table(rVal', spanwise_Chord_TSR35', beta_deg_35, profile_names, ...
    'VariableNames', {'Radius_m','Chord_m','Twist_deg','Airfoil_Profile'});
writetable(CAD_Table_35, 'CAD_Blueprint_Blade_TSR35.csv');

beta_deg_40 = rad2deg(spanwise_Twist_TSR40)';
CAD_Table_40 = table(rVal', spanwise_Chord_TSR40', beta_deg_40, profile_names, ...
    'VariableNames', {'Radius_m','Chord_m','Twist_deg','Airfoil_Profile'});
writetable(CAD_Table_40, 'CAD_Blueprint_Blade_TSR40.csv');

beta_deg_45 = rad2deg(spanwise_Twist_TSR45)';
CAD_Table_45 = table(rVal', spanwise_Chord_TSR45', beta_deg_45, profile_names, ...
    'VariableNames', {'Radius_m','Chord_m','Twist_deg','Airfoil_Profile'});
writetable(CAD_Table_45, 'CAD_Blueprint_Blade_TSR45.csv');

beta_deg_50 = rad2deg(spanwise_Twist_TSR50)';
CAD_Table_50 = table(rVal', spanwise_Chord_TSR50', beta_deg_50, profile_names, ...
    'VariableNames', {'Radius_m','Chord_m','Twist_deg','Airfoil_Profile'});
writetable(CAD_Table_50, 'CAD_Blueprint_Blade_TSR50.csv');

save('results_theoretical_capped.mat');

% VITERNA EXTEND
function [alpha_ext, Cl_ext, Cd_ext] = viterna_extend(alpha_in, Cl_in, Cd_in, AR)
    [Cl_max, idx_stall] = max(Cl_in);
    alpha_stall = alpha_in(idx_stall);
    Cd_stall    = Cd_in(idx_stall);
    Cd_90 = min(1.11 + 0.018 * AR, 2.01);
    B1 = Cd_90; A1 = B1/2;
    B2 = (Cl_max - Cd_90.*sin(alpha_stall).*cos(alpha_stall)) ...
          .*sin(alpha_stall)./cos(alpha_stall).^2;
    A2 = (Cd_stall - Cd_90.*sin(alpha_stall).^2)/cos(alpha_stall);
    alpha_end = alpha_in(end);
    if alpha_end < pi/2
        alpha_pos = linspace(alpha_end, pi/2, 50)';
        Cl_pos    = A1.*sin(2.*alpha_pos) + A2.*(cos(alpha_pos).^2./sin(alpha_pos));
        Cd_pos    = B1.*sin(alpha_pos).^2 + B2.*cos(alpha_pos);
        alpha_full = [alpha_in; alpha_pos(2:end)];
        Cl_full    = [Cl_in;    Cl_pos(2:end)];
        Cd_full    = [Cd_in;    Cd_pos(2:end)];
    else
        alpha_full = alpha_in; Cl_full = Cl_in; Cd_full = Cd_in;
    end
    alpha_start = alpha_in(1);
    if alpha_start > -pi/2
        alpha_neg = linspace(-pi/2, alpha_start, 50)';
        a_mirror  = -alpha_neg;
        Cl_neg    = -(A1.*sin(2.*a_mirror) + A2.*(cos(a_mirror).^2./sin(a_mirror)));
        Cd_neg    =   B1.*sin(a_mirror).^2 + B2.*cos(a_mirror);
        alpha_full = [alpha_neg(1:end-1); alpha_full];
        Cl_full    = [Cl_neg(1:end-1);    Cl_full];
        Cd_full    = [Cd_neg(1:end-1);    Cd_full];
    end
    [alpha_ext, ia] = unique(alpha_full);
    Cl_ext = Cl_full(ia);
    Cd_ext = Cd_full(ia);
end