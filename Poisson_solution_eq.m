%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                     SUBROUTINE                    %%%%%%%%%%%%
%------------------------  SOLVE POISSON EQUATION  -----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%  EQUILIBRIUM CONDITIONS  %%%%%%%%%%%%%%%%%%%%%%%%%
% Assumption: fn,fp = 0
%----------------------- Gummel's Iteration Method -----------------------%
%%%%%%             d2V/dx2 = -q/e (p - n - NA + ND)                   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVE NONLINEAR POISSON (EQUILIBRIUM) IN NORMALIZED FORM
%
% Normalization (PC1D-style but simplified):
%   psi      = V / Vt                    (dimensionless potential)
%   n_norm   = n / ni
%   p_norm   = p / ni
%   NA_norm  = NA / ni
%   ND_norm  = ND / ni
%   NB_norm  = ND_norm - NA_norm = (ND-NA)/ni
%
% Poisson (1D) in these variables:
%   d2 psi / dx^2 = - (q*ni / (eps*Vt)) * ( p_norm - n_norm + NB_norm )
%
% In equilibrium: n_norm = exp( +psi ), p_norm = exp( -psi )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eq = Poisson_solution_eq(DL,Lps,Lns,Vt)

    %------------------------ GLOBAL / INPUTS ----------------------------%
    % Call Global Variables
    [~,~,q,~,~,e0,~,~,~] = GlobalVariables();

    % Call Input Data
    [~,~,~,~,~,material,~, NA,ND,~, Tanode, ~,~,~,~,~,~,~,~,~,~,~,~,~,~,~ ...
        ,~,~,~, ~, ~] = Inputs();

    % Call materials
    [e,mn,mp,~,~,Eg,ni,~,~,~,~,~, ~] = SemiconductorMaterials(e0,Tanode,material);

    % Call numerical input
    [toll,max_iter] = Numerical_Input();

    %---------------------------------------------------------------------%
    %--------------------------- SETUP GRID -------------------------------%
    %---------------------------------------------------------------------%
    q_e = q / e;                  % q / permittivity

    % Define device in scaled coordinates
    xp_   = Lps / DL;
    xn_   = Lns / DL;
    xend  = Lps + Lns;                                                     % Total device length in cm
    xend_ = xend / DL;

    nodes = 300;                                                           % Number of nodes (uniform grid)

    % Define grid (scaled)
    x_  = linspace(0,xend_,nodes);
    dx_ = xend_ / (nodes - 1);

    %---------------------------------------------------------------------%
    %-------------------------- DOPING PROFILE ---------------------------%
    %---------------------------------------------------------------------%
    NAA = zeros(nodes,1);
    NDD = zeros(nodes,1);
    NB_norm  = zeros(nodes,1);

    for i = 1:nodes
        if x_(i) <= (xp_+xn_)/2       % p-region, from 0 to xp
            NAA(i) = NA;
            NDD(i) = 0;
        else                 % n-region, from xp to xp+xn
            NAA(i) = 0;
            NDD(i) = ND;
        end
        NB_norm(i) = (NDD(i) - NAA(i)) / ni;
    end


    %---------------------------------------------------------------------%
    %------------------------- INITIALIZATION ----------------------------%
    %---------------------------------------------------------------------%
    it_Num    = 0;
    iteration = zeros(max_iter,1);
    Error     = zeros(max_iter,1);

    V   = zeros(nodes,1);   % potential in volts
    V_  = zeros(nodes,1);   % normalized potential V/Vt

    pold = zeros(nodes,1);
    nold = zeros(nodes,1);

    %---------------------------------------------------------------------%
    %--------------------------- BOUNDARY CONDS --------------------------%
    %---------------------------------------------------------------------%
    % Dirichlet BCs
    V_left  =  Vt * log(ni / NA);
    V_right =  Vt * log(ND / ni);
    
    V(1)     = V_left;
    V(nodes) = V_right;
    
    % Equilibrium carrier BCs at contacts
    pold(1)     = NA;
    pold(nodes) = ni^2 / ND;
    nold(1)     = ni^2 / NA;
    nold(nodes) = ND;

    % Scaled BCs
    Vl_    = V_left / Vt;
    Vr_    = V_right / Vt;
    V_(1)  = Vl_;
    V_(nodes) = Vr_;

    % Physical x (unscaled) – can be computed once
    x = x_ * DL;

    %---------------------------------------------------------------------%
    %---------------------------- ITERATION ------------------------------%
    %---------------------------------------------------------------------%
    factor = -DL^2 * dx_^2 * q_e * ni / Vt;

    Error_N = inf;   % initialize > toll

    while (Error_N > toll)
        if it_Num >= max_iter
            fprintf(2,'************************************\n');
            fprintf(2,'Max number of iterations is reached!\n');
            fprintf(2,'Solution cannot be found!\n');
            fprintf(2,'Recheck::::::::::::::::::\n');
            fprintf(2,'a) Number of grid nodes !\n');
            fprintf(2,'b) Number of iterations !\n');
            fprintf(2,'************************************\n');
            break;
        end

        % Normalize potential
        V_    = V / Vt;
        Vold_ = V_;

        % Quasi-linearization step for Poisson's equation
        for i = 2:(nodes-1)
            pold(i) = exp(-Vold_(i));   % p/ni
            nold(i) = exp( Vold_(i));   % n/ni

            V_(i) = 0.5 * ( V_(i-1) + V_(i+1) ...
                        - factor * ( pold(i) - nold(i) + NB_norm(i) ) );
        end

        % Unscale
        V = V_ * Vt;

        % Error (absolute max difference in normalized potential)
        Error_N            = max(abs(V_ - Vold_));
        it_Num             = it_Num + 1;
        iteration(it_Num)  = it_Num;
        Error(it_Num)      = Error_N;
    end

    % Trim iteration/error vectors to actual length
    iteration = iteration(1:it_Num);
    Error     = Error(1:it_Num);

    %---------------------------------------------------------------------%
    %------------------ EQUILIBRIUM CARRIERS / OUTPUT --------------------%
    %---------------------------------------------------------------------%
    Veq   = V;
    Vnorm = Veq / Vt;

    peq = ni .* exp(-Vnorm);
    neq = ni .* exp( Vnorm);
    %-----------------------------------------------------------------%
    %----------------------- CHARGE DENSITY --------------------------%
    %-----------------------------------------------------------------%
    % Physical charge density ρ(x) = q (p - n + ND - NA)  [C/cm^3]
    rho = q * (peq - neq + NDD(:) - NAA(:));       % column vectors

    % Normalized charge density (to ni)
    n_norm  = neq / ni;
    p_norm  = peq / ni;
    ND_norm = NDD(:) / ni;
    NA_norm = NAA(:) / ni;
    rho_norm = p_norm - n_norm + ND_norm - NA_norm;

    fprintf('\n----- CHARGE MONITOR -----\n');
    fprintf(' Max |rho|      = %e C/cm^3\n', max(abs(rho)));
    fprintf(' Max |rho_norm| = %e (normalized)\n', max(abs(rho_norm)));

    % Net total charge (integral over device) – should be ~0 in eq.
    Dx = x(2) - x(1);
    Q_total = sum(rho) * Dx;    % [C/cm^2] total sheet charge
    fprintf(' Total integrated charge Q = %e C/cm^2\n', Q_total);
    fprintf('---------------------------\n');

    %---------------------------------------------------------------------%
    %--------------------- ELECTROSTATIC POTENTIAL -----------------------%
    %---------------------------------------------------------------------%
    % --- Compute electric field ---
    E = zeros(nodes,1);
    
    % dx in physical units (you already have this but recompute to be explicit)
    dx = x(2) - x(1);
    
    % central differences for interior nodes
    for i = 2:nodes-1
        E(i) = -(Veq(i+1) - Veq(i-1)) / (2*dx);
    end
    
    % forward difference at left boundary
    E(1) = -(Veq(2) - Veq(1)) / dx;
    
    % backward difference at right boundary
    E(nodes) = -(Veq(nodes) - Veq(nodes-1)) / dx;
    

    %-----------------------------------------------------------------%
    %----------------- SCHARFETTER–GUMMEL CURRENTS -------------------%
    %-----------------------------------------------------------------%
    % Mobilities (use your own values or read from Inputs())
    mu_n = mn;                       % [cm^2/V/s] 
    mu_p = mp;                       % [cm^2/V/s]  

    Dx = x(2) - x(1);                % physical spacing [cm]
    pref_n = q * mu_n * Vt / Dx;     % q Dn / Dx
    pref_p = q * mu_p * Vt / Dx;     % q Dp / Dx

    % SG edge currents between node i and i+1
    Jn_edge = zeros(nodes-1,1);      % J_n at i+1/2
    Jp_edge = zeros(nodes-1,1);      % J_p at i+1/2

    for i = 1:nodes-1
        psiL = Veq(i)   / Vt;
        psiR = Veq(i+1) / Vt;
        dpsi = psiR - psiL;

        Bp = bernoulliF(dpsi);       % B(+Δψ)
        Bm = bernoulliF(-dpsi);      % B(−Δψ)

        % Electron SG flux: Jn(i+1/2)
        %   Jn =  q Dn/Dx [ n_{i+1} B(Δψ) − n_i B(−Δψ) ]
        Jn_edge(i) = pref_n * ( neq(i+1)*Bp - neq(i)*Bm );

        % Hole SG flux: Jp(i+1/2)
        %   Jp = -q Dp/Dx [ p_{i+1} B(−Δψ) − p_i B(Δψ) ]
        Jp_edge(i) = -pref_p * ( peq(i+1)*Bm - peq(i)*Bp );
    end

    % Map edge currents to nodes for plotting (simple averaging)
    Jn = zeros(nodes,1);
    Jp = zeros(nodes,1);

    Jn(1)       = Jn_edge(1);
    Jn(2:end-1) = 0.5*(Jn_edge(1:end-1) + Jn_edge(2:end));
    Jn(end)     = Jn_edge(end);

    Jp(1)       = Jp_edge(1);
    Jp(2:end-1) = 0.5*(Jp_edge(1:end-1) + Jp_edge(2:end));
    Jp(end)     = Jp_edge(end);

    %-----------------------------------------------------------------%
    %------------ CLASSICAL DRIFT & DIFFUSION COMPONENTS -------------%
    %-----------------------------------------------------------------%

    % Carrier gradients dn/dx and dp/dx on nodes
    dn_dx = zeros(nodes,1);
    dp_dx = zeros(nodes,1);

    for i = 2:nodes-1
        dn_dx(i) = (neq(i+1) - neq(i-1)) / (2*Dx);   % [cm^-4]
        dp_dx(i) = (peq(i+1) - peq(i-1)) / (2*Dx);
    end
    dn_dx(1)     = (neq(2)   - neq(1))   / Dx;
    dn_dx(nodes) = (neq(end) - neq(end-1)) / Dx;
    dp_dx(1)     = (peq(2)   - peq(1))   / Dx;
    dp_dx(nodes) = (peq(end) - peq(end-1)) / Dx;

    % Electron drift & diffusion (node-based, classical)
    Jn_drift =  q * mu_n .* neq .* E;         % q μ_n n E        [A/cm^2]
    Jn_diff  =  q * mu_n * Vt .* dn_dx;       % q D_n dn/dx      [A/cm^2]
    Jn_tot_classic = Jn_drift + Jn_diff;      % ≈ 0 at equilibrium

    % Hole drift & diffusion (note the sign in drift)
    Jp_drift = -q * mu_p .* peq .* E;         % -q μ_p p E       [A/cm^2]
    Jp_diff  =  q * mu_p * Vt .* dp_dx;       % q D_p dp/dx      [A/cm^2]
    Jp_tot_classic = Jp_drift + Jp_diff;      % ≈ 0 at equilibrium

    % For reference: SG totals on nodes (already computed)
    Jn_tot_SG = Jn;     % from SG edge flux
    Jp_tot_SG = Jp;

    %-----------------------------------------------------------------%
    %--------------------------- CURRENT MONITOR ----------------------%
    %-----------------------------------------------------------------%

    % Total SG current at every edge (i -> i+1)
    J_edge_total = Jn_edge + Jp_edge;   % size: nodes-1

    % PC1D-style reporting:
    J_left  = J_edge_total(1);          % current at left contact (node 1/2)
    J_right = J_edge_total(end);        % current at right contact (node N-1/N)

    % A robust device current: average of all edge currents
    J_device = mean(J_edge_total);      % should be ~0 in equilibrium

    % Show in command window (just like PC1D)
    fprintf('\n----- Current Monitor (SG) -----\n');
    fprintf('  Left contact current   : %+e A/cm^2\n', J_left);
    fprintf('  Right contact current  : %+e A/cm^2\n', J_right);
    fprintf('  Average device current : %+e A/cm^2\n', J_device);
    fprintf('--------------------------------\n');
    %-----------------------------------------------------------------%
    %---------------------- DOPING FOR PLOTS -------------------------%
    %-----------------------------------------------------------------%
    ND_profile = NDD;    % donors [cm^-3]
    NA_profile = NAA;    % acceptors [cm^-3]
    %---------------------------------------------------------------------%
    %-----------------------------------------------------------------%
    %-------------------------- ENERGY BANDS --------------------------%
    %-----------------------------------------------------------------%
    psi_eq = Vnorm; % Ef - Ei = psi * Vt (in eV, since 1 eV ↔ 1 V for q)
    % at equilibrium we set EF as reference → fn = fp = 0
    fn_eq = zeros(size(psi_eq));
    fp_eq = zeros(size(psi_eq));

    Ei    = -Veq;            % Ei(x) in eV (since Veq is in volts)
    
    Ec    = Ei + Eg/2;       % conduction band [eV]
    Ev    = Ei - Eg/2;       % valence band [eV]
    Ef    = zeros(size(x));  % Fermi level at 0 eV (equilibrium)


    % You can shift everything so that bulk intrinsic or Fermi at left is 0
    shift = Ef(round(nodes/2));   % or Ef(1), or mean(Ef)
    Ec_plot = Ec - shift;
    Ei_plot = Ei - shift;
    Ev_plot = Ev - shift;
    Ef_plot = Ef - shift;
    %----------------------------- MONITORS ------------------------------%
    %---------------------------------------------------------------------%
    if Error_N < toll
        fprintf(1,'***********************************\n');
        fprintf(1,'****** Solution is converged! *****\n');
        fprintf(1,'Proceeding to non-equilibrium case!\n');
        fprintf(1,'***********************************\n');
    end
    %---------------------------------------------------------------------%
    %---------------------------- STRUCTURES -----------------------------%
    %---------------------------------------------------------------------%
    eq = struct();
    % base
    eq.x         = x;
    eq.Veq       = Veq;
    eq.Vnorm     = Vnorm;
    eq.peq       = peq;
    eq.neq       = neq;
    eq.n_norm    = n_norm;
    eq.p_norm    = p_norm;
    eq.NA_norm   = NA_norm;
    eq.ND_norm   = ND_norm;
    eq.NB_norm   = NB_norm;
    eq.rho       = rho;
    eq.rho_norm  = rho_norm;
    eq.Q_total   = Q_total;

    % currents / field
    eq.E               = E;
    eq.Jn_SG           = Jn_tot_SG;
    eq.Jp_SG           = Jp_tot_SG;
    eq.Jn_drift        = Jn_drift;
    eq.Jn_diff         = Jn_diff;
    eq.Jn_tot_classic  = Jn_tot_classic;
    eq.Jp_drift        = Jp_drift;
    eq.Jp_diff         = Jp_diff;
    eq.Jp_tot_classic  = Jp_tot_classic;
    eq.J_left          = J_left;
    eq.J_right         = J_right;
    eq.J_device        = J_device;

    % doping
    eq.NA_profile = NA_profile;
    eq.ND_profile = ND_profile;

    % bands
    eq.Ec_plot = Ec_plot;
    eq.Ei_plot = Ei_plot;
    eq.Ev_plot = Ev_plot;
    eq.Ef_plot = Ef_plot;

    % dimensionless PC1D-style
    eq.psi_eq = psi_eq;
    eq.fn_eq  = fn_eq;
    eq.fp_eq  = fp_eq;

    % convergence
    eq.iteration = iteration;
    eq.Error     = Error;
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                   END SUBROUTINE                  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B = bernoulliF(eta)
% Numerically safe Bernoulli function
%   B(eta) = eta / (exp(eta) - 1)
% with a series expansion for small |eta|.

    B = zeros(size(eta));
    small = abs(eta) < 1e-6;

    % series expansion: 1 - eta/2 + eta^2/12 - ...
    B(small)  = 1 - eta(small)/2 + (eta(small).^2)/12;

    % direct formula for larger arguments
    B(~small) = eta(~small) ./ (exp(eta(~small)) - 1);
end
