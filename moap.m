clear; clc; close all;
%% THAM SO
mstar_rel = 0.067;
eps_static = 13.18;
eps_high = 10.89;
LO_phonon_meV = 36.25;
piezo_kappa2 = 0.006;

U0_meV = 228;
Lz_nm = 15;
manning_prefactor = 6;
domain_nm = 60;
Nz = 1502;
n_states = 6;

Nd_sheet = 1e15;
doping_width_nm = 2;
B_T = 5;
E_kVcm = 25;
kx_inv_m = 0;
T_K = 50;

max_outer = 100;
mix = 0.5;
tol_VH_relative = 1e-3;
tol_EF_meV = 1e-3;

E0_kVcm = 4.5;
a0_nm = 7.5;
transitions = [1 2];
photon_orders = [1 2];
photon_energy_meV = 5:0.01:150;
qz_max_inv_nm = 2;
Nqz = 101;
v_max = 2;

e = 1.602176634e-19;
hbar = 1.054571817e-34;
kB = 1.380649e-23;
eps0 = 8.8541878128e-12;
m0 = 9.1093837015e-31;
meV = e*1e-3;
nm = 1e-9;

%% SCHRODINGER-POISSON TU HOP
z_m = linspace(-domain_nm/2,domain_nm/2,Nz).'*nm;
z_nm = z_m/nm;
dz_m = z_m(2)-z_m(1);
z_in = z_m(2:end-1);
u = z_m/(Lz_nm*nm);
Vconf_J = U0_meV*meV*manning_prefactor.*(sech(u).^4-sech(u).^2);

Nd_z_m3 = zeros(size(z_m));
width_m = doping_width_nm*nm;
mask = abs(z_m) <= width_m/2;
covered = trapz(z_m,double(mask));
Nd_z_m3(mask) = Nd_sheet/covered;

mstar = mstar_rel*m0;
magnetic_in = (hbar*kx_inv_m+e*B_T.*z_in).^2/(2*mstar);
electric_in = -e*E_kVcm*1e5.*z_in;
n_in = numel(z_in);
offH = (-hbar^2/(2*mstar)/dz_m^2)*ones(n_in-1,1);
offD = ones(n_in-1,1);
D2 = spdiags([[offD;0],-2*ones(n_in,1),[0;offD]],[-1 0 1],n_in,n_in);

VH_J = zeros(size(z_m));
EF_old_J = Inf;
VH_error_meV = Inf;
VH_error_relative = Inf;
EF_error_meV = Inf;

for iterations = 1:max_outer
    Veff_in = Vconf_J(2:end-1)+VH_J(2:end-1)+magnetic_in+electric_in;
    H = spdiags([[offH;0], 2*hbar^2/(2*mstar)/dz_m^2+Veff_in(:), [0;offH]], [-1 0 1],n_in,n_in);
    [Psi_in,E_J] = lowest_eigenpairs(H,n_states);
    Psi_in = bsxfun(@rdivide,Psi_in,sqrt(sum(abs(Psi_in).^2,1)*dz_m));
    Psi = zeros(numel(z_m),n_states);
    Psi(2:end-1,:) = Psi_in;

    
    [EF_J,subband_sheet_m2] = solve_fermi_level(Nd_sheet,E_J,T_K,mstar_rel,kB,m0,hbar);
    

    n_z_m3 = (abs(Psi).^2)*subband_sheet_m2(:);
    rhs = (e^2/(eps_static*eps0)).*(Nd_z_m3-n_z_m3);
    VH_new_J = zeros(size(z_m));
    VH_new_J(2:end-1) = D2\(dz_m^2.*rhs(2:end-1));
    VH_next_J = mix*VH_new_J+(1-mix)*VH_J;
    VH_error_meV = max(abs(VH_next_J-VH_J))/meV;
    VH_error_relative = max(abs(VH_next_J-VH_J))/max(abs(VH_J));
    EF_error_meV = abs(EF_J-EF_old_J)/meV;
    if VH_error_relative <= tol_VH_relative && EF_error_meV <= tol_EF_meV
        VH_J = VH_next_J;
        break
    end
    VH_J = VH_next_J;
    EF_old_J = EF_J;
end

% Giai lan cuoi voi the Hartree da hoi tu.
Veff_in = Vconf_J(2:end-1)+VH_J(2:end-1)+magnetic_in+electric_in;
H = spdiags([[offH;0], 2*hbar^2/(2*mstar)/dz_m^2+Veff_in(:), [0;offH]], [-1 0 1],n_in,n_in);
[Psi_in,E_J] = lowest_eigenpairs(H,n_states);
Psi_in = bsxfun(@rdivide,Psi_in,sqrt(sum(abs(Psi_in).^2,1)*dz_m));
Psi = zeros(numel(z_m),n_states);
Psi(2:end-1,:) = Psi_in;

[EF_J,subband_sheet_m2] = solve_fermi_level(Nd_sheet,E_J,T_K,mstar_rel,kB,m0,hbar);

n_z_m3 = (abs(Psi).^2)*subband_sheet_m2(:);
Vmag_J = (hbar*kx_inv_m+e*B_T.*z_m).^2/(2*mstar);
Velec_J = -e*E_kVcm*1e5.*z_m;
Veff_J = Vconf_J+VH_J+Vmag_J+Velec_J;

E_meV = E_J(:).'/meV;
EF_meV = EF_J/meV;
Vconf_meV = Vconf_J/meV;
VH_meV = VH_J/meV;
Veff_meV = Veff_J/meV;
charge_integral_m2 = trapz(z_m,n_z_m3);

%% TINH THANH PHAN CHUYEN MUC
qz_inv_m = linspace(0,qz_max_inv_nm,Nqz) / nm;
wqz = zeros(size(qz_inv_m));
wqz([1 end]) = [qz_inv_m(2) - qz_inv_m(1), qz_inv_m(end) - qz_inv_m(end-1)] / 2;
wqz(2:end-1) = (qz_inv_m(3:end) - qz_inv_m(1:end-2)) / 2;
phase = exp(1i * (z_m * qz_inv_m));

n_transitions = size(transitions,1);
initial_state = zeros(n_transitions,1);
deltaE_J = zeros(n_transitions,1);
dipole_m = zeros(n_transitions,1);
K_qz = zeros(n_transitions,1);
for it = 1:n_transitions
    i = transitions(it,1);
    f = transitions(it,2);
    product = conj(Psi(:,i)) .* Psi(:,f);
    Iqz = trapz(z_m,bsxfun(@times,product,phase),1);
    initial_state(it) = i;
    deltaE_J(it) = E_J(f)-E_J(i);
    dipole_m(it) = trapz(z_m,product.*z_m);
    K_qz(it) = sum(abs(Iqz).^2 .* wqz);
end
deltaE_meV = deltaE_J/meV;

%% TINH CONG SUAT HAP THU
energy_meV = photon_energy_meV(:).';
Ephot_J = energy_meV*meV;
n_energy = numel(Ephot_J);
n_orders = numel(photon_orders);
optical_emission = zeros(n_orders,n_energy);
optical_absorption = zeros(n_orders,n_energy);
piezoelectric = zeros(n_orders,n_energy);

kBT = kB*T_K;
kappa_LO = 1/(1/eps_high-1/eps_static);
hw0_J = LO_phonon_meV*meV;
N0 = 1/expm1(hw0_J/kBT);
E0_Vm = E0_kVcm*1e5;
a0_m = a0_nm*nm;
qd_inv_m = sqrt(Nd_sheet*e^2/(Lz_nm*nm*eps_static*eps0*kBT));

for iE = 1:n_energy
    Eph_J = Ephot_J(iE);
    Omega = Eph_J/hbar;
    P0_LO = E0_Vm^4*e^4*mstar^2*(hw0_J/hbar)*sqrt(eps_static) /(256*pi^3*hbar^6*eps0*kappa_LO*Omega^2);
    P0_PE = E0_Vm^4*e^4*mstar^2*piezo_kappa2*kBT*sqrt(eps_static) /(256*pi^3*hbar^7*eps_static*eps0*Omega^2);
    for it = 1:n_transitions
        eps_alpha_J = E_J(initial_state(it));
        transition_factor = abs(dipole_m(it))^2*K_qz(it);
        for io = 1:n_orders
            ell = photon_orders(io);
            order_factor = a0_m^(2*ell)/(2^(2*ell)*factorial(ell)^2);
            Pem = P0_LO*transition_factor*order_factor*(N0+1) * taylor_bessel(eps_alpha_J,EF_J, deltaE_J(it)+hw0_J-ell*Eph_J,ell,qd_inv_m, hbar,kB,T_K,mstar_rel,m0,v_max);
            Pab = P0_LO*transition_factor*order_factor*N0 * taylor_bessel(eps_alpha_J,EF_J, deltaE_J(it)-hw0_J-ell*Eph_J,ell,qd_inv_m, hbar,kB,T_K,mstar_rel,m0,v_max);
            optical_emission(io,iE) = optical_emission(io,iE)+ Pem;
            optical_absorption(io,iE) = optical_absorption(io,iE)+ Pab;
            piezoelectric(io,iE) = piezoelectric(io,iE) + P0_PE*transition_factor*order_factor * taylor_bessel(eps_alpha_J,EF_J, deltaE_J(it)-ell*Eph_J,ell,qd_inv_m, hbar,kB,T_K,mstar_rel,m0,v_max);
        end
    end
end

%% FWHM TRUC TIEP
empty_metric = struct('peak',NaN,'peak_x',NaN,'fwhm',NaN,'hwhm',NaN, ...
    'left_x',NaN,'right_x',NaN,'left_hwhm',NaN, ...
    'right_hwhm',NaN,'half_level',NaN);
metrics_optical_emission = repmat(empty_metric,n_orders,1);
metrics_optical_absorption = repmat(empty_metric,n_orders,1);
metrics_piezoelectric = repmat(empty_metric,n_orders,1);
for io = 1:n_orders
    metrics_optical_emission(io) = profile_fwhm( ...
        energy_meV,optical_emission(io,:));
    metrics_optical_absorption(io) = profile_fwhm( ...
        energy_meV,optical_absorption(io,:));
    metrics_piezoelectric(io) = profile_fwhm( ...
        energy_meV,piezoelectric(io,:));
end


%% CAC HAM TINH TOAN TOI THIEU
function [Psi,E] = lowest_eigenpairs(H,n_states)
    H = (H+H')/2;
    opts = struct('issym',true,'tol',1e-10,'maxit',4000);
    [Psi,D] = eigs(H,n_states,'smallestreal',opts);
    E = real(diag(D));
    [E,idx] = sort(E,'ascend');
    Psi = real(Psi(:,idx));
end

function [EF,Ni] = solve_fermi_level(Nd_sheet,E,T_K,mstar_rel,kB,m0,hbar)
    kBT = kB*T_K;
    A = mstar_rel*m0*kBT/(pi*hbar^2);
    ei = E(:)/kBT;
    softplus = @(x) max(x,0)+log1p(exp(-abs(x)));
    eta = fzero(@(eta) sum(softplus(eta-ei))-Nd_sheet/A, [min(ei)-50,min(ei)+Nd_sheet/A]);
    EF = eta*kBT;
    Ni = A*softplus(eta-ei);
end

function Tval = taylor_bessel(eps_alpha_J,EF_J,D_J,ell,qd,hbar,kB,T_K,mstar_rel,m0,v_max)
    kBT = kB*T_K;
    me = mstar_rel*m0;
    absD = abs(D_J);
    Tval = 0;
    for v = 1:v_max
        LambdaK = besselk(0,v*absD/(2*kBT),0);
        LambdaJ = (hbar^2*qd^2/(2*me))*LambdaK;
        Tval = Tval+(ell+1)*exp(-v*(eps_alpha_J+D_J/2-EF_J)/kBT)* qd^(2*ell-2)*LambdaJ;
    end
end

function m = profile_fwhm(x,y)
    x = x(:).';
    y = y(:).';
    m = struct('peak',NaN,'peak_x',NaN,'fwhm',NaN,'hwhm',NaN, ...
        'left_x',NaN,'right_x',NaN,'left_hwhm',NaN, ...
        'right_hwhm',NaN,'half_level',NaN);

    [m.peak,ip] = max(y);
    if m.peak <= 0, return; end
    m.peak_x = x(ip);
    m.half_level = m.peak/2;
    h = m.half_level;
    il = find(y(1:ip) <= h,1,'last');
    ir = ip-1+find(y(ip:end) <= h,1,'first');
    m.left_x = x(il);
    m.right_x = x(ir);
    m.left_hwhm = m.peak_x-m.left_x;
    m.right_hwhm = m.right_x-m.peak_x;
    m.fwhm = m.right_x-m.left_x;
    m.hwhm = m.fwhm/2;
end
