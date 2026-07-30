function n_z = electron_density_from_subbands(Psi, Ni)
%ELECTRON_DENSITY_FROM_SUBBANDS n(z)=sum_i Ni |psi_i(z)|^2.
    n_z = (abs(Psi).^2) * Ni(:);
end
