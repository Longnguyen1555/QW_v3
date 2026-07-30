function Nd_z = delta_doping_profile(z, Nd_sheet, width)
%DELTA_DOPING_PROFILE Finite-width representation of a centered delta layer.
% The returned volume density satisfies integral Nd(z) dz = Nd_sheet.

    Nd_z = zeros(size(z));
    if Nd_sheet <= 0
        return;
    end

    if width <= 0
        [~, idx] = min(abs(z));
        dz = z(2) - z(1);
        Nd_z(idx) = Nd_sheet / dz;
        return;
    end

    mask = abs(z) <= width / 2;
    covered = trapz(z, double(mask));
    if covered <= 0
        [~, idx] = min(abs(z));
        dz = z(2) - z(1);
        Nd_z(idx) = Nd_sheet / dz;
    else
        Nd_z(mask) = Nd_sheet / covered;
    end
end
