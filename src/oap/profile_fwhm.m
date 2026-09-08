function m = profile_fwhm(x, y)
%PROFILE_FWHM FWHM of the dominant interior local peak.
%
% The routine leaves the spectrum unchanged, ignores non-finite samples,
% and uses the nearest half-height crossing on each side of the selected
% local peak.  A boundary maximum is reported as unresolved.

    m = empty_metrics();
    if numel(x) ~= numel(y)
        return;
    end

    x = x(:).';
    y = real(y(:).');
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);

    if numel(x) < 3 || any(diff(x) <= 0)
        return;
    end

    [global_peak, global_idx] = max(y);
    m.peak = global_peak;
    m.peak_x = x(global_idx);
    if ~isfinite(global_peak) || global_peak <= 0
        return;
    end

    if y(1) == global_peak || y(end) == global_peak
        m.is_boundary_peak = true;
        return;
    end

    local_peak = false(size(y));
    local_peak(2:end-1) = ...
        y(2:end-1) >= y(1:end-2) & ...
        y(2:end-1) >= y(3:end) & ...
        (y(2:end-1) > y(1:end-2) | y(2:end-1) > y(3:end));
    candidates = find(local_peak);
    if isempty(candidates)
        return;
    end

    [peak, selected_rel] = max(y(candidates));
    idx = candidates(selected_rel);
    m.peak = peak;
    m.peak_x = x(idx);
    m.half_max = peak/2;

    m.left_x = nearest_left_crossing(x, y, idx, m.half_max, local_peak);
    m.right_x = nearest_right_crossing(x, y, idx, m.half_max, local_peak);
    if isfinite(m.left_x) && isfinite(m.right_x)
        m.fwhm = m.right_x - m.left_x;
        m.hwhm = m.fwhm/2;
    end
end

function m = empty_metrics()
    m = struct('peak', NaN, 'peak_x', NaN, 'half_max', NaN, ...
        'left_x', NaN, 'right_x', NaN, 'fwhm', NaN, 'hwhm', NaN, ...
        'is_boundary_peak', false);
end

function xc = nearest_left_crossing(x, y, idx, half, local_peak)
    xc = NaN;
    for j = idx-1:-1:1
        if y(j) <= half
            xc = linear_crossing(x(j), y(j), x(j+1), y(j+1), half);
            return;
        end
        if local_peak(j)
            return;
        end
    end
end

function xc = nearest_right_crossing(x, y, idx, half, local_peak)
    xc = NaN;
    for j = idx+1:numel(y)
        if y(j) <= half
            xc = linear_crossing(x(j-1), y(j-1), x(j), y(j), half);
            return;
        end
        if local_peak(j)
            return;
        end
    end
end

function xc = linear_crossing(x1,y1,x2,y2,target)
    if y2 == y1
        xc = (x1+x2)/2;
    else
        xc = x1 + (target-y1)*(x2-x1)/(y2-y1);
    end
end
