function m = profile_fwhm(x, y)
%PROFILE_FWHM FWHM of the dominant peak using linear interpolation.

    x = x(:).';
    y = real(y(:).');
    y(~isfinite(y)) = 0;

    [peak, idx] = max(y);
    m = struct('peak',peak,'peak_x',NaN,'fwhm',NaN,'hwhm',NaN, ...
               'left_x',NaN,'right_x',NaN);

    if isempty(idx) || peak <= 0
        return;
    end

    m.peak_x = x(idx);
    half = peak/2;

    il = find(y(1:idx) <= half, 1, 'last');
    if isempty(il) || il == idx
        left_x = NaN;
    else
        left_x = linear_crossing(x(il), y(il), ...
                                 x(il+1), y(il+1), half);
    end

    ir_rel = find(y(idx:end) <= half, 1, 'first');
    if isempty(ir_rel) || ir_rel == 1
        right_x = NaN;
    else
        ir = idx + ir_rel - 1;
        right_x = linear_crossing(x(ir-1), y(ir-1), ...
                                  x(ir), y(ir), half);
    end

    m.left_x = left_x;
    m.right_x = right_x;
    if isfinite(left_x) && isfinite(right_x)
        m.fwhm = right_x-left_x;
        m.hwhm = m.fwhm/2;
    end
end

function xc = linear_crossing(x1,y1,x2,y2,target)
    if y2 == y1
        xc = (x1+x2)/2;
    else
        xc = x1 + (target-y1)*(x2-x1)/(y2-y1);
    end
end
