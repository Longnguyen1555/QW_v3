function m = profile_fwhm(x, y)
%PROFILE_FWHM Width of a resolved peak at half its maximum.

    x = x(:).'; y = y(:).';
    m = struct('peak',NaN,'peak_x',NaN,'fwhm',NaN,'hwhm',NaN, ...
        'left_x',NaN,'right_x',NaN,'left_hwhm',NaN, ...
        'right_hwhm',NaN,'half_level',NaN,'resolved',false);

    if numel(x) ~= numel(y)
        error('QW:FWHM:SizeMismatch','x and y must have the same size.');
    end
    if numel(x) < 3 || ~isreal(x) || any(~isfinite(x)) || any(diff(x) <= 0)
        error('QW:FWHM:InvalidGrid','x must be finite and strictly increasing.');
    end
    if ~isreal(y) || any(~isfinite(y)) || any(y < 0), return; end

    [m.peak, ip] = max(y);
    if m.peak <= 0, return; end
    m.peak_x = x(ip); m.half_level = m.peak/2;

    % Without sampled half-height shoulders, the width is grid-limited.
    h = m.half_level;
    if ip == 1 || ip == numel(y) || y(ip-1) < h || y(ip+1) < h
        return;
    end
    il = find(y(1:ip) <= h,1,'last');
    ir = ip-1 + find(y(ip:end) <= h,1,'first');
    if isempty(il) || isempty(ir) || il == ip || ir == ip, return; end

    m.left_x = crossing(x,y,il,h);
    m.right_x = crossing(x,y,ir-1,h);
    m.left_hwhm = m.peak_x-m.left_x;
    m.right_hwhm = m.right_x-m.peak_x;
    m.fwhm = m.right_x-m.left_x;
    m.hwhm = m.fwhm/2;
    m.resolved = true;
end

function xc = crossing(x,y,i,h)
    if y(i+1) == y(i)
        xc = (x(i)+x(i+1))/2;
    else
        xc = x(i)+(h-y(i))*(x(i+1)-x(i))/(y(i+1)-y(i));
    end
end
