%% ========================================================================
function varargout = est_timescale_u(varargin)
    % functions to estimate the timescale for an ice core.  Uses a model in
    % which the thinning rate is proportional to the horizontal velocity and
    % the glacier is frozen to its bed. 
    [varargout{1:nargout}] = feval(varargin{:});
end

%% ========================================================================
function T = age_for_ref_model(z, H, bdot_ref) % age for reference model code
    % calculate the age of a set of points at height z, for an ice cap 
    % of thickness H and accumulation bdot_ref
    slowness = @(x) 1./vertical_velocity(x, H, bdot_ref); % define model
    T1 = NaN(size(z)); % preallocate
    [zs, ind] = sort(z, 'descend');
    if(zs(1)) == H % Boundary condition; if deepest point in depth = ice cap height, velocity = 0
        T1(1) = 0;
    else
        T1(1) = integral(slowness, H, zs(1));
    end
    for k = 2:length(z)
        T1(k) = T1(k-1) + integral(slowness, zs(k-1), zs(k));
    end
    T = zeros(size(z));
    T(ind) = T1;
end

%% ========================================================================
function w = vertical_velocity(z, H, bdot_ref)
    n = 3;
    edot = -bdot_ref/(H.*(1 - 1/(n + 2)));  % or, -bdot_ref/H * (n+1)/(n+2)
    mu = (1 - z/H);
    w = H * edot*(1 - (1/(n+2)) - (mu - (1/(n+2))*mu.^(n + 2)));
end
