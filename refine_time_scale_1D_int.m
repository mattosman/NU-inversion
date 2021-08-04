function [ConstraintMisfit, R2, m_full]=refine_time_scale_1D_int(TS, Bdot0, core)
% TS = structure containing the idealized an idealized age-depth scale for
%   an ice cap with reference accumulation (Bdot0, below) of height max(TS.depth)
% Bdot0 = the reference accumulation;
% core = structure containing the tie-point constaints' depths, ages, and uncertainties;

% 1-D version of the script refine_time_scale from the FE_solution directory.
T = TS.age;
depth = TS.depth(:);
if ~isfield(TS, 'accum')
    TS.accum = 1; 
end

% this is the depth ratio for the assumed accumulation rate
% as a function of time
depthRatio = diff(depth(:)) ./ diff(TS.accum*T(:)); % this is fraction of the original annual thickness at this depth

% for all flow bands, calculate the centerline parameters (Mc)
Tc = (T(1:end-1) + T(2:end))/2; % average age between delta_z
Tc = [0; Tc(:)];

Mc = struct('T', Tc, 'depthRatio', [1; depthRatio],'depth', [0; (depth(1:end-1)+depth(2:end))/2]); % average depth
Mc = index_struct(Mc, isfinite(Tc));

% For a set of ages (in the core), calculate the depths by interpolation into Mc:
dt = 1;
M.age = (0:dt:max(core.age))';  % 0 to max core age point
M.depth = interp1(TS.age(isfinite(TS.age)), TS.depth(isfinite(TS.age)), M.age); % interpolates the annual layers to a depth corresponding onto the idealized depth-scale
M = index_struct(M, isfinite(M.depth)); 

% calculate the initial layer thickness
M.L0 = zeros(size(M.age)) + dt*Bdot0; % assumed to be reference accumulation
% calculate the final layer thickness
M.L1 = M.L0.*interp1(Mc.T, Mc.depthRatio, M.age); % 
 
% now try to fit to the core
% model parameters are scale values that multiply bdot0*depthRatio
core = index_struct(core, core.age ~= 0);
dL1  = zeros(1, length(core.age));
core_ind = zeros(length(M.age), 1);

these = M.age(:) < core.age(end); % find all annual layers more recent than 1st tie point
core_ind(these) = length(core.age); % identify all years with corresponding closest tie point number
dL1(end) = sum(these.*M.L1(:)); % thickness up to first tie point 
for k = 0:length(core.age)-2
    these = ((M.age(:) >= core.age(end-k)) & M.age(:) < core.age(end-k-1));
    core_ind(these) = length(core.age)-k-1; % keeps index of which age constraint it is (1 to 55)
    dL1(end-k-1) = sum( ((M.age(:) >= core.age(end-k)) & M.age(:) < core.age(end-k-1)) .* M.L1(:) ); % thickness between tie points (ice eqiv)
end
core_ind(M.age == core.age(1)) = 1;

G = triu(repmat(dL1, length(core.age), 1));
% M contains the annually-interpolated idealized age-depth scale down to lowest (oldest) assigned tie point
ConstraintVal = fminbnd(@(x)chi2_of_core_int(x, G, core, M), 10^-3, 10^3); % fminbnd - minimizing lambda wrt R2 - this leaves the determind R2 value as the new lamda
[R2, m_best, Gc, ConstraintMisfit, sigma_m] = chi2_of_core_int(ConstraintVal, G, core, M); % constrain bdot value
m_full.scale = m_best(core_ind);
m_full.sigma = sigma_m(core_ind);
m_full.age   = M.age;
m_full.depth = cumsum([0; M.L1(1:end-1).*m_full.scale(1:end-1)]);

