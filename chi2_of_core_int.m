function [R2, m, Gc, ConstraintMisfit, sigma_m] = chi2_of_core_int(ConstraintVal, G, core, M)

Cvals_core = core.sigma_depth(:).^2; % depth uncertainties squared (estimated using the age-uncertainties

% values for the uniformity constraint
ConstraintVector = ConstraintVal./abs(diff([core.age; 0])); % divided by the difference between the # years between tie points
Gc               = speye(length(core.age));  % sparse Identity matrix
d                = [core.depth; ones(length(ConstraintVector),1)];

G0      = [G; Gc];
C       = spdiags([Cvals_core; ConstraintVector(:)], 0, size(G0, 1), size(G0, 1));
T       = chol(C); % chol decomposition of a diagonal matrix produces an "upper triangle" that is still diagonal
m       = (T\G0)\(T\d);
r       = core.depth - G*m;
C_core  = spdiags(Cvals_core, 0, size(G, 1), size(G, 1));
chi2    = r'*(C_core\r);

R2 = (chi2-(length(core.depth)-2))^2;

ConstraintMisfit = sum(abs(diff([core.age(:); 0])).*(m-1).^2); % eq 10, 2nd term without lambda

if nargout == 5
    ginv    = (G0'*(C\G0))\G0';
    Cm      = ginv*(C\ginv');
    sigma_m = sqrt(diag(Cm));
end


