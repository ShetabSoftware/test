function [spacing, info] = asp_geometry_spacing_limit(name, lambda, radomeRadiusLambda, gratingTargetDB)
%ASP_GEOMETRY_SPACING_LIMIT Largest usable spacing for a geometry.
%
%   [SPACING, INFO] = ASP_GEOMETRY_SPACING_LIMIT(NAME, LAMBDA, RADOMER, TARGETDB)
%
%   Returns the largest spacing parameter (side length or ring radius,
%   depending on the geometry) that satisfies BOTH constraints that bound a
%   real CRPA:
%
%     1. the array fits inside a radome of radius RADOMER wavelengths;
%     2. the peak grating response over the visible region stays at or below
%        TARGETDB (default -3 dB), so that a null steered at the spoofer
%        does not simultaneously land on a second, distinct arrival
%        direction.
%
%   Sizing every candidate geometry this way is the only fair basis for an
%   element-count comparison.  Comparing a 3-element array and a 4-element
%   array at "the same lambda/2 spacing" compares two arrays with different
%   apertures, different null widths and different ambiguity margins, and
%   attributes the whole difference to the element count.

if nargin < 3 || isempty(radomeRadiusLambda), radomeRadiusLambda = 0.5;  end
if nargin < 4 || isempty(gratingTargetDB),    gratingTargetDB    = -3.0; end

lo = 0.05;
hi = 1.20;

feasible = @(s) checkFeasible(name, lambda, s, radomeRadiusLambda, gratingTargetDB);

if ~feasible(lo)
    error('asp_geometry_spacing_limit:infeasible', ...
        'Geometry "%s" is infeasible even at the smallest spacing.', name);
end

for it = 1:40
    mid = 0.5*(lo+hi);
    if feasible(mid)
        lo = mid;
    else
        hi = mid;
    end
end

spacing = lo;
[pos, meta] = asp_array_geometry(name, lambda, spacing);
info.meta = meta;
info.pos = pos;
info.spacing = spacing;
info.radomeRadiusLambda = radomeRadiusLambda;
info.gratingTargetDB = gratingTargetDB;
info.limitedBy = 'grating';
if meta.apertureLambda >= radomeRadiusLambda - 1e-3
    info.limitedBy = 'radome';
end

end

% -------------------------------------------------------------------------
function tf = checkFeasible(name, lambda, s, radomeR, targetDB)
pos = asp_array_geometry(name, lambda, s);
aperture = max(sqrt(sum(pos.^2,1)))/lambda;
if aperture > radomeR + 1e-9
    tf = false;
    return;
end
amb = asp_array_ambiguity(pos, lambda, 201);
tf = amb.peakLevelDB <= targetDB;
end
