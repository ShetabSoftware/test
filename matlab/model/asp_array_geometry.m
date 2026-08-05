function [pos, meta] = asp_array_geometry(name, lambda, spacingLambda)
%ASP_ARRAY_GEOMETRY Element phase-centre positions for named array layouts.
%
%   [POS, META] = ASP_ARRAY_GEOMETRY(NAME, LAMBDA, SPACINGLAMBDA)
%
%   POS is 3 x N in metres, columns are element phase centres, z = 0 (planar
%   array in the local horizontal plane).  SPACINGLAMBDA is the nominal
%   inter-element spacing in wavelengths for the linear/regular layouts and
%   the ring RADIUS in wavelengths for the circular layouts.
%
%   Supported NAME values:
%     'tri3'   equilateral triangle, side = spacing.  The paper's array.
%     'sq4'    2 x 2 square, side = spacing.
%     'circ4'  4 elements on a circle of radius = spacing, at 45/135/225/315
%     'y4'     equilateral triangle on a circle of radius = spacing, plus a
%              centre element
%     'circ7'  6 on a ring + 1 centre (classic 7-element CRPA)
%     'circ8'  8 on a ring
%     'ula2'   2-element baseline, for sanity checks
%
%   META reports the quantities that actually decide array performance:
%     .maxBaselineLambda  largest inter-element distance in wavelengths.
%                         Sets the angular width of a spatial null; the
%                         null half-width is roughly lambda/(2*D).
%     .apertureLambda     radius of the smallest enclosing circle, i.e. the
%                         radome size you have to pay for.
%     .nUniqueBaselines   number of distinct baseline LENGTHS.  Redundant
%                         baselines waste aperture for DOA estimation but
%                         cost nothing for the projection-based nulling in
%                         this design, which never uses the manifold.
%     .isotropy           ||sum(p*p')/ (sum of squares) - I/2||, a measure of
%                         how azimuth-independent the beam/null width is.
%                         0 means the second moment of the aperture is
%                         circularly symmetric.
%     .ambiguous          true if any baseline exceeds lambda/2, in which
%                         case the array manifold is not injective over the
%                         visible region and a null steered at the spoofer
%                         also lands on a mirror direction that may contain
%                         a satellite.
%
%   ON THE AMBIGUITY FLAG
%   ---------------------
%   For a planar array the visible region is the full unit disk of direction
%   cosines (u,v), u^2+v^2 <= 1, because a signal from elevation E and
%   azimuth A maps to (cos E cos A, cos E sin A).  The manifold element
%   exp(j*2*pi*(x*u + y*v)/lambda) is therefore unambiguous over the visible
%   region only if every |x| and |y| separation stays within lambda/2.  This
%   caps a 4-element square at 9.5 cm sides at L1 and is why commercial
%   CRPAs are physically the size they are.  Exceeding it is sometimes done
%   deliberately (a "thinned" array) but only with a manifold-aware
%   algorithm that can resolve the ambiguity, which the projection method
%   here cannot.

if nargin < 3 || isempty(spacingLambda)
    spacingLambda = 0.5;
end

d = spacingLambda * lambda;

switch lower(name)
    case 'ula2'
        pos = [0 d; 0 0; 0 0];

    case 'tri3'
        % Equilateral triangle of side d, centred on its centroid.
        R = d/sqrt(3);
        ang = [90 210 330] * pi/180;
        pos = [R*cos(ang); R*sin(ang); zeros(1,3)];

    case 'sq4'
        pos = [0 d 0 d; 0 0 d d; 0 0 0 0];
        pos = pos - mean(pos,2);

    case 'circ4'
        ang = (45:90:315) * pi/180;
        pos = [d*cos(ang); d*sin(ang); zeros(1,4)];

    case 'y4'
        ang = [90 210 330] * pi/180;
        pos = [ [d*cos(ang); d*sin(ang); zeros(1,3)], [0;0;0] ];

    case 'circ7'
        ang = (0:60:300) * pi/180;
        pos = [ [d*cos(ang); d*sin(ang); zeros(1,6)], [0;0;0] ];

    case 'circ8'
        ang = (0:45:315) * pi/180;
        pos = [d*cos(ang); d*sin(ang); zeros(1,8)];

    otherwise
        error('asp_array_geometry:name', 'Unknown geometry "%s".', name);
end

n = size(pos,2);

% --- metadata
dist = zeros(n);
for i = 1:n
    for j = 1:n
        dist(i,j) = norm(pos(:,i) - pos(:,j));
    end
end
offDiag = dist(~eye(n));

meta.nElements        = n;
meta.positions        = pos;
meta.maxBaselineLambda = max(offDiag) / lambda;
meta.minBaselineLambda = min(offDiag) / lambda;
meta.apertureLambda   = max(sqrt(sum(pos.^2,1))) / lambda;
meta.nUniqueBaselines = numel(uniquetol_local(offDiag/lambda, 1e-6));
meta.ambiguous        = meta.maxBaselineLambda > 0.5 + 1e-9;

% Second moment isotropy: for a planar aperture, sum_i p_i p_i^T should be
% proportional to the 2x2 identity for the null/beam width to be independent
% of azimuth to second order.
M = pos(1:2,:) * pos(1:2,:)';
if trace(M) > 0
    meta.isotropy = norm(M/trace(M) - eye(2)/2, 'fro');
else
    meta.isotropy = 0;
end

end

% -------------------------------------------------------------------------
function u = uniquetol_local(v, tol)
% uniquetol is not available in all Octave builds.
v = sort(v(:));
u = v(1);
for k = 2:numel(v)
    if abs(v(k) - u(end)) > tol
        u(end+1,1) = v(k); %#ok<AGROW>
    end
end
end
