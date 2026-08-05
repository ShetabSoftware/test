function [a, tauOffset] = asp_steering(azDeg, elDeg, pos, lambda, c)
%ASP_STEERING Narrowband array response and the matching propagation delays.
%
%   [A, TAUOFFSET] = ASP_STEERING(AZDEG, ELDEG, POS, LAMBDA, C)
%
%   AZDEG, ELDEG  source direction(s), degrees.  Azimuth is measured from
%                 +x toward +y; elevation from the xy-plane toward +z.
%                 May be vectors of equal length, giving A as N x nDir.
%   POS           3 x N element positions in metres
%   LAMBDA        carrier wavelength in metres
%   C             speed of light (only needed for TAUOFFSET)
%
%   A          N x nDir narrowband array response
%   TAUOFFSET  N x nDir extra propagation delay, in seconds, of each element
%              relative to the coordinate origin
%
%   SIGN CONVENTION - stated once, used everywhere
%   ----------------------------------------------
%   Let u be the unit vector pointing FROM the array TOWARD the source.  An
%   element displaced toward the source sees the wavefront EARLIER, so its
%   propagation delay relative to the origin is
%
%       tau_i = -(u . p_i)/c                                            (1)
%
%   and the complex baseband signal at element i is
%
%       r_i(t) = s(t - tau_i) * exp(-j*2*pi*fc*tau_i)
%              = s(t + u.p_i/c) * exp(+j*2*pi*(u.p_i)/lambda)           (2)
%
%   so the narrowband array response is
%
%       a_i = exp(+j*2*pi*(u . p_i)/lambda)                             (3)
%
%   TAUOFFSET returns tau_i from (1).  The signal generator applies BOTH (1)
%   and (3), which is what separates a wideband-correct model from a
%   narrowband one.
%
%   WHY THE DELAY TERM MATTERS
%   --------------------------
%   The original scripts applied only (3).  Dropping (1) makes the array
%   perfectly narrowband by construction, so the simulation can report null
%   depths that hardware cannot reach.  The physical floor is
%
%       null depth >= 20*log10(2*pi*(B/sqrt(12))*std(tau))              (4)
%
%   For a 4-element circ4 array at L1 (D = 0.27 m, std(tau) ~ 0.25 ns) and
%   B = 16.368 MHz this is about -46 dB, which is below the estimator's own
%   accuracy and therefore not the binding constraint - but the SAME formula
%   applied to inter-channel group-delay mismatch IS binding: 1 ns of
%   uncalibrated cable/filter skew alone limits the null to about -34 dB.

if nargin < 5 || isempty(c)
    c = 299792458;
end

azDeg = azDeg(:).';
elDeg = elDeg(:).';
if numel(azDeg) ~= numel(elDeg)
    error('asp_steering:dim', 'azDeg and elDeg must have the same length.');
end

ce = cosd(elDeg);
u  = [ce .* cosd(azDeg); ce .* sind(azDeg); sind(elDeg)];   % 3 x nDir

proj = pos.' * u;                 % N x nDir, equals (u . p_i)

a         = exp(1i * 2*pi * proj / lambda);
tauOffset = -proj / c;

end
