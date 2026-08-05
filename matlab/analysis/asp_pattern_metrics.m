function met = asp_pattern_metrics(f, scn, hRef)
%ASP_PATTERN_METRICS Spatial performance of a weight vector, in dB.
%
%   MET = ASP_PATTERN_METRICS(F, SCN, HREF)
%
%   F      nAnt x 1 beamformer weights
%   SCN    scenario from asp_scenario (supplies the TRUE signature vectors)
%   HREF   quiescent weight vector for the "before" reference
%          (default: ones/sqrt(N))
%
%   ALL GAINS USE ONE CONSISTENT, UNAMBIGUOUS DEFINITION
%   ----------------------------------------------------
%       G(a) = |f' a|^2 / ||f||^2                                       (1)
%
%   With a phase-only steering vector (||a||^2 = N) this has two properties
%   that make it the only sensible normalisation:
%     * G = 1 (0 dB) for a single antenna element, so dB values read
%       directly as "gain over one antenna";
%     * E[G] = 1 over directions with E[a a'] = I, so 0 dB is also the
%       average gain over the sky.  A beamformer cannot have positive
%       average gain; it can only redistribute.
%
%   The original scripts used
%       NullDepth_dB = 10*log10(|f'*b|^2 / |h'*b|^2)
%   which measures the null relative to whatever gain the arbitrary
%   quiescent vector h happened to have toward the spoofer.  If h is
%   unlucky and already has 6 dB of loss toward the spoofer, the reported
%   null depth is 6 dB better than the truth; if h is lucky, 6 dB worse.
%   The quantity is not a property of the beamformer at all.  Both are
%   reported below so the numbers can be compared to the original scripts,
%   but MET.nullGainDB is the one to specify a product against.

n = numel(f);
if nargin < 3 || isempty(hRef)
    hRef = ones(n,1)/sqrt(n);
end

fPow = real(f'*f);
hPow = real(hRef'*hRef);

gainDb = @(w, a, wp) 10*log10(max(abs(w'*a)^2, realmin) / max(wp, realmin));

isAuth  = strcmp(scn.kinds, 'auth');
isSpoof = strcmp(scn.kinds, 'spoof');
isMp    = strcmp(scn.kinds, 'spoofmp');

met.nAnt = n;

% --- spoofer
if any(isSpoof)
    b = scn.bTrue;
    met.nullGainDB      = gainDb(f, b, fPow);
    met.quiescentGainDB = gainDb(hRef, b, hPow);
    met.suppressionDB   = met.nullGainDB - met.quiescentGainDB;
else
    met.nullGainDB = NaN; met.quiescentGainDB = NaN; met.suppressionDB = NaN;
end

% --- spoofer ground bounce
if any(isMp)
    bm = scn.aTrue(:, find(isMp,1));
    met.mpGainDB = gainDb(f, bm, fPow);
else
    met.mpGainDB = NaN;
end

% --- authentic satellites
idxAuth = find(isAuth);
met.authGainDB     = zeros(1, numel(idxAuth));
met.authGainQuiDB  = zeros(1, numel(idxAuth));
for k = 1:numel(idxAuth)
    met.authGainDB(k)    = gainDb(f,    scn.aTrue(:,idxAuth(k)), fPow);
    met.authGainQuiDB(k) = gainDb(hRef, scn.aTrue(:,idxAuth(k)), hPow);
end
met.authGainMeanDB  = 10*log10(mean(10.^(met.authGainDB/10)));
met.authGainMinDB   = min(met.authGainDB);
met.authGainMedianDB = median(met.authGainDB);

% Gain relative to the quiescent (non-adaptive) beam.  This is the number
% that answers "what did nulling cost my satellites", and it is the correct
% comparison because the quiescent beam is what the product would output if
% the detector declared no threat.  The ABSOLUTE gain figures depend on the
% quiescent beam's own shape - h = ones/sqrt(N) is a zenith-pointing beam
% with +10*log10(N) dB at zenith and correspondingly less elsewhere - so
% absolute numbers alone are easy to misread.
met.authGainQuiMeanDB = 10*log10(mean(10.^(met.authGainQuiDB/10)));
met.authGainChangeDB  = met.authGainDB - met.authGainQuiDB;
met.authGainChangeMeanDB = met.authGainMeanDB - met.authGainQuiMeanDB;
met.authGainChangeWorstDB = min(met.authGainChangeDB);

% Number of satellites pushed below single-antenna performance.  This is
% the metric that decides whether the receiver still has a usable geometry,
% and it is the metric that most cleanly separates 3 elements from 4.
met.nAuthBelowUnity = sum(met.authGainDB < 0);
met.nAuthBelowMinus6 = sum(met.authGainDB < -6);

% Post-mitigation spoof-to-authentic ratio.  This, not the null depth, is
% what determines whether the receiver acquires the right peak.
if any(isSpoof)
    saprOutDB = met.nullGainDB - met.authGainMeanDB + scn.cfg.saprDB;
    met.saprOutDB = saprOutDB;
    met.saprImprovementDB = scn.cfg.saprDB - saprOutDB;
else
    met.saprOutDB = NaN;
    met.saprImprovementDB = NaN;
end

end
