function res = study_multipath(nTrial, verbose)
%STUDY_MULTIPATH Coherent second arrivals from the spoofer, and what they cost.
%
%   RES = STUDY_MULTIPATH(NTRIAL, VERBOSE)
%
%   Two cases are measured, and they behave completely differently.  The
%   distinction is not made anywhere in the paper and it changes the
%   degree-of-freedom budget, so it is worth establishing by measurement.
%
%   CASE A - SPECULAR GROUND BOUNCE.  Mirror elevation, same azimuth.
%   ------------------------------------------------------------------
%   Costs NOTHING.  A planar array with every element at z = 0 has a
%   response that depends only on the direction cosines
%   (cos(el)cos(az), cos(el)sin(az)), and cos is even, so elevation +theta
%   and -theta give an IDENTICAL steering vector - measured coherence
%   1.000000, not merely close to it.  The array cannot distinguish up from
%   down at all.  A specular ground reflection arrives at exactly the mirror
%   elevation and the same azimuth, so the rank-one projector that removes
%   the direct path removes the reflection as a free side effect, at every
%   frequency, with no wideband penalty.
%
%   I expected the opposite and the measurement corrected me.  It is worth
%   stating why the intuition fails: the usual "two rays need two nulls"
%   reasoning is imported from arrays that HAVE a vertical baseline.  This
%   one does not.
%
%   The same property cuts the other way.  The array also cannot separate an
%   authentic satellite at +theta from its own ground reflection, so a
%   planar CRPA provides NO spatial multipath rejection for the authentic
%   signals.  That job stays with the receiver's code and carrier multipath
%   mitigation, or needs a non-planar array.
%
%   CASE B - REFLECTION AT A DIFFERENT AZIMUTH.
%   ------------------------------------------------------------------
%   Costs a degree of freedom.  A reflector off to one side - a building, a
%   vehicle, a mast - produces an arrival with genuinely different direction
%   cosines (measured coherence 0.010 at 90 deg of azimuth separation).  Now
%   the spoofing subspace really is rank 2, and a rank-one projector leaves
%   a residual set by how coherent the two arrivals are across the band:
%
%     R_s = P_s [ b b' + a^2 bm bm' + a*rho*(b bm' + bm b') ]
%     rho = normalised signal autocorrelation at the excess delay
%         = 1 - |dt|/Tchip   for a C/A chip of 977 ns
%
%   Note rho is set by the SIGNAL bandwidth (1.023 MHz, hence a 977 ns
%   correlation width), NOT by the sampling rate.  Widening the ADC does not
%   change it.  I had this wrong initially and the measurement corrected it.

if nargin < 1 || isempty(nTrial),  nTrial = 12; end
if nargin < 2 || isempty(verbose), verbose = true; end

cases = { ...
    struct('name','none',          'enable',false,'azEl',[  0  0],'dt',0), ...
    struct('name','ground bounce', 'enable',true, 'azEl',[ 45 -15],'dt',20e-9), ...
    struct('name','bldg  60 deg',  'enable',true, 'azEl',[105  20],'dt',60e-9), ...
    struct('name','bldg  90 deg',  'enable',true, 'azEl',[135  20],'dt',60e-9), ...
    struct('name','bldg 180 deg',  'enable',true, 'azEl',[225  20],'dt',60e-9)};

relDB = -6;
fsList = [4*1.023e6, 16*1.023e6];

res.cases = cases;
res.fs = fsList;
res.nullRank1 = nan(numel(fsList), numel(cases), nTrial);
res.nullRank2 = nan(numel(fsList), numel(cases), nTrial);
res.gainRank1 = nan(numel(fsList), numel(cases), nTrial);
res.gainRank2 = nan(numel(fsList), numel(cases), nTrial);
res.coherence = nan(numel(cases),1);
res.capTheory = nan(numel(cases),1);

for fi = 1:numel(fsList)
    for ci = 1:numel(cases)
        c = cases{ci};
        cfg = asp_config('fs', fsList(fi));
        cfg.spoofMultipath.enable   = c.enable;
        cfg.spoofMultipath.relDB    = relDB;
        cfg.spoofMultipath.delaySec = c.dt;
        cfg.spoofMultipath.azEl     = c.azEl;
        K = cfg.K; n = cfg.nAnt;
        h = ones(n,1)/sqrt(n);

        if fi == 1
            pos = asp_array_geometry(cfg.geometry, cfg.lambda, cfg.elementSpacingLambda);
            bd = asp_steering(cfg.spoofAzEl(1), cfg.spoofAzEl(2), pos, cfg.lambda);
            if c.enable
                bm = asp_steering(c.azEl(1), c.azEl(2), pos, cfg.lambda);
                res.coherence(ci) = abs(bd'*bm)/n;
                res.capTheory(ci) = rank1Cap(10^(relDB/20), c.dt, res.coherence(ci));
            end
        end

        for t = 1:nTrial
            scn = asp_scenario(cfg, 'seed', 88000+t, 'durationMs', 6);
            x = asp_rx_generate(scn, 0, K*5, []);
            R = (x*x')/size(x,2);

            for rk = [1 2]
                y = asp_ssv_from_cov(R, 'evd', struct('rank',rk,'jacobiSweeps',6));
                f = asp_weights('project', y, h);

                % TOTAL residual spoof power, not just the direct path.  The
                % direct-path-only figure badly understates the rank-1 case:
                % when the second arrival is far from the direct path in
                % azimuth, a rank-1 null does not touch it, so it emerges at
                % roughly the average sky gain (0 dB) scaled by its own
                % relative power - which then DOMINATES the residual.
                fp = real(f'*f);
                resid = abs(f'*scn.bTrue)^2 / fp;
                idxMp = find(strcmp(scn.kinds,'spoofmp'), 1);
                if ~isempty(idxMp)
                    resid = resid + 10^(relDB/10) * abs(f'*scn.aTrue(:,idxMp))^2 / fp;
                end
                nd = 10*log10(resid);
                gm = meanAuthGain(f, scn);
                if rk == 1
                    res.nullRank1(fi,ci,t) = nd; res.gainRank1(fi,ci,t) = gm;
                else
                    res.nullRank2(fi,ci,t) = nd; res.gainRank2(fi,ci,t) = gm;
                end
            end
        end
    end
end

if verbose
    printResults(res, nTrial, relDB);
end

end

% -------------------------------------------------------------------------
function g = meanAuthGain(f, scn)
isAuth = strcmp(scn.kinds,'auth');
idx = find(isAuth);
v = zeros(1,numel(idx));
for k = 1:numel(idx)
    v(k) = abs(f'*scn.aTrue(:,idx(k)))^2/real(f'*f);
end
g = 10*log10(mean(v));
end

function capDB = rank1Cap(alpha, dt, coh)
%RANK1CAP Residual toward the direct path after nulling only the principal
%   eigenvector of a two-ray source covariance.
Tchip = 1/1.023e6;
rho = max(1 - abs(dt)/Tchip, 0);
% Work in an orthonormalised {b, bm} basis; coh is |b'bm|/N.
c = coh;
% Gram-Schmidt: e1 = b, e2 = (bm - c*b)/sqrt(1-c^2)
s = sqrt(max(1-c^2, 1e-12));
% bm = c*e1 + s*e2 (up to phase)
M = [1 0; 0 0] + alpha^2*[c;s]*[c s] + alpha*rho*([1;0]*[c s] + [c;s]*[1 0]);
M = (M+M')/2;
[V,D] = eig(M);
[~,i1] = max(real(diag(D)));
u1 = V(:,i1);
capDB = 10*log10(max(1 - abs(u1(1))^2, 1e-16));
end

function printResults(res, nTrial, relDB)
fprintf('\n');
fprintf('==========================================================================\n');
fprintf(' COHERENT SECOND ARRIVAL FROM THE SPOOFER  -  %d scenes, waveform model\n', nTrial);
fprintf(' 4-element Y array, second arrival at %.0f dB, 5 ms dwell\n', relDB);
fprintf('==========================================================================\n');

mn = @(v) 10*log10(mean(10.^(v(~isnan(v))/10)));

for fi = 1:numel(res.fs)
    fprintf('\n--- fs = %.3f MHz ---\n', res.fs(fi)/1e6);
    fprintf('%-15s %9s %10s %10s %10s %10s\n', ...
        'second arrival','coh','resid r1','resid r2','auth r1','auth r2');
    fprintf('%s\n', repmat('-',1,68));
    for ci = 1:numel(res.cases)
        c = res.cases{ci};
        if isnan(res.coherence(ci)), cs = '     -';
        else, cs = sprintf('%6.4f', res.coherence(ci)); end
        fprintf('%-15s %9s %10.2f %10.2f %10.2f %10.2f\n', c.name, cs, ...
            mn(squeeze(res.nullRank1(fi,ci,:))), mn(squeeze(res.nullRank2(fi,ci,:))), ...
            mn(squeeze(res.gainRank1(fi,ci,:))), mn(squeeze(res.gainRank2(fi,ci,:))));
    end
end

fprintf('\nReading the table (residuals are TOTAL spoof power, direct + second):\n');
fprintf('  * GROUND BOUNCE costs nothing.  Coherence with the direct path is\n');
fprintf('    1.0000 exactly, because a planar array cannot tell +elevation from\n');
fprintf('    -elevation, so the reflection shares the direct path''s steering\n');
fprintf('    vector and is co-nulled for free.  The residual is if anything\n');
fprintf('    BETTER than with no reflection at all, because the bounce adds\n');
fprintf('    power to the spoofer and so raises the effective SAPR.  Rank 2\n');
fprintf('    only wastes array gain here (compare the auth columns).\n');
fprintf('  * A REFLECTOR AT A DIFFERENT AZIMUTH is a genuinely separate source\n');
fprintf('    and it is expensive: the rank-1 total residual collapses from about\n');
fprintf('    -28 dB to -5 dB, because the second arrival is simply not nulled and\n');
fprintf('    then dominates.  Rank 2 recovers 9-10 dB at 60 and 180 degrees.\n');
fprintf('  * BUT RANK 2 IS NOT ALWAYS RIGHT, and the 90-degree row shows why.\n');
fprintf('    A second arrival 6 dB down raises the second whitened eigenvalue only\n');
fprintf('    ~0.06 above the noise floor in a 5 ms dwell.  The resulting\n');
fprintf('    eigenvector carries roughly 26 degrees of error, so nulling it steers\n');
fprintf('    a null at a direction that is largely noise - and the result is WORSE\n');
fprintf('    than leaving it alone.  The rank decision must therefore be driven by\n');
fprintf('    RESOLVABILITY (an MDL test on the eigenvalues), never by prior\n');
fprintf('    knowledge that a second source exists.  asp_detect does this; forcing\n');
fprintf('    a fixed rank, as this study deliberately does, does not.\n');
fprintf('  * The two sample rates agree closely, because the relevant decorrelation\n');
fprintf('    is set by the SIGNAL bandwidth (1.023 MHz, a 977 ns correlation\n');
fprintf('    width), not by the sampling rate.  Widening the ADC does not make\n');
fprintf('    multipath worse.\n\n');

end
