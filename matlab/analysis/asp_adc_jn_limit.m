function jnLimitDB = asp_adc_jn_limit(wl, backoffDB, allowedLossDB)
%ASP_ADC_JN_LIMIT Interference-to-noise headroom of a quantised array receiver.
%
%   jnLimitDB = ASP_ADC_JN_LIMIT(WL, BACKOFFDB, ALLOWEDLOSSDB)
%
%   WL             ADC / datapath word length in bits (including sign)
%   BACKOFFDB      AGC backoff: composite RMS below full scale, in dB
%   ALLOWEDLOSSDB  tolerated post-nulling SNR loss, in dB (default 1.0)
%
%   Returns the largest interference-to-thermal-noise ratio, in dB, that the
%   converter supports before quantisation noise costs more than
%   ALLOWEDLOSSDB of post-nulling C/N0.
%
%   DERIVATION
%   ----------
%   An AGC in an anti-jam receiver holds the COMPOSITE signal (thermal noise
%   plus interference) at a fixed backoff from full scale, because that is
%   what prevents clipping.  As the interferer grows, the AGC therefore
%   pushes the thermal noise DOWN.  In LSBs,
%
%       sigma_thermal = A_FS * 10^(-backoff/20) * 10^(-JN/20)            (1)
%
%   with A_FS = 2^(WL-1).
%
%   The beamformer removes the interferer, because the interferer is
%   spatially rank-one and coherent across channels.  It does NOT remove
%   quantisation noise, because each converter quantises independently, so
%   quantisation noise is spatially white and lands in the passband of every
%   beam.  After nulling, the noise floor is therefore
%
%       sigma_total^2 = sigma_thermal^2 + q^2/12,   q = 1 LSB             (2)
%
%   and the post-nulling loss is
%
%       L = 10*log10(1 + (1/12)/sigma_thermal^2)   [dB]                   (3)
%
%   Setting L = ALLOWEDLOSSDB and solving (1) and (3) for JN gives the
%   answer.  This is the single most important converter specification for a
%   CRPA and it is the reason 12-bit SDR transceivers are adequate for
%   spoofing (SAPR 0-30 dB) but marginal for jamming (J/N 50-90 dB).
%
%   Worked values, backoff = 14 dB, 1 dB allowed loss:
%       12 bits  ->  ~57 dB J/N     (AD9361 raw; ~51 dB at 10.5 bit ENOB)
%       14 bits  ->  ~69 dB J/N
%       16 bits  ->  ~81 dB J/N

if nargin < 3 || isempty(allowedLossDB)
    allowedLossDB = 1.0;
end

% From (3): sigma_thermal^2 = (1/12) / (10^(L/10) - 1)
lossLin = 10^(allowedLossDB/10);
if lossLin <= 1
    jnLimitDB = -Inf;
    return;
end
sigmaThermal = sqrt((1/12) / (lossLin - 1));

% From (1): JN = 20*log10(A_FS) - backoff - 20*log10(sigma_thermal)
aFs = 2^(wl-1);
jnLimitDB = 20*log10(aFs) - backoffDB - 20*log10(sigmaThermal);

end
