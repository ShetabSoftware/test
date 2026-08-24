function fmt = fx_fmt(wl, fl, varargin)
%FX_FMT Construct a signed two's-complement fixed-point format descriptor.
%
%   fmt = FX_FMT(WL, FL)
%   fmt = FX_FMT(WL, FL, 'round', RMODE, 'overflow', OMODE)
%
%   WL  total word length in bits, including the sign bit.
%   FL  fractional length in bits.  The integer word length is WL-FL-1.
%
%   Represented values are  q * 2^-FL  for integer q in [-2^(WL-1), 2^(WL-1)-1],
%   i.e. the real range is [-2^(WL-FL-1), 2^(WL-FL-1) - 2^-FL].
%
%   RMODE (default 'convergent'):
%     'convergent' round-half-to-even.  The ONLY rounding mode that is safe
%                  in front of a long coherent accumulator, because it has
%                  zero mean error.  See note below.
%     'nearest'    round-half-away-from-zero (MATLAB round()).  Mean error is
%                  zero for symmetric data but non-zero for one-sided data.
%     'floor'      truncation toward -Inf.  Mean error is -0.5 LSB.
%     'trunc'      truncation toward zero (magnitude truncation).
%
%   OMODE (default 'saturate'):
%     'saturate'   clamp to the representable range and report the event.
%     'wrap'       two's-complement wraparound.  Never use this on a signal
%                  path in a nulling receiver: a single wrap turns a strong
%                  interferer into broadband noise that no beamformer can
%                  remove.
%
%   WHY CONVERGENT ROUNDING IS THE DEFAULT
%   --------------------------------------
%   The spatial covariance estimator accumulates K ~ 1e4 products per
%   coherent interval.  If each product carries a rounding bias of e LSB,
%   every entry of R picks up the SAME bias K*e, so R acquires a rank-one
%   error term  K*e*ones(N)  whose eigenvector is the all-ones vector.  The
%   all-ones vector is the array response of a source at boresight.  A
%   truncating implementation therefore synthesises a phantom source at
%   zenith and steers a null into the satellites.  With WL=16, K=16368 and
%   'floor' the phantom sits at roughly 10*log10(K/12) = 31 dB above the
%   per-sample quantisation floor, which is enough to dominate a
%   0 dB-SAPR spoofer.  Convergent rounding removes the bias to first order;
%   carrying full-precision products into a wide accumulator removes it
%   exactly and is what the hardware actually does (see FX_COV_ACCUM).

if nargin < 2
    error('fx_fmt:args', 'fx_fmt requires WL and FL.');
end
if wl < 2 || wl ~= floor(wl)
    error('fx_fmt:wl', 'WL must be an integer >= 2.');
end
if fl ~= floor(fl)
    error('fx_fmt:fl', 'FL must be an integer.');
end

fmt.wl       = wl;
fmt.fl       = fl;
fmt.iwl      = wl - fl - 1;
fmt.round    = 'convergent';
fmt.overflow = 'saturate';

for k = 1:2:numel(varargin)
    key = lower(varargin{k});
    val = varargin{k+1};
    switch key
        case 'round'
            fmt.round = lower(val);
        case 'overflow'
            fmt.overflow = lower(val);
        otherwise
            error('fx_fmt:key', 'Unknown option "%s".', key);
    end
end

fmt.lsb    = 2^(-fl);
fmt.maxint = 2^(wl-1) - 1;
fmt.minint = -2^(wl-1);
fmt.maxval = fmt.maxint * fmt.lsb;
fmt.minval = fmt.minint * fmt.lsb;

end
