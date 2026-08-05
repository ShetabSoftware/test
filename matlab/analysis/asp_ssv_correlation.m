function [rho, sinThetaDB] = asp_ssv_correlation(yHat, bTrue)
%ASP_SSV_CORRELATION Normalised inner product between estimated and true SSV.
%
%   [RHO, SINTHETADB] = ASP_SSV_CORRELATION(YHAT, BTRUE)
%
%   RHO         |yHat' * bTrue| / (||yHat|| * ||bTrue||), in [0,1].  This is
%               the quantity plotted in Figure 3 of the paper.
%   SINTHETADB  20*log10(sqrt(1 - rho^2)), the null depth this estimate can
%               achieve, in dB.
%
%   THE CONVERSION MATTERS AND THE PAPER DOES NOT MAKE IT
%   ----------------------------------------------------
%   A rank-one projector built from yHat leaves a residual gain toward the
%   true direction of exactly sin^2(theta) = 1 - rho^2 (relative to the mean
%   sky gain).  So Figure 3 of the paper, which tops out near rho = 0.99,
%   is implicitly a statement that the achievable null depth saturates
%   around
%
%       20*log10(sqrt(1 - 0.99^2)) = -17 dB
%
%   not the 25 dB quoted from the single simulation instance in Figure 5.
%   Reading rho = 0.9 at 0 dB SAPR off the same figure gives -7.2 dB.  Any
%   product specification has to be written against this curve, not against
%   one favourable geometry.

if isempty(yHat) || isempty(bTrue)
    rho = NaN; sinThetaDB = NaN;
    return;
end

yHat = yHat(:);
bTrue = bTrue(:);

ny = norm(yHat);
nb = norm(bTrue);
if ny < eps || nb < eps
    rho = 0;
else
    rho = abs(yHat' * bTrue) / (ny*nb);
end
rho = min(rho, 1);

sinThetaDB = 20*log10(max(sqrt(max(1 - rho^2, 0)), realmin));

end
