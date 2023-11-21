function uce = unconderr(sigmat, mut, logoutput)
% Computes the unconditional error (or logarithm of unconditional error if
% logoutput == 1). sigmat and mut need to be arrays with the same size or
% either one needs to be a scalar. 

if nargin < 3
    logoutput = 0;
end

nm = numel(mut); ns = numel(sigmat);
if ~(ns == 1 || nm == 1 || all(size(sigmat)==size(mut)))
    error('Incorrect inputs for unconderr.');
end
uce = NaN((nm>ns)*size(mut)+(nm<=ns)*size(sigmat));
negind = (mut < 0);
mutn = mut(negind);
mutp = mut(~negind);
if ns == 1 || nm == 1
    sigmatn = sigmat;
    sigmatp = sigmat;
else
    sigmatn = sigmat(negind);
    sigmatp = sigmat(~negind);
end

rn = []; rp = [];
if logoutput
    if ~isempty(mutn)
        ln0 = log_normcdf_fast(mutn./sigmatn);
        ln1 = mutn + 0.5*sigmatn.^2 + log_normcdf_fast(-(mutn+sigmatn.^2)./sigmatn);
        ln2 = log(2) + mutn + 0.5*sigmatn.^2 + log_normcdf_fast(-sigmatn);
        mn = max([ln0,ln1,ln2],[],2);
        rn = mn + log(exp(ln0-mn) + exp(ln1-mn) - exp(ln2-mn));
    end
    if ~isempty(mutp)
        lp1 = log_normcdf_fast(-mutp./sigmatp);
        lp2 = mutp + 0.5*sigmatp.^2 + log_normcdf_fast(-(mutp+sigmatp.^2)./sigmatp);
        mp = min(lp1,lp2);
        rp = mp + log(exp(lp1-mp) - exp(lp2-mp));
    end
else
    if ~isempty(mutn)
        % the equations below are slightly reformulated to avoid some numerical errors
        rn = normcdf_fast(mutn./sigmatn) ...
            + exp(mutn + 0.5*sigmatn.^2 + log_normcdf_fast(-(mutn+sigmatn.^2)./sigmatn))...
            - 2*exp(mutn + 0.5*sigmatn.^2 + log_normcdf_fast(-sigmatn));
    end
    if ~isempty(mutp)
        rp = normcdf_fast(-mutp./sigmatp) - ...
            exp(mutp + 0.5*sigmatp.^2 + log_normcdf_fast(-(mutp+sigmatp.^2)./sigmatp));
    end
end
if nm == 1
    uce = [rn rp];
else
    uce(negind) = rn;
    uce(~negind) = rp;
end
end

