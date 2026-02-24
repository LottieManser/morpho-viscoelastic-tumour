function growthRate = modelB2GrowthRate(nutrient, r, necroticRadius, radialStress, params)
%% Returns the growthRate (1/gamma * d(gamma)/dt) for Model B2:
% "Model B with stress-regulated necrotic remodelling/compaction"
%
% Option 1 (recommended): keep the intuitive biological split:
% - viable rim: proliferation is nutrient-driven and inhibited by stress
% - necrotic core: volumetric remodelling/resorption is also mechanically regulated
%
% In this model:
%   if r >= necroticRadius (viable):
%       growthRate = k  * (c - cHat) * n(sigma_r)
%   if r <  necroticRadius (necrotic):
%       growthRate = kN * (c - cHat) * nN(sigma_r)
%
% with defaults:
%   kN  = k   (if not provided)
%   nN  = nFun (same stress gate) (if not provided)

    % --- defaults for necrotic kinetics ---
    if ~isfield(params,'kN')
        params.kN = params.k;                 % necrotic remodelling rate scale
    end
    if ~isfield(params,'nNFun')
        % If you don't provide a separate function handle, reuse nFun
        % You can override by setting params.nNFun = @(sigma,params) ...
        params.nNFun = @(sigma,p) nFun(sigma,p);
    end

    % --- viable growth: as before ---
    growthRateGrowing = params.k * (nutrient - params.cHat) .* nFun(radialStress, params);

    % --- necrotic remodelling: stress-regulated too ---
    growthRateNecrotic = params.kN * (nutrient - params.cHat) .* params.nNFun(radialStress, params);

    % Identify viable tissue (outside necrotic radius)
    growingMask = r >= necroticRadius;

    % Combine
    growthRate = growthRateGrowing .* growingMask + growthRateNecrotic .* (1 - growingMask);
end