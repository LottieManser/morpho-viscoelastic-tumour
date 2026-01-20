function [muVals, varsigmaVals] = staticMaterialProperties(params)
    % Constant material parameters
    muVals = params.mu0;
    varsigmaVals = params.varsigma0;
end
