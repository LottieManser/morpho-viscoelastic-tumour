function [muVals, varsigmaVals] = dynamicMaterialProperties(density, params)

    muVals = params.mu0 * (1+density);% ./ (1 - density ) ;
    varsigmaVals = params.varsigma0;
end