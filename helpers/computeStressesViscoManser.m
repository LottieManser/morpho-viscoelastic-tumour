function [radialStress, hoopStress, bulkStress, viscoelasticStretch] = computeStressesViscoManser(RMinusB,r,drdt,growthStretch,dgrowthStretchdt,params)
%% Compute the radial and hoop stresses from the radial coordinates and growth
%% stretches at the material points in RMinusB.

    % Compute the radial stress first. ELASTIC
    elasticIntegrand = 2*growthStretch.*(r.^6 - growthStretch.^6.*(RMinusB+params.B).^6)./r.^7;
    % Ignore the value at R=0, which is NaN.
    elasticIntegrand(1) = 0;

    % Compute the radial stress first. VISCO
    viscoIntegrand = 12*growthStretch.^2.*((growthStretch.*drdt-r.*dgrowthStretchdt)./(r.^4)).*(RMinusB+params.B).^2;
    % Ignore the value at R=0, which is NaN.
    viscoIntegrand(1) = 0;

    % The integrand is poorly computed around R = 0, so we use a two-term
    % Taylor expansion of the integrand around this point. We use the Taylor
    % expansion until we reach a threshold of R = params.radialStressIntegrandThreshold.
    growthStretchPrime = gradient(growthStretch, RMinusB);
    growthStretchPrimePrime = gradient(growthStretchPrime, RMinusB);
    dotgrowthStretchPrime = gradient(dgrowthStretchdt, RMinusB);
    dotgrowthStretchPrimePrime = gradient(dotgrowthStretchPrime, RMinusB);

    approxElasticIntegrand = 2*(-3/2 * growthStretchPrime(1)/growthStretch(1) + ... 
                     (3/80 * (growthStretchPrime(1)/growthStretch(1))^2 - ...
                     6/5 * growthStretchPrimePrime(1)/growthStretch(1))*(RMinusB+params.B));
    approxViscoIntegrand = -3*( ...
    ((growthStretch(1).^2).*dotgrowthStretchPrime(1) - growthStretch(1).*growthStretchPrime(1).*dgrowthStretchdt(1))./(growthStretch(1).^3) + ...
        (RMinusB+params.B)*(13*(growthStretchPrime(1).^2).*dgrowthStretchdt(1) - 8*growthStretch(1).*growthStretchPrimePrime(1).*dgrowthStretchdt(1) - 13*growthStretch(1).*growthStretchPrime(1).*dotgrowthStretchPrime(1) + 8*(growthStretch(1).^2).*dotgrowthStretchPrimePrime(1))./(10*growthStretch(1).^3)...
        );

    elasticIntegrandComposite = approxElasticIntegrand.*(RMinusB<=(params.radialStressIntegrandThreshold-params.B)) + ...
                        elasticIntegrand.*(RMinusB>(params.radialStressIntegrandThreshold-params.B));
    viscoIntegrandComposite = approxViscoIntegrand.*(RMinusB<=(params.radialStressIntegrandThreshold-params.B)) + ...
                        viscoIntegrand.*(RMinusB>(params.radialStressIntegrandThreshold-params.B));
    elasticIntegral = cumtrapz(RMinusB,elasticIntegrandComposite);
    viscoIntegral = cumtrapz(RMinusB,viscoIntegrandComposite);

    % Anchor the viscoelastic stress the same way the elastic-only version does:
    % subtract the end-values of the integrals, and use only kappa for the boundary correction.
    radialStress = params.mu*(elasticIntegral - elasticIntegral(end)) + params.varsigma*(viscoIntegral - viscoIntegral(end)) - params.kappa*(r(end) - params.B)/params.B;

    
    % Compute the hoop stress from the radial stress.
    viscoelasticStretch = r  ./ ((RMinusB+params.B).*growthStretch);
    % Ignore the value at R=0, which is NaN.
    viscoelasticStretch(1) = 1;

    % As before, we'll use a two-term Taylor expansion to compute the
    % elastic stretch near R=0.
    approxViscoElasticStretch = 1 - (growthStretchPrime(1)/growthStretch(1)).*(RMinusB+params.B)/4;
    viscoelasticStretch = approxViscoElasticStretch.*(RMinusB<=(params.elasticStretchIntegrandThreshold-params.B)) + ...
                        viscoelasticStretch.*(RMinusB>(params.elasticStretchIntegrandThreshold-params.B));

    hoopStress = radialStress + params.mu*(viscoelasticStretch.^2 - 1./viscoelasticStretch.^4);

    % Compute the bulk stress from the radial and hoop stress.
    bulkStress = radialStress + 2*hoopStress;

end