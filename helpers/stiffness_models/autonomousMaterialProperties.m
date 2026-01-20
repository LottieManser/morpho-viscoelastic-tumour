function [muVals, varsigmaVals] = autonomousMaterialProperties(t, params)
    % Linearly increase mu and varsigma over time
    % t = scalar current time

    finalFactor = 10; % multiply mu0 by 2 at end of simulation
    muVals = params.mu0 * (1 + (finalFactor-1) * t / params.tFinal);
    varsigmaVals = params.varsigma0 * (1 + (finalFactor-1) * t / params.tFinal);
end
