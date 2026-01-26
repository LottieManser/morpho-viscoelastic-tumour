function output = runSimViscoManser(params, progressFlag)

    if nargin < 2
        progressFlag = true;
    end

    %% Initial discretisation
    RsMinusB = zeros(params.nT, params.nR);
    RsMinusB(1,:) = linspace(0, params.B, params.nR) - params.B;
    ts = linspace(0, params.tFinal, params.nT);

    %% Preallocate
    rs = zeros(params.nT, params.nR);
    drdts = zeros(params.nT, params.nR);
    drdts(1,:) = 0;
    growthStretches = ones(params.nT, params.nR);
    dgrowthStretchesdt = zeros(params.nT, params.nR);
    growthRates = zeros(params.nT, params.nR);
    radialStresses = zeros(params.nT, params.nR);
    hoopStresses = zeros(params.nT, params.nR);
    bulkStresses = zeros(params.nT, params.nR);
    viscoelasticStretches = zeros(params.nT, params.nR);
    nutrients = zeros(params.nT, params.nR);
    necroticRadii = -Inf * ones(params.nT,1);
    densities = zeros(params.nT, params.nR);
    mus = zeros(params.nT, params.nR); mus(1,:) = params.mu0;
    varsigmas = zeros(params.nT, params.nR); varsigmas(1,:) = params.varsigma0;

    %% Progress bar
    if progressFlag
        clear('textprogressbar.m')
        textprogressbar('Simulating growth: ');
    end

    %% Time loop
    for tInd = 1:params.nT

        if progressFlag
            textprogressbar(100 * tInd / params.nT);
        end

        %% Eulerian radial coordinates
        rs(tInd,:) = computeRadialCoord(RsMinusB(tInd,:), growthStretches(tInd,:), params);

        %% Update material parameters
        switch params.materialmodel
            case 'static'
                [mu_now, varsigma_now] = staticMaterialProperties(params);
                mus(tInd,:) = mu_now; varsigmas(tInd,:) = varsigma_now;
            case 'autonomous'
                [mu_now, varsigma_now] = autonomousMaterialProperties(ts(tInd), params);
                mus(tInd,:) = mu_now; varsigmas(tInd,:) = varsigma_now;
            case 'dynamic'
                [mus(tInd,:), varsigmas(tInd,:)] = dynamicMaterialProperties(densities(tInd), params);
        end

        %% Compute stresses (with relaxation)
        if tInd == 1
            radialStressPrev = [];
            dt = 0;
        else
            radialStressPrev = radialStresses(tInd-1,:);
            dt = ts(tInd) - ts(tInd-1);
        end

        [radialStresses(tInd,:), ...
         hoopStresses(tInd,:), ...
         bulkStresses(tInd,:), ...
         viscoelasticStretches(tInd,:)] = ...
            computeStressesViscoManser( ...
                RsMinusB(tInd,:), rs(tInd,:), drdts(tInd,:), ...
                growthStretches(tInd,:), dgrowthStretchesdt(tInd,:), ...
                mus(tInd,:), varsigmas(tInd,:), ...
                params, dt, radialStressPrev);
        
        %% Nutrient
        nutrients(tInd,:) = computeNutrient(rs(tInd,:), params);

        %% Necrotic core
        necroticRadii(tInd) = computeNecroticRadius(necroticRadii(1:tInd-1), ...
                                                    rs(tInd,:), nutrients(tInd,:), params);

        %% Growth rate
        switch params.model
            case 'A'
                growthRates(tInd,:) = modelAGrowthRate(nutrients(tInd,:), params);
            case 'B'
                growthRates(tInd,:) = modelBGrowthRate(nutrients(tInd,:), rs(tInd,:), ...
                                                      necroticRadii(tInd), radialStresses(tInd,:), params);
            case 'B2'
                growthRates(tInd,:) = modelB2GrowthRate(nutrients(tInd,:), rs(tInd,:), ...
                                                       necroticRadii(tInd), radialStresses(tInd,:), params);
            case 'C'
                growthRates(tInd,:) = modelCGrowthRate(nutrients(tInd,:), rs(tInd,:), ...
                                                      necroticRadii(tInd), radialStresses(tInd,:), params);
            case 'D'
                growthRates(tInd,:) = modelDGrowthRate(nutrients(tInd,:), rs(tInd,:), ...
                                                      necroticRadii(tInd), radialStresses(tInd,:), ...
                                                      hoopStresses(tInd,:), params);
            case 'E'
                growthRates(tInd,:) = modelEGrowthRate(nutrients(tInd,:), rs(tInd,:), ...
                                                      necroticRadii(tInd), radialStresses(tInd,:), ...
                                                      hoopStresses(tInd,:), params);
        end

        %% Density
        densities(tInd,:) = 1 ./ (growthStretches(tInd,:).^3);

        %% Advance to next time step
        if tInd < params.nT
            RsMinusB(tInd+1,:) = interp1(rs(tInd,:), RsMinusB(tInd,:), linspace(rs(tInd,1), rs(tInd,end), params.nR));
            rRemapped = interp1(RsMinusB(tInd,:), rs(tInd,:), RsMinusB(tInd+1,:));
            growthStretchRemapped = interp1(RsMinusB(tInd,:), growthStretches(tInd,:), RsMinusB(tInd+1,:));
            growthRateRemapped = interp1(RsMinusB(tInd,:), growthRates(tInd,:), RsMinusB(tInd+1,:));
            drdtRemapped = interp1(RsMinusB(tInd,:), drdts(tInd,:), RsMinusB(tInd+1,:));
            dgrowthStretchesdtRemapped = interp1(RsMinusB(tInd,:), dgrowthStretchesdt(tInd,:), RsMinusB(tInd+1,:));

            drdts(tInd+1,:) = drdtRemapped;
            dgrowthStretchesdt(tInd+1,:) = dgrowthStretchesdtRemapped;

            dtLocal = ts(tInd+1) - ts(tInd);
            growthStretches(tInd+1,:) = growthStretchRemapped + dtLocal * growthStretchRemapped .* growthRateRemapped;
            dgrowthStretchesdt(tInd+1,:) = (growthStretches(tInd+1,:) - growthStretches(tInd,:)) / dtLocal;
        end
    end

    if progressFlag
        textprogressbar('\nDone.');
    end

    %% Pack outputs
    output = struct();
    output.ts = ts;
    output.RsMinusB = RsMinusB;
    output.rs = rs;
    output.drdts = drdts;
    output.growthStretches = growthStretches;
    output.dgrowthStretchesdt = dgrowthStretchesdt;
    output.growthRates = growthRates;
    output.radialStresses = radialStresses;
    output.hoopStresses = hoopStresses;
    output.bulkStresses = bulkStresses;
    output.viscoelasticStretches = viscoelasticStretches;
    output.nutrients = nutrients;
    output.necroticRadii = necroticRadii;
    output.params = params;
    output.densities = densities;
    output.mus = mus;
    output.varsigmas = varsigmas;

end
