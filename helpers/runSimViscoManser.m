function output = runSimViscoManser(params, progressFlag)
% runSimViscoManser
% Finite-deformation SLS (Huber–Tsakmakis Model A) implementation.
%
% REQUIREMENTS:
%   computeStressesSLS must have signature:
%   [radialStress, hoopStress, bulkStress, elasticStretch, alphaE, ...
%    radialStressEq, radialStressOv, hoopStressEq, hoopStressOv] = ...
%      computeStressesSLS(RMinusB, r, drdt, gamma, dGammadt, muInf, muV, gA, params, dt, alphaEPrev)
%
% and you need: computeRadialCoord, computeNutrient, computeNecroticRadius,
% and the growth-rate functions.

    if nargin < 2
        progressFlag = true;
    end

    %% Initial discretisation
    RsMinusB = zeros(params.nT, params.nR);
    RsMinusB(1,:) = linspace(0, params.B, params.nR) - params.B;
    ts = linspace(0, params.tFinal, params.nT);

    %% Preallocate (Eulerian radius, kinematics, growth)
    rs = zeros(params.nT, params.nR);
    drdts = zeros(params.nT, params.nR);
    drdts(1,:) = 0;

    growthStretches = ones(params.nT, params.nR);           % gamma(R,t)
    dgrowthStretchesdt = zeros(params.nT, params.nR);       % dot{gamma}
    growthRates = zeros(params.nT, params.nR);

    %% Preallocate stresses and derived fields
    radialStresses = zeros(params.nT, params.nR);
    hoopStresses   = zeros(params.nT, params.nR);
    bulkStresses   = zeros(params.nT, params.nR);

    % NEW: split stresses (equilibrium vs Maxwell/overstress)
    radialStressesEq = zeros(params.nT, params.nR);
    radialStressesOv = zeros(params.nT, params.nR);
    hoopStressesEq   = zeros(params.nT, params.nR);
    hoopStressesOv   = zeros(params.nT, params.nR);

    % In SLS, elastic stretch alpha = r/(R*gamma)
    elasticStretches = zeros(params.nT, params.nR);

    % Internal Maxwell-branch elastic stretch alpha_e(R,t)
    alphaEs = ones(params.nT, params.nR);

    %% Nutrient / necrosis / density
    nutrients = zeros(params.nT, params.nR);
    necroticRadii = -Inf * ones(params.nT,1);
    densities = zeros(params.nT, params.nR);

    %% Material parameters for SLS
    muInfs = zeros(params.nT, params.nR); muInfs(1,:) = params.muInf0;
    muVs   = zeros(params.nT, params.nR); muVs(1,:)   = params.muV0;
    gAs    = zeros(params.nT, params.nR); gAs(1,:)    = params.gA0;

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

        %% Update material parameters (SLS)
        switch params.materialmodel
            case 'static'
                muInfs(tInd,:) = params.muInf0;
                muVs(tInd,:)   = params.muV0;
                gAs(tInd,:)    = params.gA0;

            case 'autonomous'
                error('materialmodel="autonomous" not implemented for SLS. Provide autonomousMaterialPropertiesSLS.');

            case 'dynamic'
                error('materialmodel="dynamic" not implemented for SLS. Provide dynamicMaterialPropertiesSLS.');

            otherwise
                error('Unknown params.materialmodel: %s', params.materialmodel);
        end

        %% Time step + previous internal state
        if tInd == 1
            dt = 0;
            alphaEPrev = [];
        else
            dt = ts(tInd) - ts(tInd-1);
            alphaEPrev = alphaEs(tInd-1,:);
        end

        %% Compute stresses (SLS)
        [radialStresses(tInd,:), ...
         hoopStresses(tInd,:), ...
         bulkStresses(tInd,:), ...
         elasticStretches(tInd,:), ...
         alphaEs(tInd,:), ...
         radialStressesEq(tInd,:), ...
         radialStressesOv(tInd,:), ...
         hoopStressesEq(tInd,:), ...
         hoopStressesOv(tInd,:)] = ...
            computeStressesSLS( ...
                RsMinusB(tInd,:), rs(tInd,:), drdts(tInd,:), ...
                growthStretches(tInd,:), dgrowthStretchesdt(tInd,:), ...
                muInfs(tInd,:), muVs(tInd,:), gAs(tInd,:), ...
                params, dt, alphaEPrev);

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

            otherwise
                error('Unknown params.model: %s', params.model);
        end

        %% Density
        densities(tInd,:) = 1 ./ (growthStretches(tInd,:).^3);

        %% Advance to next time step (remap fields onto new R-grid)
        if tInd < params.nT

            % Remap R grid to keep uniform sampling in Eulerian space
            RsMinusB(tInd+1,:) = interp1(rs(tInd,:), RsMinusB(tInd,:), ...
                                         linspace(rs(tInd,1), rs(tInd,end), params.nR));

            % Remap primary fields
            rRemapped                 = interp1(RsMinusB(tInd,:), rs(tInd,:), RsMinusB(tInd+1,:));
            growthStretchRemapped     = interp1(RsMinusB(tInd,:), growthStretches(tInd,:), RsMinusB(tInd+1,:));
            growthRateRemapped        = interp1(RsMinusB(tInd,:), growthRates(tInd,:), RsMinusB(tInd+1,:));
            drdtRemapped              = interp1(RsMinusB(tInd,:), drdts(tInd,:), RsMinusB(tInd+1,:));
            dgrowthStretchesdtRemapped= interp1(RsMinusB(tInd,:), dgrowthStretchesdt(tInd,:), RsMinusB(tInd+1,:));

            % Remap internal state alphaE
            alphaERemapped            = interp1(RsMinusB(tInd,:), alphaEs(tInd,:), RsMinusB(tInd+1,:));

            drdts(tInd+1,:)              = drdtRemapped;
            dgrowthStretchesdt(tInd+1,:) = dgrowthStretchesdtRemapped;

            alphaEs(tInd+1,:)            = alphaERemapped;

            % Advance growth stretch (explicit Euler)
            dtLocal = ts(tInd+1) - ts(tInd);
            growthStretches(tInd+1,:) = growthStretchRemapped + ...
                                        dtLocal * growthStretchRemapped .* growthRateRemapped;

            dgrowthStretchesdt(tInd+1,:) = (growthStretches(tInd+1,:) - growthStretches(tInd,:)) / dtLocal;

            %#ok<NASGU>
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

    % Split components (for debugging/plots)
    output.radialStressesEq = radialStressesEq;
    output.radialStressesOv = radialStressesOv;
    output.hoopStressesEq   = hoopStressesEq;
    output.hoopStressesOv   = hoopStressesOv;

    output.elasticStretches = elasticStretches;
    output.alphaEs = alphaEs;

    output.nutrients = nutrients;
    output.necroticRadii = necroticRadii;

    output.params = params;
    output.densities = densities;

    output.muInfs = muInfs;
    output.muVs = muVs;
    output.gAs = gAs;
end