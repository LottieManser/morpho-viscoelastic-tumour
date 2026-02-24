addpath(genpath('./helpers'))

%% Specify model parameters
params = struct();
params.model = 'B';
params.tFinal = 6;
params.nR = 2000;
params.nT = 3000;
params.kappa = 100;
params.k = 1;
params.sigmaHat = -100;
params.beta = 1;
params.cInf = 5;
params.cHat = 4;
params.lambda = 0.9;
params.D = 1;
params.L = sqrt(params.D * params.cInf / params.lambda);
params.B = params.L;
params.bHat = sqrt(6) * params.L; % threshold radius for nutrient perfusion
params.T = 1 / (params.k * params.cInf); % timescale
params.radialStressIntegrandThreshold = 0.05;
params.elasticStretchIntegrandThreshold = 0.05;

%% Copy params for elastic-like and relaxing SLS runs
params1 = params;
params2 = params;

%% Elastic case (pure neo-Hookean)
params1.label = 'elastic';
params1.materialmodel = 'static';
params1.muInf0 = 1;   % this is the elastic shear modulus
params1.muV0   = 0;   % NO Maxwell branch
params1.gA0    = 1;   % irrelevant when muV0 = 0

%% Viscoelastic SLS case
params2.label = 'SLS';
params2.materialmodel = 'static';
params2.muInf0 = 0;
params2.muV0   = 1;
params2.gA0    = 50;

%% Run simulations
output1 = runSimViscoManser(params1);
output2 = runSimViscoManser(params2);

%% Plot results
plot_spheroid_pair(output1, output2, params1.label, params2.label, params.nT, false);
plot_boundary_stress_pair(output1, output2, params1.label, params2.label);
plot_evolution_pair(output1, output2, params1.label, params2.label);

%% Save reduced output (optionally include alphaE)
output = output2;
reducedOutput = struct();
reducedOutput.ts = output.ts;
reducedOutput.radii = output.rs(:,end);
reducedOutput.rFinal = output.rs(end,:);
reducedOutput.growthRatesFinal = output.growthRates(end,:);
reducedOutput.growthStretchesFinal = output.growthStretches(end,:);
reducedOutput.nutrientsFinal = output.nutrients(end,:);
reducedOutput.radialStressCentre = output.radialStresses(:,1);
reducedOutput.hoopStressCentre = output.hoopStresses(:,1);
reducedOutput.radialStressFinal = output.radialStresses(end,:);
reducedOutput.hoopStressFinal = output.hoopStresses(end,:);
reducedOutput.necroticRadii = output.necroticRadii;

% SLS internal variable (useful diagnostics)
reducedOutput.alphaECentre = output.alphaEs(:,1);
reducedOutput.alphaEFinal  = output.alphaEs(end,:);

reducedOutput.params = params;

save('reducedOutputSLS.mat','reducedOutput');