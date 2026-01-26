addpath(genpath('./helpers'))

%% Specify model parameters
params = struct();
params.model = 'B';
params.tFinal = 30;
params.nR = 3000;
params.nT = 5000;
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



%% Copy params for elastic and viscoelastic runs
params1 = params;
params2 = params;

%% Elastic case
params1.label = 'elastic';
params1.mu0 = 1;
params1.varsigma0 = 5;
params1.materialmodel = 'static';
params1.tau = 1000; % relaxation timescale

%% Viscoelastic case
params2.label = 'viscoelastic';
params2.mu0 = 1;
params2.varsigma0 = 5;
params2.materialmodel = 'static';
params2.tau = 0.01; % relaxation timescale

%% Run simulations
output1 = runSimViscoManser(params1);
output2 = runSimViscoManser(params2);

%% Plot results
plot_spheroid_pair(output1, output2, params1.label, params2.label, params.nT, false);
plot_boundary_stress_pair(output1, output2, params1.label, params2.label);
%elastic_and_growth_stretches(output1, output2, params1.label, params2.label);
plot_evolution_pair(output1, output2, params1.label, params2.label);

%% Save reduced output
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
reducedOutput.params = params;
reducedOutput.necroticRadii = output.necroticRadii;

save('reducedOutputVisco.mat','reducedOutput');
