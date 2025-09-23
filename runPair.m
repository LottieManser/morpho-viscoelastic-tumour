addpath(genpath('./helpers'))

%% Specify model parameters and setup. All quantities are in SI units.
params = struct();

% Model identifier, one of 'A', 'B', 'C', 'D', or 'E', corresponding to the
% accompanying publication.
params.model = 'E';

% Time and discretisation.
params.tFinal = 20; % Final time.
params.nR = 2000; % Number of points in the spatial discretisation.
params.nT = 3000; % Number of points in the temporal discretisation.

% Material parameters.
params.kappa = 0; % Spring constant in the radial stress boundary condition.

% Growth rate parameters.
params.k = 1; % The basic growth rate constant.
params.sigmaHat = -1; % Threshold below which growth is arrested due to compressive stress, if included in the model.
params.beta = 1; % Scale factor for the local argument of n, if included in the model.

% Nutrient parameters.
params.cInf = 1; % Concentration of nutrient at the boundary.
params.cHat = 0.8; % Nutrient threshold for necrosis.
params.lambda = 0.9; % Rate of nutrient consumption by tissue.
params.D = 1; % Diffusion coefficient of nutrient inside the tumour.
params.L = sqrt(params.D * params.cInf/params.lambda); % The diffusive lengthscale.
params.bHat = sqrt(6)*params.L; % The threshold radius of the tumour for full/partial perfusion.

% Initial configuration.
params.B = params.L; % Initial (and unstressed) radius of the spheroid.

% Timescale
params.T = 1 / (params.k * params.cInf);

% Constant thresholds for Taylor expanding integrals.
params.radialStressIntegrandThreshold = 0.05; % Radial threshold for using Taylor expansion of radial stress integrand.
params.elasticStretchIntegrandThreshold = 0.05; % Radial threshold for using Taylor expansion of of elastic stretch integrand.


params1 = params;
params2 = params;

%% Pick run specific parameters

params1.mu = 1; % Shear modulus of the tumour material.
params1.varsigma = 0;
params2.mu = 0.5; % Shear modulus of the tumour material.
params2.varsigma = 0.5;


%% Run the simulation.
outputElastic = runSimViscoManser(params1);
outputVisco = runSimViscoManser(params2);

%% Plotting.
%plot_spheroid_pair(outputElastic, outputVisco, 'Elastic', 'Viscoelastic', 1);
plot_spheroid_pair(outputElastic, outputVisco, 'Elastic', 'Viscoelastic',params.nT);
plot_evolution_pair(outputElastic, outputVisco,'Elastic', 'Viscoelastic');
plot_boundary_stress_pair(outputElastic, outputVisco,'Elastic', 'Viscoelastic');

%plot_spheroid(output, params.nT);
%plot_evolution(output);
%make_spheroid_video(output,10,10,'spheroid_growth_E_high_res.mp4');  

%% Saving.
output = outputVisco;

% Generate reduced output.
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

% save('output.mat','output')
save('reducedOutputVisco.mat','reducedOutput')
