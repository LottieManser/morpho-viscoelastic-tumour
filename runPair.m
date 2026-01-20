addpath(genpath('./helpers'))

%% Specify model parameters and setup. All quantities are in SI units.
params = struct();

% Model identifier, one of 'A', 'B', 'C', 'D', or 'E', corresponding to the
% accompanying publication.
params.model = 'B';

% Time and discretisation.
params.tFinal = 15; % 20; % Final time.
params.nR = 3000; % Number of points in the spatial discretisation.
params.nT = 4000; % Number of points in the temporal discretisation.

% Material parameters.
params.kappa = 100; %316.2; % Spring constant in the radial stress boundary condition.

% Growth rate parameters.
params.k = 1; % The basic growth rate constant.
params.sigmaHat = -100; % Threshold below which growth is arrested due to compressive stress, if included in the model.
params.beta = 1; % Scale factor for the local argument of n, if included in the model.

% Nutrient parameters.
params.cInf = 5; % Concentration of nutrient at the boundary.
params.cHat = 4; % Nutrient threshold for necrosis.
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

params1.label = 'elastic';
params1.mu0 = 0.1;
params1.varsigma0 = 0;
%params1.jammingeffect = 0;
params1.materialmodel = 'static';
params1.model = 'B';


params2.label = 'viscoelastic';
params2.mu0 = 0.1;
params2.varsigma0 = 5;
%params2.jammingeffect = 0;
params2.materialmodel = 'static';
params2.model = 'B';


%% Run the simulation.
output1 = runSimViscoManser(params1);
output2 = runSimViscoManser(params2);

%% Plotting.
%plot_spheroid_pair(output1, output2, 'Elastic', 'Viscoelastic', 1);

plot_spheroid_pair(output1, output2, params1.label, params2.label, params.nT, false);
plot_boundary_stress_pair(output1, output2, params1.label, params2.label)
%plot_evolution_pair(output1, output2,params1.label, params2.label);
elastic_and_growth_stretches(output1, output2,params1.label, params2.label);
%plot_boundary_stress_pair(output1, output2, params1.label, params2.label);
%plot_boundary_stress_combined_pair(output1, output2, params1.label, params2.label)
%plot_material_pair(output1, output2, params1.label, params2.label);


%plot_spheroid(output, params.nT);
%plot_evolution(output);
%make_spheroid_video(output,10,10,'spheroid_growth_E_high_res.mp4');  

%% Saving
output = output2;

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
