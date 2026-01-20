function plot_boundary_stress_pair(output1, output2, label_1, label_2)
% Plot evolution of radial stress:
% Rows:
%   Row 1: simulation 1
%   Row 2: simulation 2
% Columns:
%   (1) outer boundary
%   (2) necrotic boundary
%   (3) centre

figure
tlo = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

% Time normalized
t1 = output1.ts / output1.params.T;
t2 = output2.ts / output2.params.T;

%% ---------- Helper: necrotic stress ----------
    function sigN = necroticStress(output,t)
        sigN = nan(size(t));
        for k = 1:numel(t)
            rn = output.necroticRadii(k);
            if isfinite(rn)
                [~,idx] = min(abs(output.rs(k,:) - rn));
                sigN(k) = output.radialStresses(k,idx);
            end
        end
    end

sigN1 = necroticStress(output1,t1);
sigN2 = necroticStress(output2,t2);

%% ---------- Row 1: simulation 1 ----------
nexttile(1)
plot(t1, output1.radialStresses(:,end) / output1.params.L, ...
    'k-', 'LineWidth',1.5)
title('Outer boundary')
ylabel({label_1,'$\sigma_r$'},'Interpreter','latex','FontWeight','bold')
box on

nexttile(2)
plot(t1, sigN1 / output1.params.L, 'k-', 'LineWidth',1.5)
title('Necrotic boundary')
box on

nexttile(3)
plot(t1, output1.radialStresses(:,1) / output1.params.L, ...
    'k-', 'LineWidth',1.5)
title('Centre')
box on

%% ---------- Row 2: simulation 2 ----------
nexttile(4)
plot(t2, output2.radialStresses(:,end) / output2.params.L, ...
    'Color',0.5*[1 1 1], 'LineWidth',1.5)
ylabel({label_2,'$\sigma_r$'},'Interpreter','latex','FontWeight','bold')
xlabel('$t/T$','Interpreter','latex')
box on

nexttile(5)
plot(t2, sigN2 / output2.params.L, ...
    'Color',0.5*[1 1 1], 'LineWidth',1.5)
xlabel('$t/T$','Interpreter','latex')
box on

nexttile(6)
plot(t2, output2.radialStresses(:,1) / output2.params.L, ...
    'Color',0.5*[1 1 1], 'LineWidth',1.5)
xlabel('$t/T$','Interpreter','latex')
box on

end
