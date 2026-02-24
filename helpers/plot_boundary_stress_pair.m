function plot_boundary_stress_pair(output1, output2, label_1, label_2)
% Plot evolution of radial stress.
% Rows:
%   Row 1: simulation 1 (typically elastic)
%   Row 2: simulation 2 (typically SLS)
% Columns:
%   (1) outer boundary
%   (2) necrotic boundary
%   (3) centre
%
% For simulation 2, additionally plots equilibrium and Maxwell (overstress)
% radial stress contributions if fields exist:
%   output2.radialStressesEq, output2.radialStressesOv

figure
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

% Time normalised
t1 = output1.ts / output1.params.T;
t2 = output2.ts / output2.params.T;

% Check availability of decomposed stresses in output2
hasDecomp = isfield(output2,'radialStressesEq') && isfield(output2,'radialStressesOv');

%% ---------- Helper: necrotic index and stress ----------
    function [sigN, idxN] = necroticStressAndIndex(output, tNorm)
        sigN = nan(size(tNorm));
        idxN = nan(size(tNorm));
        for k = 1:numel(tNorm)
            rn = output.necroticRadii(k);
            if isfinite(rn)
                [~,idx] = min(abs(output.rs(k,:) - rn));
                idxN(k) = idx;
                sigN(k) = output.radialStresses(k,idx);
            end
        end
    end

[sigN1, idxN1] = necroticStressAndIndex(output1, t1);
[sigN2, idxN2] = necroticStressAndIndex(output2, t2);

% If decomposed stresses exist, compute necrotic-boundary series for them too
if hasDecomp
    sigN2_eq = nan(size(t2));
    sigN2_ov = nan(size(t2));
    for k = 1:numel(t2)
        idx = idxN2(k);
        if isfinite(idx)
            sigN2_eq(k) = output2.radialStressesEq(k,idx);
            sigN2_ov(k) = output2.radialStressesOv(k,idx);
        end
    end
end

%% ---------- Row 1: simulation 1 ----------
nexttile(1)
plot(t1, output1.radialStresses(:,end) / output1.params.L, 'k-', 'LineWidth',1.5)
title('Outer boundary')
ylabel({label_1,'$\sigma_r$'},'Interpreter','latex','FontWeight','bold')
box on

nexttile(2)
plot(t1, sigN1 / output1.params.L, 'k-', 'LineWidth',1.5)
title('Necrotic boundary')
box on

nexttile(3)
plot(t1, output1.radialStresses(:,1) / output1.params.L, 'k-', 'LineWidth',1.5)
title('Centre')
box on

%% ---------- Row 2: simulation 2 ----------
% Styling choices: total = darker grey, eq/ov = dashed variants
colTot = 0.5*[1 1 1];
colEq  = 0.2*[1 1 1];
colOv  = 0.7*[1 1 1];

nexttile(4)
hold on
hTot = plot(t2, output2.radialStresses(:,end) / output2.params.L, '-', 'Color',colTot, 'LineWidth',1.5);
if hasDecomp
    hEq  = plot(t2, output2.radialStressesEq(:,end) / output2.params.L, '--', 'Color',colEq, 'LineWidth',1.2);
    hOv  = plot(t2, output2.radialStressesOv(:,end) / output2.params.L, ':',  'Color',colOv, 'LineWidth',1.8);
end
ylabel({label_2,'$\sigma_r$'},'Interpreter','latex','FontWeight','bold')
xlabel('$t/T$','Interpreter','latex')
box on

% Put legend only once
if hasDecomp
    legend([hTot,hEq,hOv], {'total','equilibrium','Maxwell'}, ...
        'Location','best', 'Box','off');
end
hold off

nexttile(5)
hold on
plot(t2, sigN2 / output2.params.L, '-', 'Color',colTot, 'LineWidth',1.5)
if hasDecomp
    plot(t2, sigN2_eq / output2.params.L, '--', 'Color',colEq, 'LineWidth',1.2)
    plot(t2, sigN2_ov / output2.params.L, ':',  'Color',colOv, 'LineWidth',1.8)
end
xlabel('$t/T$','Interpreter','latex')
box on
hold off

nexttile(6)
hold on
plot(t2, output2.radialStresses(:,1) / output2.params.L, '-', 'Color',colTot, 'LineWidth',1.5)
if hasDecomp
    plot(t2, output2.radialStressesEq(:,1) / output2.params.L, '--', 'Color',colEq, 'LineWidth',1.2)
    plot(t2, output2.radialStressesOv(:,1) / output2.params.L, ':',  'Color',colOv, 'LineWidth',1.8)
end
xlabel('$t/T$','Interpreter','latex')
box on
hold off

end