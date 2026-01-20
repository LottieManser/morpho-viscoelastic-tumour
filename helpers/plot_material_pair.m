function plot_material_pair(output1, output2, label_1, label_2)
% Plots mu and varsigma for two simulations:
% Rows = simulations, Columns = [μ field, varsigma field, μ evolution, varsigma evolution]

fields = {'mus','varsigmas'};
titles = {'$\mu$ (shear modulus)', '$\varsigma$ (viscosity)'};

% --- Global colour limits for consistency ---
clims = cell(1,numel(fields));
for fIdx = 1:numel(fields)
    allVals = [output1.(fields{fIdx})(:); output2.(fields{fIdx})(:)];
    cmin = min(allVals);
    cmax = max(allVals);
    if cmin == cmax
        cmin = cmin - 0.01*abs(cmin+1);
        cmax = cmax + 0.01*abs(cmax+1);
    end
    clims{fIdx} = [cmin, cmax];
end

% --- Layout ---
figure('Position',[200 200 1600 650])
t = tiledlayout(2,4,'TileSpacing','compact','Padding','compact');

for simIdx = 1:2
    if simIdx == 1
        out = output1; label = label_1;
    else
        out = output2; label = label_2;
    end
    
    % Geometry setup
    thetas = linspace(0,pi/2,100)';
    x = @(r) r.*cos(thetas);
    y = @(r) r.*sin(thetas);
    frame = size(out.rs,1);
    
    % μ and varsigma fields
    for fIdx = 1:2
        field = fields{fIdx};
        vals = out.(field)(frame,:);
        xs = x(out.rs(frame,:));
        ys = y(out.rs(frame,:));
        
        ax = nexttile((simIdx-1)*4 + fIdx);
        pcolor(ax, xs, ys, repmat(vals,length(thetas),1));
        shading interp; axis equal tight
        caxis(ax, clims{fIdx});
        cb = colorbar(ax);
        cb.TickLabelInterpreter = 'latex';
        title(ax, char(titles{fIdx}), 'Interpreter','latex', ...
            'FontSize',14, 'FontWeight','bold', ...
            'Units','normalized','Position',[0.5,1.05,0], ...
            'HorizontalAlignment','center');
        if fIdx == 1
            ylabel(ax, label, 'FontSize', 12)
        end
    end
    
    % --- μ evolution ---
    ax = nexttile((simIdx-1)*4 + 3);
    hold(ax,'on')
    tvals = out.ts / out.params.T;
    coreIdx = 1;
    surfIdx = size(out.rs,2);
    plot(ax, tvals, out.mus(:,coreIdx), 'k--', 'LineWidth', 1.8)
    plot(ax, tvals, out.mus(:,surfIdx), 'k-',  'LineWidth', 1.8)
    xlabel(ax, '$t/T$', 'Interpreter','latex')
    ylabel(ax, '$\mu$', 'Interpreter','latex')
    legend(ax, {'$\mu_{core}$','$\mu_{boundary}$'}, 'Interpreter','latex','Location','best')
    box(ax,'on')
    
    % --- varsigma evolution ---
    ax = nexttile((simIdx-1)*4 + 4);
    hold(ax,'on')
    plot(ax, tvals, out.varsigmas(:,coreIdx), 'Color', 0.5*[1 1 1], 'LineStyle','--', 'LineWidth', 1.8)
    plot(ax, tvals, out.varsigmas(:,surfIdx),  'Color', 0.5*[1 1 1], 'LineStyle','-',  'LineWidth', 1.8)
    xlabel(ax, '$t/T$', 'Interpreter','latex')
    ylabel(ax, '$\varsigma$', 'Interpreter','latex')
    legend(ax, {'$\varsigma_{core}$','$\varsigma_{boundary}$'}, 'Interpreter','latex','Location','best')
    box(ax,'on')
end
end
