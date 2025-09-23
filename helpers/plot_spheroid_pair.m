function plot_spheroid_pair(output1, output2, label_1, label_2, frame)
% Plots two quarter spheroids at the same frame, arranged in 2x2 tiles
% per quantity, with shared colourbar per row.

if nargin < 5
    frame = size(output1.rs,1);
end

thetas = linspace(0,pi/2,1e2)';
x = @(r) r.*cos(thetas);
y = @(r) r.*sin(thetas);

fields = {'growthStretches','growthRates','radialStresses','nutrients','densities'};
shaderStrings = {...
    '$\gamma$',...
    '$\frac{1}{\gamma}\frac{\partial\gamma}{\partial t}$',...
    '$\sigma_r$',...
    '$c$',...
    '$\rho$'...
    };

numFields = length(fields);
figure
t = tiledlayout(numFields,2,'TileSpacing','compact','Padding','compact');

for i = 1:numFields
    field = fields{i};
    
    vals1 = output1.(field)(frame,:);
    vals2 = output2.(field)(frame,:);
    
    % Global colour limits
    clim = [min([vals1(:); vals2(:)]), max([vals1(:); vals2(:)])];
    if clim(1) == clim(2)
        clim = clim + [-1 1]*max(abs(clim(1)),1)*0.01;
    end
    
    % Left tile (Output 1)
    ax1 = nexttile;
    xs1 = x(output1.rs(frame,:));
    ys1 = y(output1.rs(frame,:));
    pcolor(ax1, xs1, ys1, repmat(vals1,length(thetas),1)); shading interp
    caxis(ax1, clim); axis(ax1,'equal','tight')
    title(ax1, shaderStrings{i},'Interpreter','latex','FontSize',16)
    ylabel(ax1,label_1,'FontSize',14)
    
    % Right tile (Output 2)
    ax2 = nexttile;
    xs2 = x(output2.rs(frame,:));
    ys2 = y(output2.rs(frame,:));
    pcolor(ax2, xs2, ys2, repmat(vals2,length(thetas),1)); shading interp
    caxis(ax2, clim); axis(ax2,'equal','tight')
    ylabel(ax2,label_2,'FontSize',14)
    
    % Colourbar attached to right tile
    cb = colorbar(ax2,'eastoutside');
    set(cb,'TickLabelInterpreter','latex')
end


title(t,sprintf('Spheroids at t = %.2f',output1.ts(frame)),...
      'Interpreter','latex','FontSize',20)
set(gcf,'Position',[200 100 1000 800])
end
