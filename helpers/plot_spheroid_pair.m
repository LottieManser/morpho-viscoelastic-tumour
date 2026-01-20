function plot_spheroid_pair(output1, output2, label_1, label_2, frame, coupleColorbars)
% Plots two quarter spheroids at the same frame, arranged in 2 rows x N columns
% per quantity, with shared colourbar per column. Labels only on left-most column.

%% -------------------- USER CONTROLS --------------------
useManualStressClim      = false;
manualStressClim         = [-100.5 -99.5];   % example limits for sigma

forceZeroCenteredStress  = false;             % toggle for stress behaviour
%% -------------------------------------------------------

% Blue-white-red colormap (diverging)
n = 256;
blue  = [0 0 1];
white = [1 1 1];
red   = [1 0 0];
cmapBWR = [ ...
    linspace(blue(1),  white(1), n/2)', ...
    linspace(blue(2),  white(2), n/2)', ...
    linspace(blue(3),  white(3), n/2)' ; ...
    linspace(white(1), red(1),   n/2)', ...
    linspace(white(2), red(2),   n/2)', ...
    linspace(white(3), red(3),   n/2)' ];

% Defaults
if nargin < 5
    frame = size(output1.rs,1);
end
if nargin < 6
    coupleColorbars = true;
end

thetas = linspace(0,pi/2,100)';

% Fields to plot
fields = {'growthStretches','growthRates','radialStresses','nutrients'};
shaderStrings = { ...
    '$\gamma$', ...
    '$\frac{1}{\gamma}\frac{\partial\gamma}{\partial t}$', ...
    '$\sigma_r$', ...
    '$c$' };

numFields = numel(fields);

figure
t = tiledlayout(2,numFields,'TileSpacing','compact','Padding','compact');

% Global geometry limits
all_rs = [output1.rs(frame,:) output2.rs(frame,:)];
xlims = [0 max(all_rs)];
ylims = [0 max(all_rs)];

%% ==================== MAIN LOOP ====================
for i = 1:numFields
    field = fields{i};

    vals1 = output1.(field)(frame,:);
    vals2 = output2.(field)(frame,:);

    % Min/max bookkeeping
    if coupleColorbars
        cmin = min([vals1(:); vals2(:)]);
        cmax = max([vals1(:); vals2(:)]);
    else
        cmin1 = min(vals1(:));  cmax1 = max(vals1(:));
        cmin2 = min(vals2(:));  cmax2 = max(vals2(:));
    end

    %% --------- COLOUR LIMITS LOGIC ---------
    if ismember(i,[2,3])   % diverging-capable fields

        if i == 3
            % ---------- STRESS ----------
            if useManualStressClim
                clim1 = manualStressClim;
                clim2 = manualStressClim;

            elseif forceZeroCenteredStress
                if coupleColorbars
                    maxabs = max(abs([vals1(:); vals2(:)]));
                    if maxabs == 0, maxabs = eps; end
                    clim1 = [-maxabs maxabs];
                    clim2 = clim1;
                else
                    clim1 = [-max(abs(vals1(:))) max(abs(vals1(:)))];
                    clim2 = [-max(abs(vals2(:))) max(abs(vals2(:)))];
                end

            else
                % raw min/max
                if coupleColorbars
                    clim1 = [cmin cmax];
                    clim2 = clim1;
                else
                    clim1 = [cmin1 cmax1];
                    clim2 = [cmin2 cmax2];
                end
            end

        else
            % ---------- GROWTH RATE (always symmetric) ----------
            if coupleColorbars
                maxabs = max(abs([cmin cmax]));
                clim1 = [-maxabs maxabs];
                clim2 = clim1;
            else
                clim1 = [-max(abs(vals1(:))) max(abs(vals1(:)))];
                clim2 = [-max(abs(vals2(:))) max(abs(vals2(:)))];
            end
        end

    else
        % ---------- NON-DIVERGING FIELDS ----------
        if coupleColorbars
            clim1 = [cmin cmax];
            clim2 = clim1;
        else
            clim1 = [cmin1 cmax1];
            clim2 = [cmin2 cmax2];
        end
    end

    %% --------- DATA MATRICES ---------
    vals1_mat = repmat(vals1, numel(thetas), 1);
    vals2_mat = repmat(vals2, numel(thetas), 1);

    [R1, T1] = meshgrid(output1.rs(frame,:), thetas);
    [R2, T2] = meshgrid(output2.rs(frame,:), thetas);

    X1 = R1 .* cos(T1);   Y1 = R1 .* sin(T1);
    X2 = R2 .* cos(T2);   Y2 = R2 .* sin(T2);

    %% --------- TOP ROW ---------
    ax1 = nexttile(i);
    pcolor(ax1,X1,Y1,vals1_mat); shading(ax1,'interp');
    axis(ax1,'equal'); xlim(ax1,xlims); ylim(ax1,ylims);
    caxis(ax1,clim1);

    if i == 2 || (i == 3 && forceZeroCenteredStress)
        colormap(ax1,cmapBWR)
    else
        colormap(ax1,parula)
    end

    title(ax1,shaderStrings{i},'Interpreter','latex','FontSize',16);
    if i==1, ylabel(ax1,label_1,'FontSize',14); end

    hold(ax1,'on');
    B1 = output1.rs(frame,end);
    plot(ax1,B1*cos(thetas),B1*sin(thetas),'k-','LineWidth',1.2);

    cb1 = colorbar(ax1,'eastoutside');
    set(cb1,'TickLabelInterpreter','latex');

    %% --------- BOTTOM ROW ---------
    ax2 = nexttile(i+numFields);
    pcolor(ax2,X2,Y2,vals2_mat); shading(ax2,'interp');
    axis(ax2,'equal'); xlim(ax2,xlims); ylim(ax2,ylims);
    caxis(ax2,clim2);

    if i == 2 || (i == 3 && forceZeroCenteredStress)
        colormap(ax2,cmapBWR)
    else
        colormap(ax2,parula)
    end

    if i==1, ylabel(ax2,label_2,'FontSize',14); end

    hold(ax2,'on');
    B2 = output2.rs(frame,end);
    plot(ax2,B2*cos(thetas),B2*sin(thetas),'k-','LineWidth',1.2);

    cb2 = colorbar(ax2,'eastoutside');
    set(cb2,'TickLabelInterpreter','latex');
end

title(t, sprintf('Spheroids at t = %.2f',output1.ts(frame)), ...
    'Interpreter','latex','FontSize',20);

end
