function elastic_and_growth_stretches(output1, output2, label1, label2)
% Plot kinematics for two simulations:
%
% 4 rows x 3 columns:
% Row 1: gamma (growth stretch)
% Row 2: d gamma / dt
% Row 3: r (current radius)
% Row 4: d r / dt
%
% Columns: outer | necrotic | centre
% Colours: black = sim1, grey = sim2

if nargin < 3, label1 = 'Sim 1'; end
if nargin < 4, label2 = 'Sim 2'; end

t1 = output1.ts / output1.params.T;
t2 = output2.ts / output2.params.T;

%% ---------------- Helper: index selection ----------------
    function idx = selectIndex(output,k,which)
        switch which
            case 'outer'
                idx = size(output.rs,2);
            case 'necrotic'
                rn = output.necroticRadii(k);
                if ~isfinite(rn)
                    idx = NaN;
                    return
                end
                [~,idx] = min(abs(output.rs(k,:) - rn));
            case 'centre'
                idx = 1;
        end
    end

%% ---------------- Extract kinematics ----------------
    function [gamma,gammadot,r,rdot] = extractAll(output)

        nT = numel(output.ts);

        gamma    = nan(nT,3);
        gammadot = nan(nT,3);
        r        = nan(nT,3);
        rdot     = nan(nT,3);

        locations = {'outer','necrotic','centre'};

        for k = 1:nT
            for j = 1:3
                idx = selectIndex(output,k,locations{j});
                if isnan(idx), continue; end

                rk = output.rs(k,idx);
                g  = output.growthStretches(k,idx);

                if k == 1
                    gdot = 0;
                    rdk  = 0;
                else
                    dt = output.ts(k) - output.ts(k-1);

                    g0   = output.growthStretches(k-1,idx);
                    gdot = (g - g0)/dt;

                    r0   = output.rs(k-1,idx);
                    rdk  = (rk - r0)/dt;
                end

                gamma(k,j)    = g;
                gammadot(k,j) = gdot;
                r(k,j)        = rk;
                rdot(k,j)     = rdk;
            end
        end
    end

[g1,gd1,r1,rd1] = extractAll(output1);
[g2,gd2,r2,rd2] = extractAll(output2);

%% ---------------- Plotting ----------------
figure
tlo = tiledlayout(4,3,'TileSpacing','compact','Padding','compact');

rowLabels = {'$\gamma$ (growth)',...
             '$\dot{\gamma}$',...
             '$r$',...
             '$\dot{r}$'};

colTitles = {'Outer boundary','Necrotic boundary','Centre'};

for row = 1:4
    for col = 1:3
        ax = nexttile((row-1)*3 + col);
        hold(ax,'on')

        switch row
            case 1
                y1 = g1(:,col);  y2 = g2(:,col);
            case 2
                y1 = gd1(:,col); y2 = gd2(:,col);
            case 3
                y1 = r1(:,col);  y2 = r2(:,col);
            case 4
                y1 = rd1(:,col); y2 = rd2(:,col);
        end

        plot(ax,t1,y1,'k-','LineWidth',1.5)
        plot(ax,t2,y2,'Color',0.6*[1 1 1],'LineWidth',1.5)

        if col == 1
            ylabel(ax,rowLabels{row},'Interpreter','latex','FontWeight','bold')
        end
        if row == 1
            title(ax,colTitles{col})
        end
        if row == 4
            xlabel(ax,'$t/T$','Interpreter','latex')
        end

        legend(ax,label1,label2,'Location','best')
        box(ax,'on')
    end
end

end
