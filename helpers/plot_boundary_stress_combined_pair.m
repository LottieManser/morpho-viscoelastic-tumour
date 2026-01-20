function plot_boundary_stress_combined_pair(output1, output2, label1, label2)
% Plot stress and kinematics for two simulations:
%
% 5 rows x 3 columns:
% Row 1: total radial stress
% Row 2: elastic stress component (no pressure)
% Row 3: viscous stress component (no)
% Row 4: alpha (elastic stretch)
% Row 5: gamma (growth stretch)
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

%% ---------------- Extract all quantities ----------------
    function [sigT,sigE,sigV,alpha,gamma] = extractAll(output)

        nT = numel(output.ts);
        sigT  = nan(nT,3);
        sigE  = nan(nT,3);
        sigV  = nan(nT,3);
        alpha = nan(nT,3);
        gamma = nan(nT,3);

        locations = {'outer','necrotic','centre'};

        for k = 1:nT
            for j = 1:3
                idx = selectIndex(output,k,locations{j});
                if isnan(idx), continue; end

                r = output.rs(k,idx);
                R = output.RsMinusB(k,idx) + output.params.B;
                g = output.growthStretches(k,idx);
                a = r/(R*g);

                % time derivative of alpha
                if k==1
                    adot = 0;
                else
                    r0 = output.rs(k-1,idx);
                    R0 = output.RsMinusB(k-1,idx) + output.params.B;
                    g0 = output.growthStretches(k-1,idx);
                    a0 = r0/(R0*g0);
                    adot = (a-a0)/(output.ts(k)-output.ts(k-1));
                end

                sigT(k,j) = output.radialStresses(k,idx);
                sigE(k,j) = (a^2 - a^-4)/max(r,eps);
                sigV(k,j) = (1/(a*max(r,eps)))*adot;

                alpha(k,j) = a;
                gamma(k,j) = g;
            end
        end
    end

[sT1,sE1,sV1,a1,g1] = extractAll(output1);
[sT2,sE2,sV2,a2,g2] = extractAll(output2);

%% ---------------- Plotting ----------------
figure
tlo = tiledlayout(5,3,'TileSpacing','compact','Padding','compact');

rowLabels = {'Total stress',...
             'Elastic stress component',...
             'Viscous stress component',...
             'Alpha (elastic)',...
             'Gamma (growth)'};

colTitles = {'Outer boundary','Necrotic boundary','Centre'};

for row = 1:5
    for col = 1:3
        ax = nexttile((row-1)*3 + col);
        hold(ax,'on')

        switch row
            case 1
                y1 = sT1(:,col)/output1.params.L;
                y2 = sT2(:,col)/output2.params.L;
            case 2
                y1 = sE1(:,col)/output1.params.L;
                y2 = sE2(:,col)/output2.params.L;
            case 3
                y1 = sV1(:,col)/output1.params.L;
                y2 = sV2(:,col)/output2.params.L;
            case 4
                y1 = a1(:,col);
                y2 = a2(:,col);
            case 5
                y1 = g1(:,col);
                y2 = g2(:,col);
        end

        plot(ax,t1,y1,'k-','LineWidth',1.5)
        plot(ax,t2,y2,'Color',0.6*[1 1 1],'LineWidth',1.5)

        if col==1
            ylabel(ax,rowLabels{row},'FontWeight','bold')
        end
        if row==1
            title(ax,colTitles{col})
        end
        if row==5
            xlabel(ax,'$t/T$','Interpreter','latex')
        end

        legend(ax,label1,label2,'Location','best')
        box(ax,'on')
    end
end

end
