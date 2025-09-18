function plot_evolution_pair(output1, output2, label_1, label_2)
% Plot the evolution of spheroid radii for two outputs on the same axes.
    
    figure
    hold on
    
    % Time normalized
    t1 = output1.ts / output1.params.T;
    t2 = output2.ts / output2.params.T;
    
    % Outer radius
    plot(t1, output1.rs(:,end) / output1.params.L, 'k-', 'LineWidth', 1.5)
    plot(t2, output2.rs(:,end) / output2.params.L, 'k--', 'LineWidth', 1.5)
    
    % Necrotic radius
    plot(t1, output1.necroticRadii / output1.params.L, 'Color',0.7*[1,1,1], 'LineWidth',1)
    plot(t2, output2.necroticRadii / output2.params.L, 'Color',0.7*[1,1,1], 'LineStyle','--','LineWidth',1)
    
    % Optional nutrient-free radius (if exists)
    if any(output1.nutrients(:) == 0) || any(output2.nutrients(:) == 0)
        % Output 1
        nutrientsTemp1 = output1.nutrients; nutrientsTemp1(output1.nutrients > 0) = 1;
        [~,zeroNutrientInds1] = max(nutrientsTemp1,[],2,'linear');
        plot(t1, output1.rs(zeroNutrientInds1)/output1.params.L, 'Color',0.4*[1,1,1], 'LineWidth',1)
        
        % Output 2
        nutrientsTemp2 = output2.nutrients; nutrientsTemp2(output2.nutrients > 0) = 1;
        [~,zeroNutrientInds2] = max(nutrientsTemp2,[],2,'linear');
        plot(t2, output2.rs(zeroNutrientInds2)/output2.params.L, 'Color',0.4*[1,1,1], 'LineStyle','--','LineWidth',1)
        
        legendEntries = { ...
            strcat('Outer radius ', label_1), strcat('Outer radius ', label_2), ...
            strcat('Necrotic radius ', label_1), strcat('Necrotic radius ', label_2), ...
            strcat('Nutrient-free radius ', label_1), strcat('Nutrient-free radius ', label_2) ...
        };
    else
        legendEntries = { ...
            strcat('Outer radius ', label_1), strcat('Outer radius ', label_2), ...
            strcat('Necrotic radius ', label_1), strcat('Necrotic radius ', label_2) ...
        };
    end
    
    legend(legendEntries,'Location','southeast')
    box on
    xlabel('$t/T$','Interpreter','latex')
    ylabel('$r/L$','Interpreter','latex')
    title('Spheroid radius evolution')
end
