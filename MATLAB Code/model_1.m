
%%

% time interval
interval = [0, 200];

fig = uifigure('Name', 'Solution Curves');
axes = uiaxes(fig, "XLim", interval, "YLim", [-10, 620]);
label = uilabel(fig);

%grayScaleFig = uifigure('Name', '2D Grayscale Plot');
%grayScaleFigAxes = uiaxes(grayScaleFig);

%movieFig = uifigure('Name', '2D Grayscale Plot Movie');

%Change this to adjust maximum and minimum m values of slider
minM = 0.01;
maxM = 0.5;

mSlider = uislider(fig, 'Value', minM, 'Limits', [minM, maxM], 'Position', [10, 400, 500, 3], 'ValueChangedFcn', @(mSlider, event) updateSolution(event, axes, interval, label));

%[xA, xB] = getUnstableManifoldXs(500, 600, 0.285, 0.3, 0.079, 0.071);
%disp(xA); disp(xB);

evaluateMosquitoODEs(0.0006, axes, interval, label);

%evaluateNHabitatODE(N, mArray, CCuArray, initialConditions, interval)
%evaluateNHabitatODE(3, [0, 0.0321, 1.654 * 10^-6; 0.0321, 0, 1.542 * 10^-5; 1.654 * 10^-6, 1.542 * 10^-5, 0], [310, 300, 290], [0, 310, 85, 300, 0, 290], interval);
%evaluateNHabitatODE(2, [0, 0.006; 0.006, 0], [200, 100], [50, 500, 0, 600], interval);
%mArray = [0, 0.0321, 0.0321, 0.0321, 0.0321; 0.0321, 0, 0.0087, 0.0018, 0.0087; 0.0321, 0.0087, 0, 0.0087, 0.0018; 0.0321, 0.0018, 0.0087, 0, 0.0087; 0.0321, 0.0087, 0.0018, 0.0087, 0];
%evaluateNHabitatODE(5, mArray, [500, 100, 120, 140, 160], [100, 500, 0, 100, 0, 120, 0, 140, 0, 160], interval);

%kernel(distance, MDT, q, time, sigma)
distance = 200 * sqrt(2);
%NegExpMigrationParam = generateNegExpKernelM(distance, 75, 0.8, 5);
%LogNormalMigrationParam = generateLogNormalKernelM(distance, 45.2, 0.95, 7, 74.4);

%plotxAxBPhasePlane(minxA, maxxA, CCuA, CCuB, iterator, bI, bU, m)
%plotxAxBPhasePlane(0, 70, 100, 1100, 0.1, 0.285, 0.3, 0.079, 0.071, 0.006);

%createKernelHeatmap("Log-Normal", 50, 75, 50);
%createKernelHeatmap("Negative Exponential", 120, 50);

%mArray = [0, 0.011, 0.011, 0.011, 0.011; 0.011, 0, 0.0024, 0.0003418, 0.0024; 0.011, 0.0024, 0, 0.0024, 0.0003418; 0.011, 0.0003418, 0.0024, 0, 0.0024; 0.011, 0.0024, 0.0003418, 0.0024, 0];
%solutionDevals = evaluateNHabitatODE(5, mArray, [50, 300, 300, 300, 300], [0, 50, 100, 300, 0, 300, 0, 300, 0, 300], interval);

%solutionDevals = evaluateNHabitatODE(2, [0, 0.006; 0.006, 0], [200, 100], [50, 500, 0, 600], interval);

%plotNHabitatsInSpace(time, N, locations, CCs, ODESolutionDevals, interval, fig, axes) 
%plotNHabitatsInSpace(0, 5, [0, 0; 0, 200; 200, 0; 0, -200; -200, 0], [50, 300, 300, 300, 300], solutionDevals, interval, grayScaleFig, grayScaleFigAxes);
%timeSlider = uislider(grayScaleFig, "Value", interval(1), 'Limits', interval, 'Position', [10, 400, 500, 3], 'ValueChangingFcn', @(timeSlider, event) updateGrayScalePlot(event, 5, [0, 0; 0, 200; 200, 0; 0, -200; -200, 0], [50, 300, 300, 300, 300], solutionDevals, interval, grayScaleFig, grayScaleFigAxes));
%createMovie(5, [0, 0; 0, 200; 200, 0; 0, -200; -200, 0], [50, 300, 300, 300, 300], solutionDevals, interval, grayScaleFig, grayScaleFigAxes, movieFig);

%solutionDevals = evaluateNHabitatODE(3, [0, 0.0110, 0.00001542; 0.0110, 0, 0.0003418; 0.00001542, 0.0003418, 0], [300, 300, 300], [100, 300, 0, 300, 0, 300], interval);
%plotNHabitatsInSpace(0, 3, [0, 0; 200, 0; 600, 0], [300, 300, 300], solutionDevals, interval, grayScaleFig, grayScaleFigAxes);
%timeSlider = uislider(grayScaleFig, "Value", interval(1), 'Limits', interval, 'Position', [10, 400, 500, 3], 'ValueChangingFcn', @(timeSlider, event) updateGrayScalePlot(event, 3, [0, 0; 200, 0; 600, 0], [300, 300, 300], solutionDevals, interval, grayScaleFig, grayScaleFigAxes));
%createMovie(3, [0, 0; 200, 0; 600, 0], [300, 300, 300], solutionDevals, interval, grayScaleFig, grayScaleFigAxes, movieFig);



function evaluateMosquitoODEs(m, axes, interval, label)
    CCuA = 0; CCuB = 0;
    % parameter values
    %bI = 0.285; bU = 0.3; deltaI = 0.079; deltaU = 0.071; dA = 0.001; dB = 0.001;
    bI = 0.285; bU = 0.3; deltaI = 0.079; deltaU = 0.071; CCuA = 400; CCuB = 600;
    
    %The separatrix ratio for the above parameter values is 12.0434782609

    %Calculate dA and dB based on CCuA and CCuB
    if CCuA ~= 0
        a = CCuA/CCuB;
        dA = (bU - deltaU - m * (1 - 1/a))/CCuA;
        dB = (bU - deltaU + m * (a - 1))/CCuB;
        criticalM = (bU - deltaU)/(1 - 1/a); 
    end
    
    % {xA, yA, xB, yB} is the system order
    ODESystem = @(t, Y)[
        (bI - deltaI) * Y(1) - dA * Y(1) * (Y(1) + Y(2)) - m * (Y(1) - Y(3));
        bU * Y(2) * Y(2) / (Y(1) + Y(2)) - deltaU * Y(2) - dA * Y(2) * (Y(1) + Y(2)) - m * (Y(2) - Y(4));
        (bI - deltaI) * Y(3) - dB * Y(3) * (Y(3) + Y(4)) + m * (Y(1) - Y(3));
        bU * Y(4) * Y(4) / (Y(3) + Y(4)) - deltaU * Y(4) - dB * Y(4) * (Y(3) + Y(4)) + m * (Y(2) - Y(4));
        ];
    
    % starting populations, in order xA, yA, xB, yB
    %initialConditions = [18.7102, 183.4650, 4.3971, 217.8785];
    %initialConditions = [4.3971, 217.8785, 18.7102, 183.4650];
    %initialConditions = [14.79333, 190.2066667, 15.79333, 190.2066667];
    initialConditions = [100, CCuA, 0, CCuB];
    
    solution = ode45(ODESystem, interval, initialConditions);
    finalPopulationXA = join(["xA: ", string(round(deval(solution, interval(2), 1)))]);
    finalPopulationYA = join(["yA: ", string(round(deval(solution, interval(2), 2)))]);
    finalPopulationXB = join(["xB: ", string(round(deval(solution, interval(2), 3)))]);
    finalPopulationYB = join(["yB: ", string(round(deval(solution, interval(2), 4)))]);
    finalPopulations = join([join(["m: ", string(m)]); ""; finalPopulationXA; finalPopulationYA; finalPopulationXB; finalPopulationYB], "", 2);
    hold on;
    fplot(axes, {@(t)deval(solution,t,1), @(t)deval(solution,t,2), @(t)deval(solution,t,3), @(t)deval(solution,t,4)}, interval);
    hold off;
    xlabel(axes, "Time (Days)", "FontName", "Times New Roman");
    ylabel(axes, "Mosquito Population", "FontName", "Times New Roman");
    title(axes, "Mosquito Populations over Time", "FontName", "Times New Roman");
    legend(axes, {'xA', 'yA', 'xB', 'yB'}, 'Location', 'northeast');
    label.Text = finalPopulations;
    label.Position = [420, 200, 100, 100];
    label.WordWrap = "on";

end

function updateSolution(event, axes, interval, label) 
    evaluateMosquitoODEs(event.Value, axes, interval, label);
end

function [xA, xB] = getUnstableManifoldXs(yA, yB, bI, bU, deltaI, deltaU)
    c = yA/yB;
    manifoldRatio = (bU - bI + deltaI - deltaU)/bU;
    xB = manifoldRatio * (yA + yB) / (c + 1 - (manifoldRatio * (c + 1)));
    xA = c * xB;
end

%%




function [migrationParam] = generateNegExpKernelM(x, mdt, q, time)
    lambda = 1/mdt;
    fun = @(theta, t) lambda .* t .* exp(-lambda .* sqrt((t.*cos(theta) + x).^2 + (t .* sin(theta)).^2))./(2 .* pi .* sqrt((t.*cos(theta) + x).^2 + (t .* sin(theta)).^2));
    radius = -mdt * log(1-q);
    disp("NegExponential Radius: " + radius);
    probabilityFun = integral2(fun, 0, 2 * pi, 0, radius);
    migrationParam = probabilityFun/time;
end

%sigma is the standard deviation
function [migrationParam] = generateLogNormalKernelM(x, mdt, q, time, sigma)
    M = log(mdt^2/(sqrt(mdt^2 + sigma^2)));
    S = sqrt(log((mdt^2 + sigma^2)/sigma^2));
    fun = @(theta, r) r .* exp((-1./(2 .* S.^2)) .* (log(sqrt((r .* cos(theta) - x).^2 + (r .* sin(theta)).^2))-M).^2) ./((r .* cos(theta) - x).^2 + (r .* sin(theta)).^2);
    radius = exp(M - sqrt(2) * S * erfinv(1 - 2 * q));
    disp("LogNormal Radius: " + radius);
    probabilityFun = (1/(S * (2 * pi)^(3/2))) * integral2(fun, 0, 2 * pi, 0, radius);
    migrationParam = probabilityFun/time;
end


%Given initial conditions for xA, yA, and yB, find the smallest xB that
%ensures infection fixation.
function [xB] = getMinimumxB(xA, yA, yB, iterator, bI, bU, deltaI, deltaU, dA, dB, m) 
    interval = [0, 300];
    
    ODESystem = @(t, Y)[
        (bI - deltaI) * Y(1) - dA * Y(1) * (Y(1) + Y(2)) - m * (Y(1) - Y(3));
        bU * Y(2) * Y(2) / (Y(1) + Y(2)) - deltaU * Y(2) - dA * Y(2) * (Y(1) + Y(2)) - m * (Y(2) - Y(4));
        (bI - deltaI) * Y(3) - dB * Y(3) * (Y(3) + Y(4)) + m * (Y(1) - Y(3));
        bU * Y(4) * Y(4) / (Y(3) + Y(4)) - deltaU * Y(4) - dB * Y(4) * (Y(3) + Y(4)) + m * (Y(2) - Y(4));
        ];
    
    finalxBPopulation = 0;
    initialxBPopulation = -iterator;
    while finalxBPopulation < 50
        initialxBPopulation = initialxBPopulation + iterator;
        % starting populations, in order xA, yA, xB, yB
        initialConditions = [xA, yA, initialxBPopulation, yB];
        solution = ode45(ODESystem, interval, initialConditions);
        finalxBPopulation = deval(solution, interval(2), 3);
    end
    xB = initialxBPopulation;
end


function plotxAxBPhasePlane(minxA, maxxA, CCuA, CCuB, iterator, bI, bU, deltaI, deltaU, m)
    yA = CCuA;
    yB = CCuB;
    a = CCuA/CCuB;
    dA = (bU - deltaU - m * (1 - 1/a))/CCuA;
    dB = (bU - deltaU + m * (a - 1))/CCuB;

    xAValues = [];
    xBValues = [];
    for xA = minxA:1:maxxA
        xB = getMinimumxB(xA, yA, yB, iterator, bI, bU, deltaI, deltaU, dA, dB, m);
        xAValues = [xAValues, xA];
        xBValues = [xBValues, xB];
    end

    %Create polygon for area above dataset:
    upperPolygonXVals = [xAValues, maxxA, minxA];
    upperPolygonYVals = [xBValues, max(xBValues) + 10, max(xBValues) + 10];
    %Create polygon for area below dataset:
    lowerPolygonXVals = [xAValues, maxxA, minxA];
    lowerPolygonYVals = [xBValues, 0, 0];
    %Fill polygons
    hold on;
    fill(upperPolygonXVals, upperPolygonYVals, [0 0.4470 0.7410], lowerPolygonXVals, lowerPolygonYVals, [0.4660 0.6740 0.1880], 'LineWidth', 1.5);
    ylim([0, max(xBValues) + 10]);

    xlabel("Habitat A Initial Infected Population (xA)", "FontName", "Times New Roman");
    ylabel("Habitat B Initial Infected Population (xB)", "FontName", "Times New Roman");
    title("Infected Mosquitoes Needed for Infection Fixation", "FontName", "Times New Roman");
    hold off;
end

% Creates heatmap for mosquito distributions from two habitats separated by
% a given distance.
function createKernelHeatmap(type, distance, MDT, sigma)
    %q is the proportion of mosquitoes bounded by the habitat radius
    q = 0.80;
    if type == "Negative Exponential"
        kernelfunc = @(x, y) 1./(2 .* pi .* MDT) .* (exp(-sqrt(x.^2 + y.^2)./MDT)./sqrt(x.^2 + y.^2) + exp(-sqrt((x-distance).^2 + y.^2)./MDT)./sqrt((x-distance).^2 + y.^2));
        radius = -MDT * log(1 - q);
        xVals = linspace(radius/2, distance - radius/2);
        yVals = linspace(-radius, radius);
    elseif type == "Log-Normal"
        M = log(MDT^2/(sqrt(MDT^2 + sigma^2)));
        S = sqrt(log((MDT^2 + sigma^2)/sigma^2));
        kernelfunc = @(x, y) 1./(S .* (2 .* pi).^(3/2)) .* (exp(-1./(2 .* S.^2) .* (log(sqrt(x.^2 + y.^2))-M).^2)./(x.^2 + y.^2) + exp(-1./(2 .* S.^2) .* (log(sqrt((x-distance).^2 + y.^2))-M).^2)./((x-distance).^2 + y.^2));
        radius = exp(M - sqrt(2) * S * erfinv(1 - 2 * q));
        xVals = linspace(-radius/5, distance + radius/5);
        yVals = linspace(-radius/5, radius/5);
    end
    [X, Y] = meshgrid(xVals, yVals);
    zVals = kernelfunc(X, Y);
    imagesc(xVals, yVals, zVals);
    colormap parula;
    xlabel("Horizontal Distance", "FontName", "Times New Roman");
    ylabel("Vertical Distance", "FontName", "Times New Roman");
    title("Log-Normal Mosquito Probability Density With Two Habitats", "FontName", "Times New Roman");
    colorbar;
end



% mArray is an N X N matrix, CCuArray is a 1D N X 1 array
function solutionDevals = evaluateNHabitatODE(N, mArray, CCuArray, initialConditions, interval)
    % constant parameter values
    bI = 0.285; bU = 0.3; deltaI = 0.079; deltaU = 0.071;
    
    %Calculate dArray based on CCuArray and mArray
    dArray = zeros(1, N);
    for i=1:1:N
        mTerm = 0;
        for j=1:1:N
            mTerm = mTerm + (mArray(i, j) * (CCuArray(j) - CCuArray(i)));
        end
        dArray(i) = ((bU - deltaU) * CCuArray(i) + mTerm)/(CCuArray(i) * CCuArray(i));
    end
    disp(dArray);

    % {x1, y1, x2, y2, ..., xN, yN} is the system order    
    function dydt = ODESystem(t, y, bI, bU, deltaI, deltaU, dArray, mArray, N)
        dydt = zeros(2 * N, 1);
        for k = 1:2:2*N
            mTermX = 0;
            mTermY = 0;
            for l=1:2:2*N
                mTermX = mTermX + mArray((k+1)/2, (l+1)/2) * (y(l) - y(k));
                mTermY = mTermY + mArray((k+1)/2, (l+1)/2) * (y(l+1) - y(k+1));
            end
            dydt(k) = (bI - deltaI) * y(k) - dArray((k+1)/2) * y(k) * (y(k) + y(k+1)) + mTermX;
            dydt(k+1) = bU * y(k+1) * y(k+1) / (y(k) + y(k+1)) - deltaU * y(k+1) - dArray((k+1)/2) * y(k+1) * (y(k) + y(k+1)) + mTermY;
        end
    end

    
    
    solution = ode45(@(t, y) ODESystem(t, y, bI, bU, deltaI, deltaU, dArray, mArray, N), interval, initialConditions);
    
    solutionDevals = cell(1, 2 * N);
    legendKey = cell(1, 2 * N);
    for i=1:1: 2 * N
        solutionDevals{i} = @(t)deval(solution, t, i);
        if mod(i, 2) == 1
            legendKey{i} = join(['x', string((i + 1)/2)], '');
        else
            legendKey{i} = join(['y', string((i)/2)], '');
        end
    end


    hold on;
    fplot(solutionDevals, interval);
    hold off;
    xlabel("Time (Days)", "FontName", "Times New Roman");
    ylabel("Mosquito Population", "FontName", "Times New Roman");
    title("Mosquito Populations over Time", "FontName", "Times New Roman");
    legend(legendKey, 'Location', 'northeast');

end

function plotNHabitatsInSpace(time, N, locations, CCs, ODESolutionDevals, interval, fig, axes) 
    axis(axes, "equal");
    axes.Color = [1, 1, 1];
    scale = 1/4;
    timePoint = time;
    for i = 1:1:N
        alpha = 1 - ODESolutionDevals{i * 2 - 1}(timePoint)/(ODESolutionDevals{i * 2 - 1}(timePoint) + ODESolutionDevals{i * 2}(timePoint));
        rectangle(axes, 'Position', [locations(i, 1) - CCs(i)*scale, locations(i, 2) - CCs(i)*scale, CCs(i)*scale*2, CCs(i)*scale*2], 'Curvature', [1,1], 'FaceColor', [alpha, alpha, alpha]);
    end

    xlabel(axes, "Position in X (m)", "FontName", "Times New Roman");
    ylabel(axes, "Position in Y (m)", "FontName", "Times New Roman");
    title(axes, "Regional Infection Densities", "FontName", "Times New Roman");
end


function updateGrayScalePlot(event, N, locations, CCs, ODESolutionDevals, interval, fig, axes)
    plotNHabitatsInSpace(event.Value, N, locations, CCs, ODESolutionDevals, interval, fig, axes);
end

function createMovie(N, locations, CCs, ODESolutionDevals, interval, grayScaleFig, grayScaleFigAxes, fig)
    v = VideoWriter('three_habitats_video_1.avi');
    open(v);
    for time = 1:1:interval(2)
        plotNHabitatsInSpace(time, N, locations, CCs, ODESolutionDevals, interval, grayScaleFig, grayScaleFigAxes)
        writeVideo(v, getframe(grayScaleFigAxes));
    end
    close(v);
end
