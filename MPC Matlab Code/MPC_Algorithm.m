
% Constants
g = 9.81; % m/Iner
rho = Inf;  % m
Width = 3; % m
Mass = 2050; % kg
Inertia = 3344; % kg-m^2
FStiff = 1050;   % Front wheel cornering stiffness
RStiff_In = 95;   % Rear wheel inner cornering stiffness
RStiff_Out = 8000;   % Rear wheel outer cornering stiffness
% mu = .3;  % Slippery
Vr = 14.0; % m/s
% Variables
Cx       = .01;
S0  = [0 0 0 0 0]';
inputMin = [-pi/3 1]';
inputMax = [pi/3 1]';
Period  = 0.1;
Ts  = 25;
Np   = 45;
stateMin = [];
stateMax = [];
% Geometry
a = 2;
b = 2;
c = 1.8;
% Constraints
SM = [0 0 0 1 0; 0 0 0 -1 0; 0 0 0 0 1;0 0 0 0 -1;1/Vr a/Vr 0 0 -1;-1/Vr -a/Vr 0 0 1;1/Vr -b/Vr 0 0 0;-1/Vr b/Vr 0 0 0;];
ICM = [];
eqICM = [0 1];
SCV = [Width/2;Width/2;pi/3;pi/3;pi/6;pi/6;pi/6;pi/6;];
ICV = [];
eqICV = 1;

obstacle1.duration  = [3 4];
obstacle1.width     = [1 -0.5];
obstacle2.duration  = [12 13];
obstacle2.width     = [.2 -0.8];
obstacle3.duration  = [19 20];
obstacle3.width     = [0.8 -0.1];
obstacles = [obstacle1 obstacle2 obstacle3];
[SMTimeVarying, SCVTimeVarying, timeVaryingSteps] = CreateObstacle(obstacles, Period);

showPlots = true;
Psi_r = 1/rho;

% State Space Modeling
states = ["y_dot", "Phi_dot", "e_Phi", "e_y", "delta"];
A_c = [(FStiff + RStiff_In)/(Vr*Mass), (-Vr + (FStiff*a - b*RStiff_In)/(Vr*Mass)), 0, 0, -FStiff/Mass; (a*FStiff - b*RStiff_Out)/(Inertia*Vr), (a^2*FStiff + b^2*RStiff_Out)/(Inertia*Vr), 0, 0, -a*FStiff/Inertia; 0, 1, 0, 0, 0; 1, 0, Vr, 0, 0; 0, 0, 0, 0, 0;];
B_c = [0 0 0 0 1; 0 0 -Vr*Psi_r 0 0]';
C_c = [0 0 0 1 0];
D_c = [0 0];

% Discretize
[A, B, C, D] = DSS(A_c, B_c, C_c, D_c, Period);
% MPC Control
[X, U, Y] = ModelPredictiveControl(A, B, C, D, Np, Period, Ts, Cx, S0, inputMin, inputMax, stateMin, stateMax, SM, ICM, eqICM, SMTimeVarying, SCVTimeVarying, timeVaryingSteps, SCV, ICV, eqICV, showPlots);

% Initialization
filename = 'MPC_LaneKeeping.gif';
delayTime = 0.1; % Delay time between frames in seconds

% Plot Road
figure;
hold on
timeVec = 0:Period:Ts;
plot(timeVec, (Width/2)*ones(size(timeVec)), "r--", "DisplayName", "Boundary")
plot(timeVec, -(Width/2)*ones(size(timeVec)), "r--", "HandleVisibility", "off")
plot(timeVec, zeros(size(timeVec)), "b--", "DisplayName", "Center Line")
plotObstacles(obstacles);
legend
title("Lane Keeping Validation and Obstacle Avoidance at Speed 14.0 m/s")
xlabel("Time [sec]")
ylabel("e_y(t) [m]")

% Initialize the plot for Y data
h = plot(timeVec(1), Y(1), 'DisplayName', 'Position');
scalingFactor = 8; % Adjust this to make the arrow larger
q = quiver(timeVec(1), Y(1), 0, 0, 'r', 'MaxHeadSize', 20, 'LineWidth', 2, 'DisplayName', 'Direction');

% Create GIF
for k = 1:length(timeVec)
    % Update plot data
    set(h, 'XData', timeVec(1:k), 'YData', Y(1:k))
    
    % Update quiver (arrow) position and direction
    if k > 1
        dx = (timeVec(k) - timeVec(k-1)) * scalingFactor;
        dy = (Y(k) - Y(k-1)) * scalingFactor;
        set(q, 'XData', timeVec(k-1), 'YData', Y(k-1), 'UData', dx, 'VData', dy)
    end
    
    % Capture the frame
    drawnow
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    % Write to the GIF File
    if k == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', delayTime);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', delayTime);
    end
end


function [constraint_matrix, bounds_vector, time_vector] = CreateObstacle(obstacles, sampling_time)
    constraint_matrix = [];
    bounds_vector = [];
    time_vector = cell(size(obstacles));
    
    % Loop through each obstacle
    for i = 1:length(obstacles)
        obstacle = obstacles(i);
        
        % Check if the mean width of the obstacle is negative
        if mean(obstacle.width) < 0
            % Append e_y constraint for negative width
            constraint_matrix = [constraint_matrix; 0 0 0 -1 0];
            % Append maximum negative width to the bounds vector
            bounds_vector = [bounds_vector; -max(obstacle.width)];
        else
            % Append e_y constraint for positive width
            constraint_matrix = [constraint_matrix; 0 0 0 1 0];
            % Append minimum positive width to the bounds vector
            bounds_vector = [bounds_vector; min(obstacle.width)];
        end
        
        % Create a time vector for the current obstacle
        time_vector{i} = (obstacle.duration(1)/sampling_time) : (obstacle.duration(2)/sampling_time);
    end
end
% MPC Dual Mode %
function [stateTrajectories, inputTrajectories, outputTrajectories] = ModelPredictiveControl(A, B, C, D, horizonSteps, samplingPeriod, simulationDuration, inputCostMultiplier, initialState, inputMin, inputMax, stateMin, stateMax, linCombStateBounds, linCombInputBounds, eqCombInputBounds, timeVaryingStateFunc, timeVaryingStateBound, timeVaryingSteps, stateBound, inputBound, eqInputBound, generatePlots)
    warning off all
    
    numSteps = floor(simulationDuration / samplingPeriod);
    [numStates, numInputs] = size(B);
    
    % Dual Mode Matrices
    stateCostMatrix = C' * C;
    inputCostMatrix = inputCostMultiplier * eye(length(D));
    
    stateTransitionMatrix = zeros(numStates * horizonSteps, numInputs * horizonSteps);
    tempMatrix = [B zeros(numStates, numInputs * (horizonSteps - 1))];
    
    for i = 1:horizonSteps
        stateTransitionMatrix((i-1) * numStates + 1:i * numStates, :) = tempMatrix;
        tempMatrix = [A^(i) * B tempMatrix(:, 1:numInputs * (horizonSteps - 1))];
    end
    stateTransitionMatrix = [zeros(numStates, numInputs * horizonSteps); stateTransitionMatrix];
    clear tempMatrix
    
    controlMatrix = zeros(numStates * (horizonSteps + 1), numStates);
    for i = 1:horizonSteps + 1
        controlMatrix((i-1) * numStates + 1:i * numStates, :) = A^(i-1);
    end
    
    stateCostTemp = stateCostMatrix;
    for i = 1:horizonSteps - 1
        stateCostTemp = blkdiag(stateCostTemp, stateCostMatrix);
    end
    [steadyStateCostMatrix, ~, ~] = idare(A, B, stateCostMatrix, inputCostMatrix, [], []);
    steadyStateFeedback = (inputCostMatrix + B' * steadyStateCostMatrix * B) \ (B' * steadyStateCostMatrix * A);
    
    statePenaltyMatrix = dlyap((A - B * steadyStateFeedback)', stateCostMatrix + steadyStateFeedback' * inputCostMatrix * steadyStateFeedback);
    terminalCostMatrix = statePenaltyMatrix;
    
    totalCostMatrix = blkdiag(stateCostTemp, terminalCostMatrix);
    clear stateCostTemp
    
    inputCostTemp = inputCostMatrix;
    for i = 1:horizonSteps - 2
        inputCostTemp = blkdiag(inputCostTemp, inputCostMatrix);
    end
    totalInputCostMatrix = blkdiag(inputCostTemp, inputCostMatrix);
    clear inputCostTemp
    
    quadraticCostMatrix = stateTransitionMatrix' * totalCostMatrix * stateTransitionMatrix + totalInputCostMatrix;
    
    stateTrajectories = zeros(numStates, numSteps + 1);
    stateTrajectories(:, 1) = initialState;
    inputTrajectories = zeros(numInputs, numSteps);
    outputTrajectories(:, 1) = C * initialState;
    
    % Time Invariant Constraints
    if ~isempty(inputMin)
        inputMinMatrix = repmat(inputMin, horizonSteps, 1);
    else
        inputMinMatrix = [];
    end
    
    if ~isempty(inputMax)
        inputMaxMatrix = repmat(inputMax, horizonSteps, 1);
    else
        inputMaxMatrix = [];
    end
    
    constraintMatrix = [];
    if ~isempty(stateMin)
        stateMinMatrix = repmat(stateMin, horizonSteps + 1, 1);
        constraintMatrix = [constraintMatrix; -stateTransitionMatrix];
    end
    
    if ~isempty(stateMax)
        stateMaxMatrix = repmat(stateMax, horizonSteps + 1, 1);
        constraintMatrix = [constraintMatrix; stateTransitionMatrix];
    end
    
    if ~isempty(stateBound)
        stateBoundMatrix = [];
        for i = 1:horizonSteps + 1
            stateBoundMatrix = blkdiag(stateBoundMatrix, linCombStateBounds');
        end
        stateBoundValues = repmat(stateBound, horizonSteps + 1, 1);
        constraintMatrix = [constraintMatrix; stateBoundMatrix' * stateTransitionMatrix];
    end
    
    if ~isempty(inputBound)
        inputBoundMatrix = [];
        for i = 1:horizonSteps
            inputBoundMatrix = blkdiag(inputBoundMatrix, linCombInputBounds');
        end
        inputBoundValues = repmat(inputBound, horizonSteps, 1);
        constraintMatrix = [constraintMatrix; inputBoundMatrix'];
    else
        inputBoundValues = [];
    end
    
    if ~isempty(eqInputBound)
        eqInputBoundMatrix = [];
        for i = 1:horizonSteps
            eqInputBoundMatrix = blkdiag(eqInputBoundMatrix, eqCombInputBounds');
        end
        eqInputBoundValues = repmat(eqInputBound, horizonSteps, 1);
        eqConstraintMatrix = [eqInputBoundMatrix'];
    else
        eqInputBoundValues = [];
        eqConstraintMatrix = [];
    end
    
    eqConstraintValues = [eqInputBoundValues];
    
    for i = 1:numSteps
        %% Time Varying Constraints
        if ~isempty(stateMin)
            currentStateMin = -stateMinMatrix + controlMatrix * stateTrajectories(:, i);
        else
            currentStateMin = [];
        end
        
        if ~isempty(stateMax)
            currentStateMax = stateMaxMatrix - controlMatrix * stateTrajectories(:, i);
        else
            currentStateMax = [];
        end
        
        if ~isempty(stateBound)
            currentStateBound = stateBoundValues - stateBoundMatrix' * controlMatrix * stateTrajectories(:, i);
        else
            currentStateBound = [];
        end
        
        if ~isempty(timeVaryingSteps)
            varyingConstraintMatrix = [];
            varyingConstraintValues = [];
            for stepIdx = 1:length(timeVaryingSteps)
                stepMatches = ismember(i:i + horizonSteps - 1, timeVaryingSteps{stepIdx});
                varyingStateMatrix = [];
                varyingStateMatrix = blkdiag(varyingStateMatrix, zeros(numStates, 1));
                for match = stepMatches
                    if match
                        vec = timeVaryingStateFunc(stepIdx, :)';
                    else
                        vec = zeros(numStates, 1);
                    end
                    varyingStateMatrix = blkdiag(varyingStateMatrix, vec);
                end
                varyingConstraintMatrix = [varyingConstraintMatrix; (varyingStateMatrix' * stateTransitionMatrix)];
                timeVaryingBoundTemp = [0 stepMatches] .* timeVaryingStateBound(stepIdx);
                timeVaryingBoundTemp = reshape(timeVaryingBoundTemp, [], 1);
                varyingConstraintValues = [varyingConstraintValues; timeVaryingBoundTemp - varyingStateMatrix' * controlMatrix * stateTrajectories(:, i)];
            end
        else
            varyingConstraintMatrix = constraintMatrix;
            varyingConstraintValues = currentStateBound;
        end
        varyingConstraintMatrix = [constraintMatrix; varyingConstraintMatrix];
        
        constraints = [currentStateMin; currentStateMax; currentStateBound; inputBoundValues; varyingConstraintValues];
        
        quadraticCostVector = (stateTrajectories(:, i)' * controlMatrix' * totalCostMatrix * stateTransitionMatrix)';
        
        % Quadratic Programming Cost Optimization 
        optimalInputs = quadprog(quadraticCostMatrix, quadraticCostVector', varyingConstraintMatrix, constraints, eqConstraintMatrix, eqConstraintValues, inputMinMatrix, inputMaxMatrix, [], optimset("Display", "off"));
        inputTrajectories(:, i) = optimalInputs(1:numInputs);
        
        % State Update
        stateTrajectories(:, i + 1) = A * stateTrajectories(:, i) + B * inputTrajectories(:, i);
        outputTrajectories(:, i + 1) = C * stateTrajectories(:, i + 1) + D * inputTrajectories(:, 1);
    end
    
    warning on
    
    % Plot Generation
    if (generatePlots)
        figure;
        subplot(3, 1, 1)
        stairs(0:samplingPeriod:simulationDuration, outputTrajectories, 'DisplayName', 'Output')
        title("MPC DualMode | Horizon Steps = " + horizonSteps + " & \rho = " + inputCostMultiplier + " | Output")
        xlabel("t_k")
        ylabel("y_k")
        subplot(3, 1, 2)
        stairs(samplingPeriod:samplingPeriod:simulationDuration, inputTrajectories')
        title("MPC DualMode | Horizon Steps = " + horizonSteps + " & \rho = " + inputCostMultiplier + " | Input")
        xlabel("t_k")
        ylabel("u_k")
        legend(arrayfun(@(idx) sprintf('u_{%d,k}', idx), 1:numInputs, 'UniformOutput', false));
        subplot(3, 1, 3)
        stairs(0:samplingPeriod:simulationDuration, stateTrajectories')
        title("MPC DualMode | Horizon Steps = " + horizonSteps + " & \rho = " + inputCostMultiplier + " | States")
        xlabel("t_k")
        ylabel("x_k")
        legend(arrayfun(@(idx) sprintf('x_{%d,k}', idx), 1:numStates, 'UniformOutput', false));
    end
end
function [A_d, B_d, C_d, D_d] = DSS(A_cont, B_cont, C_cont, D_cont, sampling_time)
    % DiscretizeStateSpace: Discretizes a continuous-time state-space model.

    % Discretize A matrix using matrix exponential
    A_d = expm(A_cont * sampling_time);

    % Discretize B matrix using numerical integration
    % Define function handle for integration
    phi = @(tau) expm(A_cont * tau);
    
    % Numerical integration over one sampling period
    B_integral = integral(@(tau) phi(tau), 0, sampling_time, 'ArrayValued', true);
    B_d = B_integral * B_cont;

    % Directly transfer C and D matrices
    C_d = C_cont;
    D_d = D_cont;
end
function plotObstacles(Obstacles)
    durations = vertcat(Obstacles.duration);
    widths = vertcat(Obstacles.width);
    
    positions = [min(durations, [], 2) min(widths, [], 2) abs(durations(:,2) - durations(:,1)) abs(widths(:,2) - widths(:,1))];
    
    arrayfun(@(i) rectangle('Position', positions(i, :), 'Curvature', 0.3), 1:size(positions, 1));
end
