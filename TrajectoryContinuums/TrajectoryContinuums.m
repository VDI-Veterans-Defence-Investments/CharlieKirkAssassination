function trajectoryContinuumFrameMap
% Complete trajectory-continuum simulation in one MATLAB file.
%
% The file includes:
%   - Boundary manifold generation
%   - External boundary reformulation
%   - AABB acceleration
%   - Ray-based mesh containment
%   - Frame-by-frame trajectory propagation
%   - Parallel trajectory updates
%   - PCA projection
%   - 2D point-map and arrow animation

    clear;
    close all;
    clc;

    rng(7);

    %% Parameters

    nTrajectories = 500;
    maxFrames = 500;

    START.z = 1.0;
    END.z = -4.0;

    radius = 1.0;
    dt = 0.015;

    baseMargin = 1e-8;
    aabbScale = 1.5;
    maxMargin = 0.20;
    rayTolerance = 1e-10;

    %% Create original closed manifold

    nAround = 64;
    nHeight = 30;

    [V,F] = makeClosedCylinder( ...
        radius, ...
        END.z, ...
        START.z, ...
        nAround, ...
        nHeight);

    %% Reformulate boundary using external equation

    boundaryParameters.radius = radius;
    boundaryParameters.zTop = START.z;
    boundaryParameters.zBottom = END.z;

    boundaryParameters.twistAmount = 0.35;
    boundaryParameters.radialWave = 0.08;
    boundaryParameters.xScale = 1.00;
    boundaryParameters.yScale = 1.00;

    Vreformulated = reformulateBoundary( ...
        V,boundaryParameters);

    %% Build AABB acceleration data

    mesh = buildMeshAccelerationData(Vreformulated,F);

    %% Generate random START points

    P0 = generateStartPoints( ...
        nTrajectories, ...
        radius, ...
        START.z);

    %% Build main-axis projection

    projection = buildMainAxisProjection(mesh.V);

    %% Initialize parallel processing

    if isempty(gcp('nocreate'))
        parpool;
    end

    %% Run frame-by-frame simulation

    [trajectory,velocityHistory,contained,numberOfSteps] = ...
        simulateFrameByFrame( ...
            P0, ...
            mesh, ...
            maxFrames, ...
            dt, ...
            END.z, ...
            baseMargin, ...
            aabbScale, ...
            maxMargin, ...
            rayTolerance);

    %% Print results

    fprintf('\nTrajectory continuum result\n');
    fprintf('---------------------------\n');
    fprintf('Trajectories:       %d\n',nTrajectories);
    fprintf('Contained:          %d\n',sum(contained));
    fprintf('Escaped:            %d\n',sum(~contained));
    fprintf('Containment ratio:  %.4f\n',mean(contained));

    %% Animate 2D projected point map

    plotProjectedTrajectoryFrames( ...
        trajectory, ...
        velocityHistory, ...
        contained, ...
        numberOfSteps, ...
        projection, ...
        mesh, ...
        dt);
end


function P0 = generateStartPoints(n,radius,zStart)
% Generate random points inside the circular START cross-section.

    P0 = zeros(n,3);
    acceptedCount = 0;

    while acceptedCount < n

        candidates = 2*radius*(rand(2*n,2)-0.5);

        valid = sum(candidates.^2,2) <= radius^2;
        candidates = candidates(valid,:);

        numberToTake = min(size(candidates,1), ...
            n-acceptedCount);

        P0(acceptedCount+1:acceptedCount+numberToTake,1:2) = ...
            candidates(1:numberToTake,:);

        acceptedCount = acceptedCount+numberToTake;
    end

    P0(:,3) = zStart;
end


function Vnew = reformulateBoundary(V,param)
% Apply an external equation to the boundary geometry.
%
% The example equation is:
%
%   r'     = r*(1+a*sin(4*theta+2*pi*zeta))
%   theta' = theta+twist*zeta
%
% where zeta is normalized height.

    x = V(:,1);
    y = V(:,2);
    z = V(:,3);

    theta = atan2(y,x);
    r = sqrt(x.^2+y.^2);

    zeta = (z-param.zBottom)/ ...
        (param.zTop-param.zBottom);

    zeta = max(0,min(1,zeta));

    radialFactor = 1+param.radialWave* ...
        sin(4*theta+2*pi*zeta);

    rNew = r.*radialFactor;

    thetaNew = theta+param.twistAmount*zeta;

    xNew = param.xScale*rNew.*cos(thetaNew);
    yNew = param.yScale*rNew.*sin(thetaNew);
    zNew = z;

    Vnew = [xNew,yNew,zNew];
end


function mesh = buildMeshAccelerationData(V,F)
% Precompute triangle data and AABBs.

    A = V(F(:,1),:);
    B = V(F(:,2),:);
    C = V(F(:,3),:);

    mesh.V = V;
    mesh.F = F;

    mesh.A = A;
    mesh.B = B;
    mesh.C = C;

    mesh.faceAABBMin = min(cat(3,A,B,C),[],3);
    mesh.faceAABBMax = max(cat(3,A,B,C),[],3);

    mesh.globalAABBMin = min(V,[],1);
    mesh.globalAABBMax = max(V,[],1);
end


function [trajectory,velocityHistory,contained,numberOfSteps] = ...
    simulateFrameByFrame( ...
        P0, ...
        mesh, ...
        maxFrames, ...
        dt, ...
        endZ, ...
        baseMargin, ...
        aabbScale, ...
        maxMargin, ...
        rayTolerance)
% Advance all trajectories one frame at a time.
%
% Trajectories are independent and update in parallel.
% Each individual trajectory advances sequentially.

    nTrajectories = size(P0,1);

    trajectory = nan(nTrajectories,3,maxFrames);
    velocityHistory = nan(nTrajectories,3,maxFrames);

    trajectory(:,:,1) = P0;

    contained = true(nTrajectories,1);
    active = true(nTrajectories,1);
    numberOfSteps = ones(nTrajectories,1);

    for frame = 2:maxFrames

        activeIndices = find(active);

        if isempty(activeIndices)
            break;
        end

        nActive = numel(activeIndices);

        nextPositions = nan(nActive,3);
        nextVelocities = nan(nActive,3);
        nextInside = false(nActive,1);
        nextActive = false(nActive,1);

        parfor j = 1:nActive

            trajectoryIndex = activeIndices(j);

            previousPoint = squeeze( ...
                trajectory(trajectoryIndex,:,frame-1));

            currentTime = (frame-1)*dt;

            velocity = trajectoryVelocity( ...
                previousPoint,currentTime);

            velocity = velocity(:)';
            velocityMagnitude = norm(velocity);

            if velocityMagnitude < eps
                nextPositions(j,:) = previousPoint;
                nextVelocities(j,:) = [0 0 0];
                continue;
            end

            direction = velocity/velocityMagnitude;
            nextPoint = previousPoint+dt*direction;

            distanceMoved = norm(nextPoint-previousPoint);

            margin = min( ...
                maxMargin, ...
                baseMargin+aabbScale*distanceMoved);

            isInside = pointInsideMeshAABB( ...
                nextPoint, ...
                mesh, ...
                margin, ...
                rayTolerance);

            nextPositions(j,:) = nextPoint;
            nextVelocities(j,:) = direction;
            nextInside(j) = isInside;
            nextActive(j) = isInside && nextPoint(3) > endZ;
        end

        % Commit the completed frame.
        for j = 1:nActive

            trajectoryIndex = activeIndices(j);

            trajectory(trajectoryIndex,:,frame) = ...
                nextPositions(j,:);

            velocityHistory(trajectoryIndex,:,frame) = ...
                nextVelocities(j,:);

            contained(trajectoryIndex) = ...
                contained(trajectoryIndex) && nextInside(j);

            numberOfSteps(trajectoryIndex) = frame;
            active(trajectoryIndex) = nextActive(j);
        end
    end
end


function velocity = trajectoryVelocity(p,t)
% Define the trajectory arrow/vector field here.

    p = p(:)';

    velocity = [ ...
        0.035*sin(1.4*p(2))+0.012*cos(t), ...
        0.035*cos(1.2*p(1))+0.012*sin(t), ...
       -1.0];
end


function inside = pointInsideMeshAABB(P,mesh,margin,rayTolerance)
% Test point containment using AABB filtering and ray casting.

    globalMin = mesh.globalAABBMin-margin;
    globalMax = mesh.globalAABBMax+margin;

    if any(P < globalMin) || any(P > globalMax)
        inside = false;
        return;
    end

    rayEndX = globalMax(1)+2*margin+1;

    rayMin = [P(1),P(2)-margin,P(3)-margin];
    rayMax = [rayEndX,P(2)+margin,P(3)+margin];

    candidates = ...
        mesh.faceAABBMin(:,1) <= rayMax(1) & ...
        mesh.faceAABBMax(:,1) >= rayMin(1) & ...
        mesh.faceAABBMin(:,2) <= rayMax(2) & ...
        mesh.faceAABBMax(:,2) >= rayMin(2) & ...
        mesh.faceAABBMin(:,3) <= rayMax(3) & ...
        mesh.faceAABBMax(:,3) >= rayMin(3);

    candidateIndices = find(candidates);

    if isempty(candidateIndices)
        inside = false;
        return;
    end

    rayOrigin = P(:)';
    rayOrigin(2) = rayOrigin(2)+1e-11;
    rayOrigin(3) = rayOrigin(3)+2e-11;

    rayDirection = [1 0 0];
    crossings = 0;

    for q = 1:numel(candidateIndices)

        f = candidateIndices(q);

        A = mesh.A(f,:);
        B = mesh.B(f,:);
        C = mesh.C(f,:);

        edge1 = B-A;
        edge2 = C-A;

        h = cross(rayDirection,edge2);
        determinant = dot(edge1,h);

        if abs(determinant) < rayTolerance
            continue;
        end

        inverseDeterminant = 1/determinant;
        s = rayOrigin-A;

        u = inverseDeterminant*dot(s,h);

        if u < -rayTolerance || u > 1+rayTolerance
            continue;
        end

        qVector = cross(s,edge1);
        v = inverseDeterminant*dot(rayDirection,qVector);

        if v < -rayTolerance || u+v > 1+rayTolerance
            continue;
        end

        rayDistance = inverseDeterminant*dot(edge2,qVector);

        if rayDistance > rayTolerance
            crossings = crossings+1;
        end
    end

    inside = mod(crossings,2) == 1;
end


function projection = buildMainAxisProjection(V)
% Use PCA to find the main manifold axis and 2D projection plane.

    origin = mean(V,1);
    X = V-origin;

    covarianceMatrix = (X'*X)/size(X,1);

    [eigenvectors,eigenvalues] = eig(covarianceMatrix);
    [~,order] = sort(diag(eigenvalues),'descend');

    mainAxis = eigenvectors(:,order(1));
    mainAxis = mainAxis/norm(mainAxis);

    if mainAxis(3) < 0
        mainAxis = -mainAxis;
    end

    reference = [0;0;1];

    if abs(dot(reference,mainAxis)) > 0.90
        reference = [1;0;0];
    end

    e1 = cross(mainAxis,reference);
    e1 = e1/norm(e1);

    e2 = cross(mainAxis,e1);
    e2 = e2/norm(e2);

    projection.origin = origin(:);
    projection.mainAxis = mainAxis;
    projection.e1 = e1;
    projection.e2 = e2;
end


function Q = projectPoints2D(P,projection)
% Project 3D points into the plane normal to the main manifold axis.

    if isvector(P)
        P = reshape(P,1,3);
    end

    relative = P-projection.origin';

    Q = [ ...
        relative*projection.e1, ...
        relative*projection.e2];
end


function plotProjectedTrajectoryFrames( ...
    trajectory, ...
    velocityHistory, ...
    contained, ...
    numberOfSteps, ...
    projection, ...
    mesh, ...
    dt)
% Animate the projected point map and arrow field.

    nTrajectories = size(trajectory,1);
    nFrames = size(trajectory,3);

    boundary2D = projectPoints2D(mesh.V,projection);

    figure('Color','w');

    for frame = 1:nFrames

        clf;
        hold on;
        axis equal;
        grid on;

        scatter( ...
            boundary2D(:,1), ...
            boundary2D(:,2), ...
            4, ...
            [0.65 0.75 0.85], ...
            'filled');

        currentPoints = [];
        currentVelocities = [];
        currentColors = [];

        for i = 1:nTrajectories

            if frame > numberOfSteps(i)
                continue;
            end

            p = squeeze(trajectory(i,:,frame));
            v = squeeze(velocityHistory(i,:,frame));

            if any(~isfinite(p))
                continue;
            end

            q = projectPoints2D(p,projection);

            currentPoints = [currentPoints;q]; %#ok<AGROW>

            if any(~isfinite(v))
                v = [0 0 0];
            end

            velocity2D = [ ...
                dot(v,projection.e1), ...
                dot(v,projection.e2)];

            currentVelocities = ...
                [currentVelocities;velocity2D]; %#ok<AGROW>

            if contained(i)
                currentColors = ...
                    [currentColors;0.00 0.65 0.15]; %#ok<AGROW>
            else
                currentColors = ...
                    [currentColors;0.85 0.05 0.02]; %#ok<AGROW>
            end
        end

        if ~isempty(currentPoints)

            scatter( ...
                currentPoints(:,1), ...
                currentPoints(:,2), ...
                24, ...
                currentColors, ...
                'filled');

            arrowScale = 0.12;

            quiver( ...
                currentPoints(:,1), ...
                currentPoints(:,2), ...
                arrowScale*currentVelocities(:,1), ...
                arrowScale*currentVelocities(:,2), ...
                0, ...
                'Color',[0.10 0.10 0.10], ...
                'LineWidth',0.6, ...
                'MaxHeadSize',1.8);
        end

        xlabel('Projection axis e_1');
        ylabel('Projection axis e_2');

        title(sprintf( ...
            'Frame %d of %d    Time = %.3f', ...
            frame,nFrames,(frame-1)*dt));

        drawnow;
        pause(0.02);
    end
end


function [V,F] = makeClosedCylinder( ...
    radius,zBottom,zTop,nAround,nHeight)
% Create a closed triangulated cylinder.

    theta = linspace(0,2*pi,nAround+1);
    theta(end) = [];

    zLevels = linspace(zBottom,zTop,nHeight+1);

    V = zeros((nHeight+1)*nAround,3);

    for iz = 1:nHeight+1
        for it = 1:nAround

            index = (iz-1)*nAround+it;

            V(index,:) = [ ...
                radius*cos(theta(it)), ...
                radius*sin(theta(it)), ...
                zLevels(iz)];
        end
    end

    F = zeros(2*nHeight*nAround,3);
    face = 1;

    % Cylinder side.
    for iz = 1:nHeight
        for it = 1:nAround

            itNext = mod(it,nAround)+1;

            a = (iz-1)*nAround+it;
            b = (iz-1)*nAround+itNext;
            c = iz*nAround+itNext;
            d = iz*nAround+it;

            F(face,:) = [a b c];
            face = face+1;

            F(face,:) = [a c d];
            face = face+1;
        end
    end

    % Cap centers.
    bottomCenter = size(V,1)+1;
    topCenter = size(V,1)+2;

    V(bottomCenter,:) = [0 0 zBottom];
    V(topCenter,:) = [0 0 zTop];

    Fcaps = zeros(2*nAround,3);
    face = 1;

    % Bottom cap.
    for it = 1:nAround

        itNext = mod(it,nAround)+1;

        Fcaps(face,:) = [bottomCenter,itNext,it];
        face = face+1;
    end

    % Top cap.
    topOffset = nHeight*nAround;

    for it = 1:nAround

        itNext = mod(it,nAround)+1;

        a = topOffset+it;
        b = topOffset+itNext;

        Fcaps(face,:) = [topCenter,a,b];
        face = face+1;
    end

    F = [F;Fcaps];
end

