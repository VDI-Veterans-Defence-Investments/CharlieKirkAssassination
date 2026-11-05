function plotProjectedTrajectoryFrames( ...
    trajectory, ...
    velocityHistory, ...
    contained, ...
    numberOfSteps, ...
    projection, ...
    mesh, ...
    dt)
%PLOTPROJECTEDTRAJECTORYFRAMES
%
% Creates a 2D projected point map and arrow field for every frame.
%
% The projection plane is perpendicular to the manifold main axis.

    nTrajectories = size(trajectory,1);
    nFrames = size(trajectory,3);

    % Project the boundary into the main-axis-normal plane.
    boundary2D = projectPoints2D(mesh.V,projection);

    figure('Color','w');

    for frame = 1:nFrames

        clf;
        hold on;
        axis equal;
        grid on;

        % Project manifold boundary vertices.
        scatter( ...
            boundary2D(:,1), ...
            boundary2D(:,2), ...
            3, ...
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

            q = projectPoints2D(p(:)',projection);

            currentPoints = [currentPoints;q]; %#ok<AGROW>

            if any(~isfinite(v))
                v = [0 0 0];
            end

            % Project the velocity into the same 2D plane.
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

            % Draw projected trajectory positions as a point map.
            scatter( ...
                currentPoints(:,1), ...
                currentPoints(:,2), ...
                22, ...
                currentColors, ...
                'filled');

            % Draw projected arrows.
            arrowScale = 0.08;

            quiver( ...
                currentPoints(:,1), ...
                currentPoints(:,2), ...
                arrowScale*currentVelocities(:,1), ...
                arrowScale*currentVelocities(:,2), ...
                0, ...
                'Color',[0.15 0.15 0.15], ...
                'LineWidth',0.6, ...
                'MaxHeadSize',1.8);
        end

        xlabel('Projection axis e_1');
        ylabel('Projection axis e_2');

        title(sprintf( ...
            'Frame %d / %d   Time = %.3f', ...
            frame,nFrames,(frame-1)*dt));

        legend( ...
            'Projected manifold', ...
            'Trajectory point map', ...
            'Trajectory arrows', ...
            'Location','best');

        drawnow;

        % Optional animation speed.
        pause(0.02);
    end
end
