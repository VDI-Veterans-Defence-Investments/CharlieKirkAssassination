function [fig, transformedObjectVertices] = plotCombinedManifoldOBJ( ...
    manifoldVertices, manifoldFaces, objFilename, transformMatrix, options)
%PLOTCOMBINEDMANIFOLDOBJ Plot a trajectory manifold and transformed OBJ.
%
%   [fig,VobjWorld] = plotCombinedManifoldOBJ( ...
%       Vmanifold, Fmanifold, objFilename, T)
%
% Inputs
%   manifoldVertices
%       N-by-3 vertices of the trajectory manifold.
%
%   manifoldFaces
%       M-by-3 triangular face connectivity of the manifold.
%
%   objFilename
%       Path to the OBJ file.
%
%   transformMatrix
%       4-by-4 homogeneous transformation matrix:
%
%           [ R(3,3)  t(3,1) ]
%           [ 0 0 0      1   ]
%
%       The transformation is applied as:
%
%           pWorld = T * [pObject; 1]
%
% Optional input
%   options
%       Structure with optional fields:
%
%       options.ManifoldColor  Default: [0.20 0.55 0.90]
%       options.ObjectColor    Default: [0.90 0.35 0.15]
%       options.ManifoldAlpha  Default: 0.25
%       options.ObjectAlpha    Default: 0.85
%       options.ShowEdges      Default: true
%       options.Trajectory      Optional N-by-3-by-K trajectory array
%       options.TrajectoryColor Default: [0.05 0.05 0.05]
%       options.View            Default: [35 25]
%       options.Title           Default: 'Combined trajectory geometry'
%
% Outputs
%   fig
%       Figure handle.
%
%   transformedObjectVertices
%       Transformed OBJ vertices in world coordinates.
%
% Example
%   [Vobj,Fobj] = loadOBJGeometry('tool.obj');
%
%   T = makeTransform( ...
%       [1 0 0; 0 0 -1; 0 1 0], ...
%       [0.2 -0.1 -2.0]);
%
%   [fig,VobjWorld] = plotCombinedManifoldOBJ( ...
%       V, F, 'tool.obj', T);
%
% The OBJ loader must provide:
%
%   [Vobj,Fobj] = loadOBJGeometry(filename)

    if nargin < 5 || isempty(options)
        options = struct();
    end

    validateattributes(manifoldVertices, {'numeric'}, ...
        {'2d','ncols',3,'finite','real'}, ...
        mfilename, 'manifoldVertices');

    validateattributes(manifoldFaces, {'numeric'}, ...
        {'2d','ncols',3,'integer','positive'}, ...
        mfilename, 'manifoldFaces');

    validateattributes(transformMatrix, {'numeric'}, ...
        {'size',[4 4],'finite','real'}, ...
        mfilename, 'transformMatrix');

    if transformMatrix(4,4) == 0
        error('The homogeneous transform has an invalid final element.');
    end

    if any(manifoldFaces(:) > size(manifoldVertices,1))
        error('manifoldFaces contains invalid vertex indices.');
    end

    % Load OBJ geometry.
    [objectVertices, objectFaces] = loadOBJGeometry(objFilename);

    % Convert OBJ vertices to homogeneous coordinates.
    objectVerticesHomogeneous = [ ...
        objectVertices, ...
        ones(size(objectVertices,1),1)];

    % Apply the arbitrary homogeneous transformation.
    transformedHomogeneous = ...
        (transformMatrix*objectVerticesHomogeneous')';

    % Support a general projective homogeneous transform.
    homogeneousScale = transformedHomogeneous(:,4);

    if any(abs(homogeneousScale) < eps)
        error('The transform produces invalid homogeneous coordinates.');
    end

    transformedObjectVertices = ...
        transformedHomogeneous(:,1:3)./homogeneousScale;

    % Default plotting options.
    manifoldColor = getOption(options, ...
        'ManifoldColor', [0.20 0.55 0.90]);

    objectColor = getOption(options, ...
        'ObjectColor', [0.90 0.35 0.15]);

    manifoldAlpha = getOption(options, ...
        'ManifoldAlpha', 0.25);

    objectAlpha = getOption(options, ...
        'ObjectAlpha', 0.85);

    showEdges = getOption(options, ...
        'ShowEdges', true);

    trajectoryColor = getOption(options, ...
        'TrajectoryColor', [0.05 0.05 0.05]);

    viewDirection = getOption(options, ...
        'View', [35 25]);

    plotTitle = getOption(options, ...
        'Title', 'Combined trajectory geometry');

    if showEdges
        edgeColor = [0.15 0.15 0.15];
        edgeAlpha = 0.25;
    else
        edgeColor = 'none';
        edgeAlpha = 0;
    end

    % Create figure.
    fig = figure( ...
        'Color', 'w', ...
        'Name', 'Combined manifold and OBJ geometry');

    ax = axes(fig);
    hold(ax, 'on');

    % Plot the trajectory abstract manifold.
    patch(ax, ...
        'Vertices', manifoldVertices, ...
        'Faces', manifoldFaces, ...
        'FaceColor', manifoldColor, ...
        'FaceAlpha', manifoldAlpha, ...
        'EdgeColor', edgeColor, ...
        'EdgeAlpha', edgeAlpha, ...
        'LineStyle', '-');

    % Plot transformed OBJ geometry.
    patch(ax, ...
        'Vertices', transformedObjectVertices, ...
        'Faces', objectFaces, ...
        'FaceColor', objectColor, ...
        'FaceAlpha', objectAlpha, ...
        'EdgeColor', edgeColor, ...
        'EdgeAlpha', edgeAlpha, ...
        'LineStyle', '-');

    % Optionally plot trajectory paths and current locations.
    if isfield(options, 'Trajectory') && ...
            ~isempty(options.Trajectory)

        trajectory = options.Trajectory;

        if size(trajectory,2) ~= 3
            error('options.Trajectory must have size N-by-3-by-K.');
        end

        nTrajectories = size(trajectory,1);

        for i = 1:nTrajectories

            path = squeeze(trajectory(i,:,:))';

            valid = all(isfinite(path),2);

            if nnz(valid) < 2
                continue;
            end

            plot3(ax, ...
                path(valid,1), ...
                path(valid,2), ...
                path(valid,3), ...
                'Color', trajectoryColor, ...
                'LineWidth', 0.8);
        end
    end

    % Combine all points when determining plot limits.
    allPoints = [ ...
        manifoldVertices;
        transformedObjectVertices];

    if isfield(options, 'Trajectory') && ...
            ~isempty(options.Trajectory)

        trajectory = options.Trajectory;
        trajectoryPoints = reshape(trajectory, [], 3);
        trajectoryPoints = trajectoryPoints( ...
            all(isfinite(trajectoryPoints),2), :);

        allPoints = [allPoints; trajectoryPoints];
    end

    pointMin = min(allPoints, [], 1);
    pointMax = max(allPoints, [], 1);

    center = 0.5*(pointMin + pointMax);
    halfRange = 0.5*max(pointMax - pointMin);

    if halfRange <= 0
        halfRange = 1;
    end

    halfRange = 1.10*halfRange;

    axis(ax, [ ...
        center(1)-halfRange, center(1)+halfRange, ...
        center(2)-halfRange, center(2)+halfRange, ...
        center(3)-halfRange, center(3)+halfRange]);

    axis(ax, 'equal');
    grid(ax, 'on');
    box(ax, 'on');
    view(ax, viewDirection);

    xlabel(ax, 'x');
    ylabel(ax, 'y');
    zlabel(ax, 'z');

    title(ax, plotTitle);

    camlight(ax, 'headlight');
    lighting(ax, 'gouraud');

    legend(ax, ...
        {'Trajectory manifold', 'Transformed OBJ object'}, ...
        'Location', 'best');
end


function value = getOption(options, fieldName, defaultValue)
%GETOPTION Return an option value or its default.

    if isfield(options, fieldName) && ...
            ~isempty(options.(fieldName))

        value = options.(fieldName);
    else
        value = defaultValue;
    end
end

