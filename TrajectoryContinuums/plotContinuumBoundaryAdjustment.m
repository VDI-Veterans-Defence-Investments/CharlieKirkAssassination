function [fig, VAdjusted] = plotContinuumBoundaryAdjustment( ...
    V, F, externalBoundaryFcn, parameters)
%PLOTCONTINUUMBOUNDARYADJUSTMENT
% Plot a continuum manifold boundary before and after adjustment.
%
% Inputs
%   V                  Original vertices, N-by-3
%   F                  Triangle connectivity, M-by-3
%   externalBoundaryFcn
%                      Function handle with the form:
%                      VAdjusted = externalBoundaryFcn(V, parameters)
%   parameters         Structure passed to externalBoundaryFcn
%
% Outputs
%   fig                Figure handle
%   VAdjusted          Adjusted vertices
%
% Example:
%   parameters.radius = 1.0;
%   parameters.zTop = 1.0;
%   parameters.zBottom = -4.0;
%   parameters.twistAmount = 0.35;
%   parameters.radialWave = 0.08;
%   parameters.xScale = 1.0;
%   parameters.yScale = 1.0;
%
%   [fig,VAdjusted] = plotContinuumBoundaryAdjustment( ...
%       V, F, @reformulateBoundary, parameters);

    if nargin < 4
        parameters = struct();
    end

    validateattributes(V, {'numeric'}, ...
        {'2d','ncols',3,'finite','real'}, mfilename, 'V');

    validateattributes(F, {'numeric'}, ...
        {'2d','ncols',3,'integer','positive'}, mfilename, 'F');

    if ~isa(externalBoundaryFcn, 'function_handle')
        error('externalBoundaryFcn must be a function handle.');
    end

    if any(F(:) > size(V,1))
        error('F contains vertex indices outside the range of V.');
    end

    % Apply the externally defined boundary transformation.
    VAdjusted = externalBoundaryFcn(V, parameters);

    if ~isequal(size(VAdjusted), size(V))
        error('The external function must return an N-by-3 vertex array.');
    end

    if any(~isfinite(VAdjusted(:)))
        error('The adjusted boundary contains NaN or Inf values.');
    end

    % Use common plotting limits so both geometries are comparable.
    allVertices = [V; VAdjusted];

    xyzMin = min(allVertices, [], 1);
    xyzMax = max(allVertices, [], 1);

    center = 0.5*(xyzMin + xyzMax);
    halfRange = 0.5*max(xyzMax - xyzMin);

    if halfRange <= 0
        halfRange = 1;
    end

    plotMin = center - 1.10*halfRange;
    plotMax = center + 1.10*halfRange;

    % Create figure and tiled layout.
    fig = figure( ...
        'Color', 'w', ...
        'Name', 'Continuum Manifold Boundary Adjustment');

    tiledlayout(fig, 1, 2, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    % Original boundary.
    ax1 = nexttile;
    patch(ax1, ...
        'Vertices', V, ...
        'Faces', F, ...
        'FaceColor', [0.25 0.55 0.85], ...
        'EdgeColor', [0.15 0.15 0.15], ...
        'FaceAlpha', 0.80, ...
        'EdgeAlpha', 0.25);

    configureBoundaryAxes(ax1, plotMin, plotMax);
    title(ax1, 'Original continuum boundary');
    xlabel(ax1, 'x');
    ylabel(ax1, 'y');
    zlabel(ax1, 'z');

    % Adjusted boundary.
    ax2 = nexttile;
    patch(ax2, ...
        'Vertices', VAdjusted, ...
        'Faces', F, ...
        'FaceColor', [0.90 0.35 0.20], ...
        'EdgeColor', [0.15 0.15 0.15], ...
        'FaceAlpha', 0.80, ...
        'EdgeAlpha', 0.25);

    configureBoundaryAxes(ax2, plotMin, plotMax);
    title(ax2, 'Externally adjusted boundary');
    xlabel(ax2, 'x');
    ylabel(ax2, 'y');
    zlabel(ax2, 'z');

    % Synchronize camera interaction between the two plots.
    linkprop([ax1 ax2], {'CameraPosition', ...
                         'CameraTarget', ...
                         'CameraUpVector', ...
                         'CameraViewAngle'});

    % Add an overall title.
    sgtitle(fig, 'Continuum Manifold Boundary Before and After Adjustment');
end


function configureBoundaryAxes(ax, plotMin, plotMax)
% Configure common 3D plotting properties.

    hold(ax, 'on');
    axis(ax, 'equal');
    axis(ax, [plotMin(1) plotMax(1) ...
              plotMin(2) plotMax(2) ...
              plotMin(3) plotMax(3)]);

    grid(ax, 'on');
    box(ax, 'on');
    view(ax, 3);
    camlight(ax, 'headlight');
    lighting(ax, 'gouraud');
end

