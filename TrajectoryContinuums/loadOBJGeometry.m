function [V, F, mesh] = loadOBJGeometry(filename)
%LOADOBJGEOMETRY Load a triangulated geometry from an OBJ file.
%
%   [V,F] = loadOBJGeometry(filename)
%   [V,F,mesh] = loadOBJGeometry(filename)
%
% Outputs
%   V     N-by-3 vertex coordinates
%   F     M-by-3 triangular face indices
%   mesh  Optional mesh structure containing vertices, faces, triangle data,
%         and AABB acceleration data
%
% Supported OBJ elements:
%   v x y z              Vertex positions
%   f v1 v2 v3           Triangular faces
%   f v1 v2 v3 v4        Polygon faces, triangulated automatically
%   f v/vt/vn            Position/texture/normal format
%   f v//vn              Position/normal format
%   Negative vertex indices are supported.
%
% Example:
%   [V,F,mesh] = loadOBJGeometry('myGeometry.obj');
%
%   % Use with the trajectory code:
%   Vreformulated = reformulateBoundary(V, boundaryParameters);
%   mesh = buildMeshAccelerationData(Vreformulated, F);

    if nargin ~= 1
        error('Usage: [V,F,mesh] = loadOBJGeometry(filename)');
    end

    if ~(ischar(filename) || isstring(filename))
        error('filename must be a character vector or string.');
    end

    filename = char(filename);

    if ~isfile(filename)
        error('OBJ file not found: %s', filename);
    end

    fileID = fopen(filename, 'r');

    if fileID == -1
        error('Could not open OBJ file: %s', filename);
    end

    cleanupObject = onCleanup(@() fclose(fileID));

    vertices = zeros(0,3);
    faces = zeros(0,3);

    lineNumber = 0;

    while true

        line = fgetl(fileID);

        if ~ischar(line)
            break;
        end

        lineNumber = lineNumber + 1;

        line = strtrim(line);

        % Ignore blank lines and comments.
        if isempty(line) || startsWith(line, '#')
            continue;
        end

        % Remove inline comments.
        commentIndex = strfind(line, '#');

        if ~isempty(commentIndex)
            line = strtrim(line(1:commentIndex(1)-1));
        end

        tokens = strsplit(line);

        if isempty(tokens)
            continue;
        end

        recordType = tokens{1};

        switch recordType

            case 'v'
                % Vertex position.
                if numel(tokens) < 4
                    error('Invalid vertex on line %d.', lineNumber);
                end

                vertex = str2double(tokens(2:4));

                if any(~isfinite(vertex))
                    error('Invalid vertex coordinates on line %d.', ...
                        lineNumber);
                end

                vertices(end+1,:) = vertex; %#ok<AGROW>

            case 'f'
                % Face may contain:
                % f v1 v2 v3
                % f v1/vt1 v2/vt2 v3/vt3
                % f v1//vn1 v2//vn2 v3//vn3
                if numel(tokens) < 4
                    error('Face has fewer than three vertices on line %d.', ...
                        lineNumber);
                end

                faceTokens = tokens(2:end);
                faceIndices = zeros(1, numel(faceTokens));

                for k = 1:numel(faceTokens)

                    components = strsplit(faceTokens{k}, '/');

                    vertexIndex = str2double(components{1});

                    if ~isfinite(vertexIndex) || vertexIndex == 0 || ...
                            vertexIndex ~= fix(vertexIndex)
                        error('Invalid face index on line %d.', ...
                            lineNumber);
                    end

                    % OBJ negative indices refer backward from the current
                    % last vertex.
                    if vertexIndex < 0
                        vertexIndex = size(vertices,1) + vertexIndex + 1;
                    end

                    if vertexIndex < 1 || ...
                            vertexIndex > size(vertices,1)
                        error('Face index out of range on line %d.', ...
                            lineNumber);
                    end

                    faceIndices(k) = vertexIndex;
                end

                % Fan triangulation:
                % [1 2 3 4] becomes [1 2 3] and [1 3 4].
                for k = 2:(numel(faceIndices)-1)
                    faces(end+1,:) = [ ...
                        faceIndices(1), ...
                        faceIndices(k), ...
                        faceIndices(k+1)]; %#ok<AGROW>
                end

            otherwise
                % Ignore OBJ records such as:
                % vt, vn, g, usemtl, mtllib, o, s, and similar.
                continue;
        end
    end

    V = vertices;
    F = faces;

    if isempty(V)
        error('The OBJ file contains no vertex positions.');
    end

    if isempty(F)
        error('The OBJ file contains no faces.');
    end

    % Optional mesh structure compatible with pointInsideMeshAABB.
    if nargout >= 3
        mesh = buildMeshAccelerationDataLocal(V, F);
    end
end


function mesh = buildMeshAccelerationDataLocal(V, F)
% Build triangle and axis-aligned bounding-box data.

    A = V(F(:,1),:);
    B = V(F(:,2),:);
    C = V(F(:,3),:);

    mesh.V = V;
    mesh.F = F;

    mesh.A = A;
    mesh.B = B;
    mesh.C = C;

    mesh.faceAABBMin = min(cat(3, A, B, C), [], 3);
    mesh.faceAABBMax = max(cat(3, A, B, C), [], 3);

    mesh.globalAABBMin = min(V, [], 1);
    mesh.globalAABBMax = max(V, [], 1);
end

