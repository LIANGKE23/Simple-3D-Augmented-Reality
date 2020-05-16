%% Read in the file
file_path = 'COLMAP_Files\points3D.txt';
my_file = fopen(file_path);
my_data = textscan(my_file, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'CommentStyle','#');
fclose(my_file);


%% Organize data
my_data = cell2mat(my_data);
points = my_data(:, 2:4);


%% RANSAC: Find the largest subset of 3D points that can be described by the equation of a 3D plane
% parameters
min_pts = 3;
dist_thresh = .5;

% initialize
inlier_count = 0;
global_on_plane = [];
global_normal = [];

% begin
for i = 1:10000
    % Randomly select (min_pts) from the data
    sample_idxs = randsample(length(points(:,1)), min_pts);
    samples = points(sample_idxs, :);
    
    % Fit a plane using the samples above
    % Two vectors parallel to the plane
    p_vec_1 = samples(2, :) - samples(1, :);
    p_vec_2 = samples(3, :) - samples(1, :);
    
    % The normal vector (cross product of the two parallel vectors)
    normal = p_vec_1.*p_vec_2;
    
    % Plane: ax + by + cz + d = 0
    d = -sum(normal.*samples(1, :));
    
    % Determine the points that fall on the plane
    numerator = abs(sum(points.*normal(ones(length(points(:,1)),1), :),2) + d);
    denom = sqrt(normal(1)^2 + normal(2)^2 + normal(3)^2);
    dist = numerator/denom;
    
    on_plane = dist <= dist_thresh;
    on_plane_count = sum(on_plane);
    
    % Update
    if on_plane_count > inlier_count
        plane_pts = points(on_plane, :);
        global_on_plane = on_plane;
        inlier_count = on_plane_count;
        global_normal = normal;
    end
end


%% Display
% set colors
inlier_color = [0 1 0];
outlier_color = [1 0 0];
temp_in_colors = inlier_color(ones(length(points(:,1)),1), :);
temp_out_colors = outlier_color(ones(length(points(:,1)),1), :);

colors = temp_out_colors;
colors(global_on_plane, :) = temp_in_colors(global_on_plane, :);

% plot
sz = ones(length(points(:,1)), 1);
figure, scatter3(points(:,1), points(:,2), points(:,3), sz, colors)
xlim([-20, 50])
zlim([0, 100])
hold on

syms x y z
P = [x,y,z];
realdot = @(u, v) u*transpose(v);
planefunction = realdot(P-plane_pts(1,:), global_normal);
zplane = solve(planefunction, z);
ezmesh(zplane, [-20, 50, -20, 50]), hold off

zlim([0, 100])
ylim([-20, 15])

% x = points(global_on_plane, 1);
% y = points(global_on_plane, 2);
% z = points(global_on_plane, 3);
% x1 = points(~global_on_plane, 1);
% y1 = points(~global_on_plane, 2);
% z1 = points(~global_on_plane, 3);
% figure, plot3(x, y, z, '.', 'Color', 'g')
% hold on
% plot3(x1, y1, z1, '.', 'Color', 'r')
% hold off
% xlim([-20, 50])
% zlim([0, 100])


%% Geometric Transformation
% Translation: origin x=0, y=0 lies in center of inliers
center_x = (range(plane_pts(:,1))/2) + min(plane_pts(:,1));
center_y = (range(plane_pts(:,2))/2) + min(plane_pts(:,2));
center_z = subs(planefunction, x, center_x);
center_z = subs(center_z, y, center_y);
center_z = double(solve(center_z, z));

origin = -[center_x, center_y, center_z];
translation_matrix = repmat(origin, length(points(:,1)), 1);

% Rotation matrix: make dominant plane z = 0
z_vec = -global_normal;
z_axis = z_vec/norm(z_vec);
x_vec = cross([0, 1, 0], z_axis);
x_axis = x_vec/norm(x_vec);
y_vec = cross(z_axis, x_axis);
y_axis = y_vec/norm(y_vec);
R = [x_axis; y_axis; z_axis];

% Apply rotation and translation
trans_points = points + translation_matrix;
rotPoints = R*trans_points';
rotPoints = rotPoints';

% Plot
sz = ones(length(points(:,1)), 1);
figure, scatter3(rotPoints(:,1), rotPoints(:,2), rotPoints(:,3), sz, colors)
xlim([-50, 50])
ylim([-50, 50])
zlim([-50, 50])
hold on

[x1, y1] = meshgrid(-100:1:100); % Generate x and y data
z1 = zeros(length(x1(:,1)));
surf(x1,y1,z1, 'FaceAlpha', 0.5) %Plot the surface


%% Create the object and display on plane (simple 3D box)

% 3D box on the z=0 and centered at x=0, y=0, z=0
box_points = [2,2,0; -2,2,0; -2,-2,0; 2,-2,0; 2,2,3; -2,2,3; -2,-2,3; 2,-2,3];
fill3(box_points(1:4,1), box_points(1:4,2), box_points(1:4,3), 1); %bottom
fill3(box_points(5:8,1), box_points(5:8,2), box_points(5:8,3), 1); % top
x = [box_points(2,1) box_points(6,1) box_points(7,1) box_points(3,1)];
y = [box_points(2,2) box_points(6,2) box_points(7,2) box_points(3,2)];
z = [box_points(2,3) box_points(6,3) box_points(7,3) box_points(3,3)];
fill3(x, y, z, 3);
x = [box_points(2,1) box_points(6,1) box_points(5,1) box_points(1,1)];
y = [box_points(2,2) box_points(6,2) box_points(5,2) box_points(1,2)];
z = [box_points(2,3) box_points(6,3) box_points(5,3) box_points(1,3)];
fill3(x, y, z, 3);
x = [box_points(3,1) box_points(4,1) box_points(8,1) box_points(7,1)];
y = [box_points(3,2) box_points(4,2) box_points(8,2) box_points(7,2)];
z = [box_points(3,3) box_points(4,3) box_points(8,3) box_points(7,3)];
fill3(x, y, z, 3);
x = [box_points(1,1) box_points(4,1) box_points(8,1) box_points(5,1)];
y = [box_points(1,2) box_points(4,2) box_points(8,2) box_points(5,2)];
z = [box_points(1,3) box_points(4,3) box_points(8,3) box_points(5,3)];
fill3(x, y, z, 3);
hold off


%% Transform the object points back to scene coordinates
% rotate
rot_box_pts = R'*box_points';
rot_box_pts = rot_box_pts';

% translate
trans_box_pts = rot_box_pts - translation_matrix(1:numel(box_points(:,1)), :);

% display the box points in the original scene coordinates
figure, scatter3(points(:,1), points(:,2), points(:,3), sz, colors)
xlim([-20, 50])
zlim([0, 100])
hold on

syms x y z
P = [x,y,z];
realdot = @(u, v) u*transpose(v);
planefunction = realdot(P-plane_pts(1,:), global_normal);
zplane = solve(planefunction, z);
ezmesh(zplane, [-100, 100, -100, 100])

siz = 10*ones(length(rot_box_pts(:,1)), 1);
scatter3(trans_box_pts(:,1), trans_box_pts(:,2), trans_box_pts(:,3), siz, 'b')

hold off
zlim([0, 100])
ylim([-20, 15])


%% Read in the camera and image parameters
% camera.txt
% CAMERA_ID, MODEL, WIDTH, HEIGHT, PARAMS[]
% PARAMS = [focal length (in pixels), principal point x, principal point y, ?]
camera_file_path = 'COLMAP_Files\cameras.txt';

my_camera_file = fopen(camera_file_path);
my_camera_params = textscan(my_camera_file, '%f %s %f %f %f %f %f %f', 'CommentStyle','#');
fclose(my_camera_file);

% image.txt
% IMAGE_ID, QW, QX, QY, QZ, TX, TY, TZ, CAMERA_ID, NAME, POINTS2D[] as (X, Y, POINT3D_ID)
image_file_path = 'COLMAP_Files\images.txt';

my_image_file = fopen(image_file_path);
my_image_params = textscan(my_image_file, '%f %f %f %f %f %f %f %f %f %s', 'delimiter', '\t\t\t\t\t\t\t\t\t\t', 'CommentStyle','#');
fclose(my_image_file);
to_add = [];
for i = 1:numel(my_image_params)
    temp = my_image_params{i}(1:2:end);
    if i < numel(my_image_params)
        to_add = cat(2, to_add, my_image_params{i}(2:2:end));
    else
        last_col = my_image_params{i}(2:2:end);
    end
    my_image_params{i} = temp;
end

last_col = cellfun(@str2num, last_col, 'UniformOutput', false);
to_concat = {};
for i =1:length(to_add(:,1))
    temp = [to_add(i,:), last_col{i}];
    to_concat{i} = temp;
end
my_image_params{end + 1} = to_concat;


%% Place object in images
% read in image and place object in scene
images = [];
base_dir = 'Scene\';

for i = 1:numel(my_image_params{1,10})
    % read in image
    file_name = [base_dir my_image_params{1,10}{i}];
    image = imread(file_name);
    
    % get rotation
    rotation = quat2rotm([my_image_params{1,2}(i), my_image_params{1,3}(i), my_image_params{1,4}(i), my_image_params{1,5}(i)]);
    
    % get translation
    translation = [my_image_params{1,6}(i), my_image_params{1,7}(i), my_image_params{1,8}(i)];
    
    % principal points
    pp_x = my_camera_params{1,6};
    pp_y = my_camera_params{1,7};
    
    % focal length
    f = my_camera_params{1,5};
    
    [x, y] = CameraProjection3Dto2D(rotation, translation, pp_x, pp_y, f, trans_box_pts);

    depth_sorted = sortrows([x y], [2 1]);
    plane{1}= {[depth_sorted(1,:) depth_sorted(2,:) depth_sorted(6,:) depth_sorted(5,:)]};
    plane{2} = {[depth_sorted(5,:) depth_sorted(6,:) depth_sorted(8,:) depth_sorted(7,:)]};
    plane{3} = {[depth_sorted(3,:) depth_sorted(1,:) depth_sorted(5,:) depth_sorted(7,:)]};
    plane{4} = {[depth_sorted(4,:) depth_sorted(2,:) depth_sorted(6,:) depth_sorted(8,:)]};
    plane{5} = {[depth_sorted(1,:) depth_sorted(2,:) depth_sorted(4,:) depth_sorted(3,:)]};
    plane{6} = {[depth_sorted(3,:) depth_sorted(4,:) depth_sorted(8,:) depth_sorted(7,:)]};
    colors = {'yellow'; 'magenta'; 'cyan'; 'red'; 'green'; 'blue'};
    
    for k = 1: numel(plane)
        color = colors{k};
        image = insertShape(image, 'FilledPolygon', plane{k}, 'Color', color, 'Opacity', 0.7);
    end    
    figure, imshow(image)    
end


function [u, v, lambda] = CameraProjection3Dto2D(rotation, translation, pp_x, pp_y, fl, pts)
% set up
filmPlane2pixels = [1 0 pp_x; 0 1 pp_y; 0 0 1];
perspective_proj = [fl 0 0 0; 0 fl 0 0; 0 0 1 0];
rotation(end+1, :) = [0 0 0];
rotation(:, end+1) = [0; 0; 0; 1];
translation = [translation'; 1];
translation = [[diag([1 1 1]); 0 0 0], translation];
pts = [pts, ones(length(pts(:,1)),1)];

% multiply together
matrix = filmPlane2pixels*perspective_proj*rotation*translation;
pix_locs = matrix*pts';
pix_locs = pix_locs';
u = pix_locs(:,1);
v = pix_locs(:,2);
lambda = pix_locs(:,3);
u = round(u./lambda);
v = round(v./lambda);
end
