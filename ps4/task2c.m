% Saptadeep Debnath
% EECS 568 PS4

clc,clear,close all
import gtsam.*
graph = NonlinearFactorGraph;

%% Initialize iSAM
isam2 = ISAM2();

%% open file, read initial pose and graph
vertex_points = 1661;
file = '.\Data\parking-garage.g2o';
[initial_values,graph_values] = task2a(file,vertex_points); 
sorted_graph = sortrows(graph_values,2);

%% initialize graph and initial estimates
graph = NonlinearFactorGraph;
initial = Values;

%% add prior noise
a = 1e-5;
b = 1e-5;
model = noiseModel.Diagonal.Sigmas([a;a;a;b;b;b]);

%% initialise for key '0'
graph.add(PriorFactorPose3(0, Pose3(), model)); % add directly to graph
initial.insert(initial_values(1,1),Pose3());
isam2.update(graph, initial);
result = isam2.calculateEstimate();    

%% add vertex and related edges
odometry_mat = graph_values(graph_values(:,1)+1 == graph_values(:,2),:);
c = 1; % counter
for i=1:vertex_points-1
    
    % graph re-initialised to clear the clutter
    graph = NonlinearFactorGraph;
    initial.clear();   
    
    % vertex added to initial estimates
    updated_x = result.at(i-1).x + odometry_mat(i,3);
    updated_y = result.at(i-1).y + odometry_mat(i,4);
    updated_z = result.at(i-1).z + odometry_mat(i,5);
    noisevec = odometry_mat(i,10:30);
    r_new = Rot3(quat2rotm([odometry_mat(i,9),odometry_mat(i,6),odometry_mat(i,7),odometry_mat(i,8)]));
    updated_rot = result.at(i-1).rotation.compose(r_new);
    t_new = Point3(updated_x,updated_y,updated_z);
    initial.insert(i,Pose3(updated_rot,t_new));
    
    % graph edges added to corresponding vertex
    while sorted_graph(c,2) == i
        clear r t noisevec
        picked_row = sorted_graph(c,:);
        noisevec = picked_row(10:30);
        noisemodel = noiseModel.Gaussian.Covariance(info2Cov(noisevec));
        r = Rot3(quat2rotm([picked_row(9),picked_row(6),picked_row(7),picked_row(8)]));
        t = Point3(picked_row(3),picked_row(4),picked_row(5));
        graph.add(BetweenFactorPose3(picked_row(1), picked_row(2),...
            Pose3(r,t), noisemodel));
        c = c+1;
        if c > size(sorted_graph,1)
            break
        end
    end
    
    % ISAM optimisation   
    isam2.update(graph, initial);
    result = isam2.calculateEstimate();
end 
result_final    = utilities.extractPose3(result);

%% for plotting ease
initial = Values;
for i=1:1661
    r = Rot3(quat2rotm([initial_values(i,8),initial_values(i,5),initial_values(i,6),initial_values(i,7)]));
    t = Point3(initial_values(i,2),initial_values(i,3),initial_values(i,4));
    initial.insert(initial_values(i,1),Pose3(r,t));
end    

%% plot covariance ellipses
figure();
plot3(initial_values(:,2),initial_values(:,3),initial_values(:,4),'r-');
hold on
plot3(result_final(:,10),result_final(:,11),result_final(:,12),'b-');
legend('initial estimate','optimized trajectory')
title('Incremental Optimization - 3D Parking Garage dataset')
axis equal
view(-90,90) % for top down
% view(-45,33)
% colormap hot

%% function to convert given information matrix to covariance matrix
function cov_matrix = info2Cov(info_vector)

% 2D case- there should be 6 elements to form upper triangle of 3X3 matrix
if numel(info_vector) == 6
    info_matrix =  [info_vector(1) info_vector(2) info_vector(3);
                    info_vector(2) info_vector(4) info_vector(5);
                    info_vector(3) info_vector(5) info_vector(6)];

% 3D case- there should be 21 elements to form upper triangle of 6X6 matrix
elseif numel(info_vector) == 21
    info_matrix = zeros(6,6);
    info_matrix(1,1:6) = info_vector(1:6);
    info_matrix(2,2:6) = info_vector(7:11);
    info_matrix(3,3:6) = info_vector(12:15);
    info_matrix(4,4:6) = info_vector(16:18);
    info_matrix(5,5:6) = info_vector(19:20);
    info_matrix(6,6) = info_vector(21);
    info_matrix = info_matrix + info_matrix' - diag(diag(info_matrix));
end

cov_matrix = inv(info_matrix);
cov_matrix = chol(cov_matrix,'lower');
end