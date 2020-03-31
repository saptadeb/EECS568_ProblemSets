% Saptadeep Debnath
% EECS 568 PS4

clc,clear,close all
import gtsam.*
graph = NonlinearFactorGraph;

%% open file, read initial pose and graph
vertex_points = 1661;
file = '.\Data\parking-garage.g2o';
[initial_values,graph_values] = task2a(file,vertex_points); 

%% add prior noise
a = 1e-5;
b = 1e-5;
model = noiseModel.Diagonal.Sigmas([a;a;a;b;b;b]);
graph.add(PriorFactorPose3(1, Pose3(), model))

%% add edge values
for i=1:6275
    noisevec = graph_values(i,10:30);
    r = Rot3(quat2rotm([graph_values(i,9),graph_values(i,6),graph_values(i,7),graph_values(i,8)]));
    t = Point3(graph_values(i,3),graph_values(i,4),graph_values(i,5));
    graph.add(BetweenFactorPose3(graph_values(i,1),graph_values(i,2),Pose3(r,t),...
        noiseModel.Gaussian.Covariance(info2Cov(noisevec))));
end   

%% initialize initial estimates
initial = Values;
for i=1:1661
    r = Rot3([quat2rotm([initial_values(i,8),initial_values(i,5),initial_values(i,6),initial_values(i,7)])]);
    t = Point3(initial_values(i,2),initial_values(i,3),initial_values(i,4));
    initial.insert(initial_values(i,1),Pose3(r,t));
end    

%% optimize using GaussNewtonOptimizer
optimizer       = GaussNewtonOptimizer(graph, initial);
result          = optimizer.optimizeSafely();
result_final    = utilities.extractPose3(result);

%% plot covariance ellipses
figure();
plot3(initial_values(:,2),initial_values(:,3),initial_values(:,4),'r-');
hold on
plot3(result_final(:,10),result_final(:,11),result_final(:,12),'b-');
legend('initial estimate','optimized trajectory')
title('Batch Optimization - 3D Parking Garage dataset')
axis equal
% view(-90,90) % for top down
view(-45,33)
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