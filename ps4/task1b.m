% Saptadeep Debnath
% EECS 568 PS4

clc,clear,close all
import gtsam.*
graph = NonlinearFactorGraph;

%% open file, read initial pose and graph
vertex_points = 1228;
file = '.\Data\input_INTEL_g2o.g2o';
[initial_values,graph_values] = task1a(file,vertex_points); 

%% add prior noise
a=0.3;
sigma       = diag([a;a;a]);
noise       = sigma * rand(3,1);
priorNoise  = noiseModel.Diagonal.Sigmas(noise);

%% add edge values
for i=1:1483
    noisevec = graph_values(i,6:11);
    graph.add(BetweenFactorPose2(graph_values(i,1), graph_values(i,2),...
        Pose2(graph_values(i,3), graph_values(i,4), graph_values(i,5)), ...
        noiseModel.Gaussian.Covariance(info2Cov(noisevec))));
end   
graph.add(PriorFactorPose2(1, Pose2(0, 0, 0), priorNoise)); % add directly to graph

%% initialize initial estimates
initial = Values;
for i=1:1228
    initial.insert(initial_values(i,1),...
        Pose2(initial_values(i,2),initial_values(i,3),initial_values(i,4)));
end

%% optimize using GaussNewtonOptimizer
optimizer       = GaussNewtonOptimizer(graph, initial);
result          = optimizer.optimizeSafely();
result_final    = utilities.extractPose2(result);

%% plot covariance ellipses
figure();
plot(initial_values(:,2),initial_values(:,3),'r-');
hold on
plot(result_final(:,1),result_final(:,2),'b-');
legend('initial estimate','optimized trajectory')
title('Batch Optimization - 2D INTEL dataset')

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