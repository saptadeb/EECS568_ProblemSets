% Saptadeep Debnath
% EECS 568 PS4

clc,clear,close all
import gtsam.*

%% Initialize iSAM
isam2 = ISAM2();

%% open file, read initial pose and graph
vertex_points = 1228;
file = '.\Data\input_INTEL_g2o.g2o';
[initial_values,graph_values] = task1a(file,vertex_points); 
sorted_graph    = sortrows(graph_values,2);

%% initialize graph and initial estimates
graph   = NonlinearFactorGraph;
initial = Values;

%% add prior noise
a=0.3;
sigma   = diag([a;a;a]);
noise   = sigma * rand(3,1);
priorNoise = noiseModel.Diagonal.Sigmas(noise);

%% initialise for key '0'
graph.add(PriorFactorPose2(0, Pose2(0, 0, 0), priorNoise)); % add directly to graph
initial.insert(initial_values(1,1),...
        Pose2(initial_values(1,2),initial_values(1,3),initial_values(1,4)));
isam2.update(graph, initial);
result = isam2.calculateEstimate();    

%% add vertex and related edges
c = 1;  %counter
for i=1:vertex_points-1
    
    % graph re-initialised to clear the clutter
    graph = NonlinearFactorGraph;
    initial.clear(); 
    
    % vertex added to initial estimates
    updated_x     = result.at(i-1).x + graph_values(i,3);
    updated_y     = result.at(i-1).y + graph_values(i,4);
    updated_theta = result.at(i-1).theta + graph_values(i,5);
    initial.insert(i,Pose2(updated_x,updated_y,updated_theta));
    
    % graph edges added to corresponding vertex
    while sorted_graph(c,2) == i
        picked_row  = sorted_graph(c,:);
        noisevec    = graph_values(i,6:11);
        noisemodel  = noiseModel.Gaussian.Covariance(info2Cov(noisevec));
        graph.add(BetweenFactorPose2(picked_row(1), picked_row(2),...
            Pose2(picked_row(3),picked_row(4),picked_row(5)), noisemodel));
        c = c+1;
        if c > size(sorted_graph,1)
            break
        end
    end
    
    % ISAM optimisation
    isam2.update(graph, initial);
    result = isam2.calculateEstimate();
end

% Pose2 extracted from 'result' for ease of plotting
result_poses = utilities.extractPose2(result);

%% plot the trajectories
figure();
plot(initial_values(:,2),initial_values(:,3),'r-');
hold on
plot(result_poses(:,1),result_poses(:,2),'b-');
legend('initial estimate','optimized trajectory')
title('Incremental Optimization - 2D INTEL dataset')    

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