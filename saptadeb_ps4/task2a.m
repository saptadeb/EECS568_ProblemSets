function [initial_values,graph_values] = task2a(file,vertex_points) 
% function to read a g2o file for 3D dataset and return the initial
% estimates and graph constraints+odometry

%% open file, read initial pose and graph
    fid = fopen(file);
    C_vertex = textscan(fid,'%*s %f %f %f %f %f %f %f %f',vertex_points);
    C_edge = textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
    fclose(fid);

    %% convert cell to matrix
    initial_values = cell2mat(C_vertex);
    graph_values = cell2mat(C_edge);

end