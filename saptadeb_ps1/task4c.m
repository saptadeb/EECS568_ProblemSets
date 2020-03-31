%% NA 568 - Problem Set 1
% Task 4C
% Saptadeep Debnath (saptadeb)

%%
clear,clc;

muX = 0;
varX = 1000;

Z1 = [10.6715 8.7925 10.7172 11.6302 10.4889 11.0347 10.7269 9.6966 ...
    10.2939 9.2127];
Z2 = [10.7107 9.0823 9.1449 9.3524 10.2602];

varZ1 = 1;
varZ2 = 0.64;

%% 1 (x conditioned on z1)
[muX_Z1,varX_Z1]=findPosterior(Z1,muX,varX,varZ1); % mean = 10.3255, var = 0.1000

%% 2 (x conditioned on z2)
[muX_Z2,varX_Z2]=findPosterior(Z2,muX,varX,varZ2); % mean = 9.7089, var = 0.1280

%% 3 (sensor choice)
precX_Z1 = 1/varX_Z1;    % precision = 10.0010
precX_Z2 = 1/varX_Z2;    % precision = 7.8125

%{
    It is found that when X is conditioned on Z1 it is more precise than
    when X is conditioned on Z2, even when Z2 as a sensor is more accurate
    with less variance. Since, sensor data for Z1 has more data points it
    increases the accuracy when X is conditioned on it, and finally
    computing a lower variance and hance a higher precision. 
%}

% inputing the same amount of data points for both sensors
[muX_Z1_2,varX_Z1_2]=findPosterior(Z1(randi(10,[1,5])),muX,varX,varZ1); 
% mean = 9.8047, var = 0.2000
[muX_Z2,varX_Z2]=findPosterior(Z2,muX,varX,varZ2);                      
% mean = 9.7089, var = 0.1280

precX_Z1_2 = 1/varX_Z1_2;    % precision = 5.0010
precX_Z2 = 1/varX_Z2;       % precision = 7.8125

%{
    As expected, when the same amount of randomly chosen data points are 
    given as the input to the system, sensor 2 (Z2) seems to provide a
    better precision than sensor 1 (Z1). In the end it depends on the
    application which kind of sensor to use. If in the application it is
    possible to collect more data points to recursively correct the sensor
    reading, sensor 1 proves to be the best option. But if higher accuracy
    is required over a short period of time, sensor 2 is preferred.
%}

%% function
function [mu, var] = findPosterior(data,prior_mu,prior_var,data_var)
    mu = prior_mu;       % prior mu
    var = prior_var;     % prior var
    for i=1:length(data)
        mu = (var*data(i)+mu*data_var)/(var+data_var); 
        var = (var*data_var)/(var+data_var);                 
    end
end
