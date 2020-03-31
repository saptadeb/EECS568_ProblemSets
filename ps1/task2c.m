%% NA 568 - Problem Set 1
% Task 2C
% Saptadeep Debnath (saptadeb)

%% Initialising variables
clear, clc;

p = linspace(0,9,10);
L = [0 3 6];
P_L = [];       % position being a landmark

for i=1:length(p)
    if ismember(p(i),L)
        P_L(i) = 1;
    else
        P_L(i) = 0;
    end
end

% initial belief
for i=1:length(p)
    bel(i) = 1/length(p); 
end

% sensor model
L_L = 0.8;       % probability of sensing landmark given landmark
nL_L = 1-L_L;    % probability of sensing no landmark given landmark
L_nL = 0.4;      % probability of sensing landmark given no landmark
nL_nL = 1-L_nL;  % probability of sensing no landmark given no landmark

%% first step - robot detects a landmark (i+3)
P_z = [];
for i=1:length(p)
    if P_L(i)
        P_z(i) = L_L;
    else
        P_z(i) = L_nL;
    end
end

for i=1:length(p)
    j = mod((i+2),10)+1;
    bel_dash(i) = P_z(j)*bel(i);
end
eta = sum(bel_dash); % normalizer

for i=1:length(p)
    bel(i) = inv(eta)*bel_dash(i);
end
figure(1)
clf
bar(p,bel)
title('Belief after the First sequence')
xlabel('Position')
ylabel('Probability')

%% second step - robot detects a landmark (i+3+3)
clear P_z;
for i=1:length(p)
    if P_L(i)
        P_z(i) = L_L;
    else
        P_z(i) = L_nL;
    end
end

for i=1:length(p)
    j = mod((i+5),10)+1;
    bel_dash(i) = P_z(j)*bel(i);
end
eta = sum(bel_dash); % normalizer

for i=1:length(p)
    bel(i) = inv(eta)*bel_dash(i);
end
figure(2)
clf
bar(p,bel)
title('Belief after the Second sequence')
xlabel('Position')
ylabel('Probability')

%% third step - robot doesn't detect a landmark (i+3+3+4)
clear P_z;
for i=1:length(p)
    if P_L(i)
        P_z(i) = nL_L;
    else
        P_z(i) = nL_nL;
    end
end

for i=1:length(p)
    j = mod((i+9),10)+1;
    bel_dash(i) = P_z(j)*bel(i);
end
eta = sum(bel_dash); % normalizer

for i=1:length(p)
    bel(i) = inv(eta)*bel_dash(i);
end
figure(3)
clf
bar(p,bel)
title('Belief after the Final sequence')
xlabel('Position')
ylabel('Probability')

%% Result - 

%{
    It is clearly seen from the graph that the robot's final position seems
    most likely to be position 7 (P=30%), thus can be infered that the
    robot started from position 0.
%}

