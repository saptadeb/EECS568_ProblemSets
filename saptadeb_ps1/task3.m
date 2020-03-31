%% NA 568 - Problem Set 1
% Task 3
% Saptadeep Debnath (saptadeb)

%% Part A (Generating point cloud)
clear, clc;

% initialising mean and standard deviation for range and bearing 
muR = 10.0; %in m
sdR = 0.5;  %in m
muB = 0;    %in rad
sdB = 0.25; %in rad

% (i) generating point cloud in Sensor Frame
for i=1:10000
    r(i) = muR+sdR*randn;
    b(i) = muB+sdB*randn;
end

% (ii) generating point could in Cartesian Coordinate Frame
for i=1:10000
    x(i) = r(i)*cos(b(i));
    y(i) = r(i)*sin(b(i));
end

figure(1);
clf;

subplot(121);
scatter(r,b,'.')
axis square
h = xlabel('range, m');
set(h,'interpreter','latex','fontsize',14);
h = ylabel('bearing, rad');
set(h,'interpreter','latex','fontsize',14);
h = title('Sensor Frame');
set(h,'interpreter','latex','fontsize',14);

subplot(122);
scatter(x,y,'.')
axis square
h = xlabel('X, m');
set(h,'interpreter','latex','fontsize',14);
h = ylabel('Y, m');
set(h,'interpreter','latex','fontsize',14);
h = title('Cartesian Coordinate Frame');
set(h,'interpreter','latex','fontsize',14);

%% Part B (Linearized covariance in cartesian coordinates)
clc;

% equilibrium point 
r_eq = muR;
b_eq = muB;
muP = [r_eq;b_eq];

syms a theta
f = a*cos(theta); % function of 'x'
g = a*sin(theta); % function of 'y'

j = jacobian([f;g], [a;theta]);       % calculating jacobian 
J = subs(j, [a theta], [r_eq b_eq]);  % absolute value of the jacobian 

B = subs([f;g], [a theta], [r_eq b_eq]) - J*[r_eq;b_eq]; % affine translation

muC = J*muP+B;
covP = diag([sdR^2,sdB^2]);         % covariance in polar coordinate system
covC = double(J*covP*J.');           % covariance in cartesian coordinates

%% Part C (Drawing sigma contours) 
clc;

actual_mu_xy = [mean(x);mean(y)];
actual_cov_xy = cov(x,y);

actual_mu_rb = [mean(r);mean(b)];
actual_cov_rb = cov(r,b);

figure(2);
clf;

h1 = scatter(r,b,'.');
hold on
axis square
h2 = draw_ellipse(actual_mu_rb, actual_cov_rb, 1);
set(h2,'Color','blue')
h = draw_ellipse(actual_mu_rb, actual_cov_rb, 4);
set(h,'Color','blue')
h = draw_ellipse(actual_mu_rb, actual_cov_rb, 9);
set(h,'Color','blue')
h3 = draw_ellipse(muP, covP, 1);
set(h3,'Color','red')
h = draw_ellipse(muP, covP, 4);
set(h,'Color','red')
h = draw_ellipse(muP, covP, 9);
set(h,'Color','red')
h = legend([h1 h3 h2],'data points','analytical contour','sample contour');
set(h,'interpreter','latex','fontsize',14);
h = xlabel('range, m');
set(h,'interpreter','latex','fontsize',14);
h = ylabel('bearing, rad');
set(h,'interpreter','latex','fontsize',14);
h = title('Sensor Frame');
set(h,'interpreter','latex','fontsize',14);

figure(3)
h1 = scatter(x,y,'.');
hold on
axis equal
h2 = draw_ellipse(actual_mu_xy, actual_cov_xy, 1);
set(h2,'Color','blue')
h = draw_ellipse(actual_mu_xy, actual_cov_xy, 4);
set(h,'Color','blue')
h = draw_ellipse(actual_mu_xy, actual_cov_xy, 9);
set(h,'Color','blue')
h3 = draw_ellipse(muC, covC, 1);
set(h3,'Color','red')
h = draw_ellipse(muC, covC, 4);
set(h,'Color','red')
h = draw_ellipse(muC, covC, 9);
set(h,'Color','red')
h = legend([h1 h3 h2],'data points','analytical contour','sample contour');
set(h,'interpreter','latex','fontsize',14);
h = xlabel('X, m');
set(h,'interpreter','latex','fontsize',14);
h = ylabel('Y, m');
set(h,'interpreter','latex','fontsize',14);
h = title('Cartesian Coordinate Frame');
set(h,'interpreter','latex','fontsize',14);

%% Part D (Samples lying within sigma-contours) 
clc;

% for Cartesian system
figure(3);
clf;
clear in_percent md
k = [1 4 9];
for i = 1:length(k)
    [xv,yv] = calculateEllipseXY(muC, covC, k(i));
    xv = double(xv);
    yv = double(yv);

    [in,on] = inpolygon(x,y,xv,yv);
    p_in = numel(x(in));
    p_on = numel(x(on));
    p_out = numel(x(~in));

    in_percent(i) = (p_in)*100/(p_in+p_on+p_out);
    md(i) = mean(mahal([xv yv],[x.' y.']));
    
    subplot(2,3,3+i)
    plot(xv,yv)                     % polygon
    axis square

    hold on
    plot(x(in),y(in),'r.')          % points inside
    plot(x(~in),y(~in),'b.')        % points outside
end
set(findall(gcf,'type','line'),'linewidth',2)

% for Polar system
clear in_percent md
k = [1 4 9];

for i = 1:length(k)
    [xv,yv] = calculateEllipseXY(muP, covP, k(i));
    xv = double(xv);
    yv = double(yv);

    [in,on] = inpolygon(r,b,xv,yv);
    p_in = numel(r(in));
    p_on = numel(r(on));
    p_out = numel(r(~in));

    in_percent(i) = (p_in)*100/(p_in+p_on+p_out);
    md(i) = mean(mahal([xv yv],[r.' b.']));
    
    subplot(2,3,i)
    plot(xv,yv) % polygon
    axis square

    hold on
    plot(r(in),b(in),'r.')          % points inside
    plot(r(~in),b(~in),'b.')        % points outside
end
set(findall(gcf,'type','line'),'linewidth',2)

%% Part E (Varying the noise parameters to match theoretical values)
clc;

% manually changing the standard deviations 
muR = 10.0; %in m
sdR = 0.01; %in m
muB = 0; %in rad
sdB = 0.01; %in rad

for i=1:10000
    r(i) = muR+sdR*randn;
    b(i) = muB+sdB*randn;
end

for i=1:10000
    x(i) = r(i)*cos(b(i));
    y(i) = r(i)*sin(b(i));
end

% equilibrium point 
r_eq = muR;
b_eq = muB;
muP = [r_eq;b_eq];

syms a theta
f = a*cos(theta); % function of 'x'
g = a*sin(theta); % function of 'y'

j = jacobian([f;g], [a;theta]);       % calculating jacobian 
J = subs(j, [a theta], [r_eq b_eq]);  % absolute value of the jacobian 

B = subs([f;g], [a theta], [r_eq b_eq]) - J*[r_eq;b_eq]; % affine translation

muC = J*muP+B;
covP = diag([sdR^2,sdB^2]);         % covariance in polar coordinate system
covC = double(J*covP*J.');          % covariance in cartesian coordinates

clear in_percent md
k = [1 4 9];

for i = 1:length(k)
    [xv,yv] = calculateEllipseXY(muC, covC, k(i));
    xv = double(xv);
    yv = double(yv);

    [in,on] = inpolygon(x,y,xv,yv);
    p_in = numel(x(in));
    p_on = numel(x(on));
    p_out = numel(x(~in));

    in_percent(i) = (p_in)*100/(p_in+p_on+p_out);
    md(i) = mean(mahal([xv yv],[x.' y.']));
end
set(findall(gcf,'type','line'),'linewidth',2)

%% Part F (measurements are not independent)
clear, clc;

muR = 10.0; %in m
sdR = 0.5; %in m
muB = 0; %in rad
sdB = 0.25; %in rad
muP = [muR;muB];

rho_rb = [0.1 0.5 0.9];
clear covC
for j = 1:length(rho_rb)
    covP(:,:,j) = [sdR^2 rho_rb(j)*sdB*sdR ; rho_rb(j)*sdB*sdR sdB^2];
    L = chol(covP(:,:,j)).';
    for i=1:10000
        p(i,:,j) = L*randn(2,1)+muP;
    end
    for i=1:10000
        c(i,:,j) = [p(i,1,j)*cos(p(i,2,j));p(i,1,j)*sin(p(i,2,j))];
    end
end

syms a theta
f = a*cos(theta); % function of 'x'
g = a*sin(theta); % function of 'y'

j = jacobian([f;g], [a;theta]);       % calculating jacobian 
J = subs(j, [a theta], [muR muB]);    % absolute value of the jacobian 

B = subs([f;g], [a theta], [muR muB]) - J*[muR;muB]; % affine translation

muC = J*muP+B;
for j = 1:3
    covC(:,:,j) = double(J*covP(:,:,j)*J.');% covariance in cartesian coordinates
end

actual_mu_xy = [mean(c)];
actual_mu_rb = [mean(p)];

for i=1:3
    actual_cov_xy(:,:,i) = cov(c(:,1,i),c(:,2,i));
    actual_cov_rb(:,:,i) = cov(p(:,1,i),p(:,2,i));
end

for i = 1:3
    figure();
    h1 = scatter(p(:,1,i),p(:,2,i),'.');
    hold on
    axis square
    h2 = draw_ellipse(actual_mu_rb(:,:,i), actual_cov_rb(:,:,i), 1);
    set(h2,'Color','blue')
    h = draw_ellipse(actual_mu_rb(:,:,i), actual_cov_rb(:,:,i), 4);
    set(h,'Color','blue')
    h = draw_ellipse(actual_mu_rb(:,:,i), actual_cov_rb(:,:,i), 9);
    set(h,'Color','blue')
    h3 = draw_ellipse(muP, covP(:,:,i), 1);
    set(h3,'Color','red')
    h = draw_ellipse(muP, covP(:,:,i), 4);
    set(h,'Color','red')
    h = draw_ellipse(muP, covP(:,:,i), 9);
    set(h,'Color','red')
    h = legend([h1 h3 h2],'data points','analytical contour','sample contour');
    set(h,'interpreter','latex','fontsize',14);
    h = xlabel('range, m');
    set(h,'interpreter','latex','fontsize',14);
    h = ylabel('bearing, rad');
    set(h,'interpreter','latex','fontsize',14);
    h = title('Sensor Frame');
    set(h,'interpreter','latex','fontsize',14);
    hold off;
end

for i = 1:3
    figure();
    h1 = scatter(c(:,1,i),c(:,2,i),'.');
    hold on
    axis equal
    h2 = draw_ellipse(actual_mu_xy(:,:,i), actual_cov_xy(:,:,i), 1);
    set(h2,'Color','blue')
    h = draw_ellipse(actual_mu_xy(:,:,i), actual_cov_xy(:,:,i), 4);
    set(h,'Color','blue')
    h = draw_ellipse(actual_mu_xy(:,:,i), actual_cov_xy(:,:,i), 9);
    set(h,'Color','blue')
    h3 = draw_ellipse(muC, covC(:,:,i), 1);
    set(h3,'Color','red')
    h = draw_ellipse(muC, covC(:,:,i), 4);
    set(h,'Color','red')
    h = draw_ellipse(muC, covC(:,:,i), 9);
    set(h,'Color','red')
    h = legend([h1 h3 h2],'data points','analytical contour','sample contour');
    set(h,'interpreter','latex','fontsize',14);
    h = xlabel('X, m');
    set(h,'interpreter','latex','fontsize',14);
    h = ylabel('Y, m');
    set(h,'interpreter','latex','fontsize',14);
    h = title('Cartesian Coordinate Frame');
    set(h,'interpreter','latex','fontsize',14);
    hold off;
end
