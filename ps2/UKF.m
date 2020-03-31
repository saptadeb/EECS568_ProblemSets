classdef UKF < handle
    % UKF class for state estimation of a robot within a given map
    %
    %   Author: Saptadeep Debnath
    %   Date: 02/08/2020
    properties
        mu;             % Pose Mean
        Sigma;          % Pose Covariance
        gfun;           % Motion Model Function
        hfun;           % Measruement Model Function
        M;              % Motion model noise(dynamical and function of input)
        Q;              % Sensor Noise
        kappa_g;        
        mu_pred;
        Sigma_pred;
        n;
        X;
        w;
        Y;
        mean;
        X_x;
    end
    
    methods
        function obj = UKF(sys, init)
            % motion model
            obj.gfun = sys.gfun;
            % measurement model
            obj.hfun = sys.hfun;
            % motion noise covariance
            obj.M = sys.M;
            % measurement noise covariance
            obj.Q = sys.Q;
            obj.kappa_g = init.kappa_g;
            % initial mean and covariance
            obj.mu = init.mu;
            obj.Sigma = init.Sigma;
        end
        
        function prediction(obj, u)
            mu_aug = [obj.mu', zeros(1,3),zeros(1,2)]';
            Sigma_aug = blkdiag(obj.Sigma,obj.M(u),obj.Q);
            sigma_point(obj,mu_aug,Sigma_aug,obj.kappa_g);
            obj.mean = 0;
            for j = 1:2*obj.n+1
                obj.X_x(:,j) = obj.gfun(obj.X(1:3,j) , u + obj.X(4:6,j));
                obj.mean = obj.mean + obj.w(j) * obj.X_x(:,j);
            end
            obj.mu_pred = obj.mean;
            obj.mu_pred(3) = wrapToPi(obj.mu_pred(3));
            obj.Sigma_pred = zeros(size(obj.mu,1));
            for j = 1:2*obj.n+1
                diff = obj.X_x(:,j)-obj.mu_pred;
                diff(3) = wrapToPi(diff(3));
                temp = diff*obj.w(j)*diff';
                obj.Sigma_pred = obj.Sigma_pred+temp;
            end
        end
        
        function correction(obj,z,u)
            global FIELDINFO;
            landmark_x = FIELDINFO.MARKER_X_POS(z(3));
            landmark_y = FIELDINFO.MARKER_Y_POS(z(3));
            temp = 0;
            z_hat = 0;
            for j = 1:2*obj.n+1
                z_x(:,j) = obj.hfun(landmark_x, landmark_y,obj.X_x(:,j)) + obj.X(7:8,j);
                z_x(1,:) =  wrapToPi(z_x(1,:));
                temp = obj.w(j)*z_x(:,j);                
                z_hat = z_hat+temp;
                z_hat(1) = wrapToPi(z_hat(1));
            end
            S = zeros(size(z,1) - 1);
            Sigma_xz = zeros(3,2);
            for i=1:2*obj.n+1
                diff1 = (z_x(:,i) - z_hat);
                diff1(1) = wrapToPi(diff1(1));
                temp1 = obj.w(i) * (diff1)*(diff1)';
                S = S + temp1;
                diff2 = (obj.X_x(:,i) - obj.mu_pred);
                diff2(3) = wrapToPi(diff2(3));
                temp2 = obj.w(i) * (diff2)*(diff1)';
                Sigma_xz = Sigma_xz + temp2;
            end
            K = Sigma_xz*inv(S);
            z(1) = wrapToPi(z(1));
            v = z(1:2) - z_hat;
            v(1) = wrapToPi(v(1));
            obj.mu = obj.mu_pred + K*(v);
            obj.mu(3) = wrapToPi(obj.mu(3));
            obj.Sigma = obj.Sigma_pred - K*S*K';           
        end
        
        function sigma_point(obj, mean, cov, kappa)
            obj.n = numel(mean);
            L = sqrt(obj.n + kappa) * chol(cov,'lower');
            obj.Y = mean(:,ones(1, numel(mean)));
            obj.X = [mean,obj.Y + L, obj.Y - L];
            obj.w = zeros(2 * obj.n + 1,1);
            obj.w(1) = kappa / (obj.n + kappa);
            obj.w(2:end) = 0.5 / (obj.n + kappa);
        end   
    end
end