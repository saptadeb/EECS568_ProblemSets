classdef InEKF < handle   
    properties
        mu;                 % Pose Mean
        Sigma;              % Pose Sigma
        gfun;               % Motion model function
        mu_pred;             % Mean after prediction step
        Sigma_pred;          % Sigma after prediction step
        mu_cart;
        sigma_cart;
        q;
        v;
        twst;
        Adj;
    end
    
    methods
        function obj = InEKF(sys, init)
            obj.gfun = sys.gfun;
            obj.mu = init.mu;
            obj.Sigma = init.Sigma;
            obj.q = 1*[2.5e-4 0 0; 0 2.5e-4 0; 0 0 2.5e-4];
            obj.v = 1000*diag([1,1]);
        end
        
        function w = wedge(obj, X)
            w = [0 -X(3) X(1); X(3) 0 X(2); 0 0 0];
        end
        
        function tw = twist(obj, u)
            x = obj.mu(1,3);
            y = obj.mu(2,3);
            theta = atan2(obj.mu(2,1), obj.mu(1,1));
            obj.mu_cart = [x; y; theta];
            X_p = obj.gfun(obj.mu_cart, u);
            X_pred = obj.posemat(X_p);
            tw = logm(inv(obj.mu) * X_pred);
        end
        
        function prediction(obj, u)
            obj.Adj = [obj.mu(1:2,1:2), [obj.mu(2,3); -obj.mu(1,3)]; 0 0 1];
            obj.twst = twist(obj, u);
            propagation(obj,u);
        end
        
        function propagation(obj, u)
            obj.mu_pred = obj.mu * expm(obj.twst);
            obj.Sigma_pred = obj.Sigma + obj.Adj * obj.q * obj.Adj';
        end
        
        function correction(obj, Y, Y2, landmark_ids)
            global FIELDINFO;        
            landmark_x = FIELDINFO.MARKER_X_POS(landmark_ids(1));
            landmark_y = FIELDINFO.MARKER_Y_POS(landmark_ids(1));       
            landmark_x2 = FIELDINFO.MARKER_X_POS(landmark_ids(2));
            landmark_y2 = FIELDINFO.MARKER_Y_POS(landmark_ids(2));
            
            H1 = [-1 0 landmark_y; 0 -1 -landmark_x];
            H2 = [-1 0 landmark_y2; 0 -1 -landmark_x2];
            H = [H1;H2];
            N = obj.mu_pred * blkdiag(obj.v,0) * obj.mu_pred';
            N = N(1:2,1:2);
            N_blk = blkdiag(N,N);
            S = H * obj.Sigma_pred * H' + N_blk;
            L = obj.Sigma_pred * H' * inv(S);
            b1 = [landmark_x; landmark_y; 1];
            b2 = [landmark_x2; landmark_y2; 1];
            diff1 = obj.mu_pred*Y' - b1;
            diff1 = diff1(1:2);
            diff2 = obj.mu_pred*Y2' - b2;
            diff2 = diff2(1:2);
            diff = [diff1;diff2];
            del = obj.wedge((L * diff));
            obj.mu = expm(del) * obj.mu_pred;
            obj.Sigma = (eye(3) - L * H) * obj.Sigma_pred *...
                (eye(3) - L * H)' + L * N_blk * L';
        end
        
        function H = posemat(obj,state)
            x = state(1);
            y = state(2);
            h = state(3);
            % construct a SE(2) matrix element
            H = [...
                cos(h) -sin(h) x;
                sin(h)  cos(h) y;
                     0       0 1];
        end
    end
end
