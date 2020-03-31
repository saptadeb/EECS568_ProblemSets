classdef PF < handle
    % PF class for state estimation of a robot within a given map
    %
    %   Author: Saptadeep Debnath
    %   Date: 02/08/2020    
    properties
        gfun;               %Motion model function
        hfun;               %Measurement Model Function
        Q;                  %Sensor Noise
        M;                  %Motion Model Noise (dynamical and function of input)
        n;                  %Number of Particles
        particles;          %Pose of particle
        particle_weight;    %Particle Weight
        mu;
        Sigma;
    end
    
    methods
        function obj = PF(sys, init)
            % motion model
            obj.gfun = sys.gfun;
            % measurement model
            obj.hfun = sys.hfun;
            % motion noise covariance
            obj.M = sys.M;
            % measurement noise covariance
            obj.Q = sys.Q;
            % PF parameters
            obj.mu = init.mu;
            obj.Sigma = init.Sigma;
            obj.n = init.n;
            obj.particles = init.particles;
            obj.particle_weight = init.particle_weight;
        end
        
        function prediction(obj, u)
            for i = 1:obj.n
                obj.particles(:,i) = obj.gfun(obj.particles(:,i), ...
                    chol(obj.M(u),"lower") * randn(3,1) + u );
            end
        end

        function correction(obj, z)
            global FIELDINFO;
            landmark_x = FIELDINFO.MARKER_X_POS(z(3));
            landmark_y = FIELDINFO.MARKER_Y_POS(z(3));
            weight_prob = zeros(obj.n,1);
            for i= 1:obj.n    
                z_k = obj.hfun(landmark_x,landmark_y,obj.particles(:,i));
                v = z(1:2)-z_k;
                v(1) = wrapToPi(v(1));
                weight_prob(i) = mvnpdf(v,[0;0], obj.Q);
            end
            obj.particle_weight = obj.particle_weight .* weight_prob;
            obj.particle_weight = obj.particle_weight./sum(obj.particle_weight);
            neff = 1/sum(obj.particle_weight.^2);
            if neff<obj.n/1.5
                resample(obj);
            end
            meanAndVariance(obj);
        end 
        
        function resample(obj)
            newSamples = zeros(size(obj.particles));
            newWeight = zeros(size(obj.particle_weight));
            W = cumsum(obj.particle_weight);
            r = rand/obj.n;
            count = 1;
            for j = 1:obj.n
                u = r+(j-1)/obj.n;
                while u > W(count)
                    count = count+1;
                end
                newSamples(:,j) = obj.particles(:,count);
                newWeight(j) = 1/obj.n;
            end
            obj.particles = newSamples;
            obj.particle_weight = newWeight;
        end
        
        function meanAndVariance(obj)
            obj.mu = mean(obj.particles, 2); 
            % orientation is a bit more tricky.
            sinSum = 0;
            cosSum = 0;
            for s = 1:obj.n
                cosSum = cosSum + cos(obj.particles(3,s));
                sinSum = sinSum + sin(obj.particles(3,s));
            end
            obj.mu(3) = atan2(sinSum, cosSum);     
            % Compute covariance.
            zeroMean = obj.particles - repmat(obj.mu, 1, obj.n);
            for s = 1:obj.n
                zeroMean(3,s) = wrapTo2Pi(zeroMean(3,s));
            end
            obj.Sigma = zeroMean * zeroMean' / obj.n;
        end
    end
end