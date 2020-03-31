classdef semantic_grid_map_continuous_S_CSM < handle
    properties
        % map dimensions
        range_x = [-15, 20];
        range_y = [-25, 10];
        % sensor parameters
        z_max = 30;                 % max range in meters
        n_beams = 133;              % number of beams
        % grid map paremeters
        grid_size = 0.27;
        alpha = 2 * 0.135;          % 2 * grid_size
        beta = 2 * pi/133;          % 2 * pi/n_beams
        nn = 16;                    % number of nearest neighbor search
        map;                        % map!
        pose;                       % pose data
        scan;                       % laser scan data
        m_i = [];                   % cell i
        % semantic
        num_classes = 6;
        
        l = 0.2;
        sigma = 0.1;
    end
    
    methods
        function obj = semantic_grid_map_continuous_S_CSM(pose, scan)
            % class constructor
            % construct map points, i.e., grid centroids.
            x = obj.range_x(1):obj.grid_size:obj.range_x(2);
            y = obj.range_y(1):obj.grid_size:obj.range_y(2);
            [X,Y] = meshgrid(x,y);
            t = [X(:), Y(:)];
            % a simple KDtree data structure for map coordinates.
            obj.map.occMap = KDTreeSearcher(t);
            obj.map.size = size(t,1);
            
            % -----------------------------------------------
            % To Do: 
            % map parameter initialization such as map.alpha
            obj.map.alpha = 1e-6 * ones(size(t,1),obj.num_classes+1);     
            
            for i=1:obj.num_classes+1
                alpha_sum = 0;
                for j = 1:obj.num_classes+1
                    alpha_sum = alpha_sum+obj.map.alpha(:,j);
                end
                obj.map.mean(:,i) = obj.map.alpha(:,i)./alpha_sum;

            end    
            obj.map.variance = (obj.map.mean(:,1).*(1-obj.map.mean(:,1)))./...
                (alpha_sum+1);    
            % -----------------------------------------------
            
            % set robot pose and laser scan data
            obj.pose = pose;
            obj.pose.mdl = KDTreeSearcher([pose.x, pose.y]);
            obj.scan = scan;
        end
        
        function build_ogm(obj)
            % build occupancy grid map using the binary Bayes filter.
            % we first loop over all map cells, then for each cell, we find
            % N nearest neighbor poses to build the map. Note that this is
            % more efficient than looping over all poses and all map cells
            % for each pose which should be the case in online
            % (incremental) data processing.
            for i = 1:obj.map.size
                m = obj.map.occMap.X(i,:);
                idxs = knnsearch(obj.pose.mdl, m, 'K', obj.nn);
                if ~isempty(idxs)
                    for k = idxs
                        % pose k
                        pose_k = [obj.pose.x(k),obj.pose.y(k), obj.pose.h(k)];
                        if obj.is_in_perceptual_field(m, pose_k)
                            % laser scan at kth state; convert from
                            % cartesian to polar coordinates
                            [bearing, range] = cart2pol(obj.scan{k}(1,:), obj.scan{k}(2,:));
                            z = [range' bearing' obj.scan{k}(3,:)'];
                            
                            % -----------------------------------------------
                            % To Do: 
                            % update the semantic sensor model in cell i
                            obj.continuous_S_CSM(z,i);                            
                            % -----------------------------------------------
                        end
                    end
                end
                
                % -----------------------------------------------
                % To Do: 
                % update mean and variance for each cell i
                alpha_sum = sum(obj.map.alpha(i,:));
                obj.map.mean(i,:) = obj.map.alpha(i,:)/alpha_sum;
                maxmean = max(obj.map.mean(i,:));
                obj.map.variance(i) = (maxmean*(1-maxmean))/(alpha_sum+1);  
                % -----------------------------------------------
                
            end
        end
        
        function inside = is_in_perceptual_field(obj, m, p)
            % check if the map cell m is within the perception field of the
            % robot located at pose p.
            inside = false;
            d = m - p(1:2);
            obj.m_i.range = sqrt(sum(d.^2));
            obj.m_i.phi = wrapToPi(atan2(d(2),d(1)) - p(3));
            % check if the range is within the feasible interval
            if (0 < obj.m_i.range) && (obj.m_i.range < obj.z_max)
                % here sensor covers -pi to pi!
                if (-pi < obj.m_i.phi) && (obj.m_i.phi < pi)
                    inside = true;
                end
            end
        end
        
        function continuous_S_CSM(obj, z, i)
           % -----------------------------------------------
           % To Do: 
           % implement the continuous semantic counting sensor model
           bearing_diff = abs(wrapToPi(z(:,2) - obj.m_i.phi));            
           [bearing_min, k] = min(bearing_diff);            
           if obj.m_i.range > min(obj.z_max, z(k,1) + obj.alpha/2) || bearing_min > obj.beta/2                
               %do nothing in this case            
           elseif z(k,1) < obj.z_max && abs(obj.m_i.range - z(k,1)) < obj.alpha/2                 
               [x,y] = pol2cart(z(k,2),z(k,1));
               [x_star_x, x_star_y] = pol2cart(obj.m_i.phi, obj.m_i.range);
               d = sqrt((x_star_x-x)^2+(x_star_y-y)^2);
               if d < obj.l
                K = obj.kernel(d);
                obj.map.alpha(i,z(k,3)) = obj.map.alpha(i,z(k,3))+K;          
               end
           elseif obj.m_i.range <  z(k,1) && z(k,1) < obj.z_max                
               K = obj.sampling(z(k,:));
               obj.map.alpha(i,obj.num_classes+1) = obj.map.alpha(i,obj.num_classes+1)+K;
           end
           % -----------------------------------------------
        end
                function k = kernel(obj,d)
            dbyl = d/obj.l;
            k = obj.sigma*((((2+cos(2*pi*dbyl))*(1-dbyl))/3)+((sin(2*pi*dbyl))/(2*pi)));
        end

        function K = sampling(obj,z)
            K = 0;
            r(:,2) = [obj.m_i.range + obj.l ; obj.m_i.range - obj.l ;
                obj.m_i.range + (obj.l/3) ; obj.m_i.range - (obj.l/3)];
            r(:,1) = z(2);
            [x,y] = pol2cart(r(:,1),r(:,2));
            [x_star_x, x_star_y] = pol2cart(obj.m_i.phi, obj.m_i.range);
            d = sqrt((x_star_x - x).^2+(x_star_y - y).^2);
            for i=1:length(d)
                if d(i) < obj.l
                    k = obj.kernel(d(i));
                    K = K+k;  
                end
            end
        end
    end
end