% https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method
% https://chejunkie.com/knowledge-base/simplex-nelder-mead-optimization-amoeba-search/#Inside_Contraction
classdef NelderMead < handle
    properties
        J
        n
        x
        start
        printInfo
        rho = 1
        chi = 2
        gamma = 1/2
        sigma = 1/2
        c = 0.05
        threshold = 1e-9 % I just put a really low number
    end
    methods
        function [o] = NelderMead(J,n,start,printInfo)
            o.printInfo = printInfo;
            
            o.J = J;
            o.n = n;
            o.x = o.generatePoints(start);
        end
        function [x] = generatePoints(o,start)
            x = zeros([o.n, o.n+1]);
            % set point 1
            x(:,1) = start;
            % randomly generate points 2 through n+1
            for i=2:o.n+1
                x(:,i) = x(:,1);
                % this will create a new vertex c distance in each dimension
                % it doesn't really matter, it can be chosen randomly even.
                x(i-1,i) = x(i-1,i) + o.c;
            end
        end
        % Run the algorithm one step at a time. Good if you want to graph
        % it.
        function [] = step(o)
            o.x = o.sortPoints(o.x);
            xo = o.calculateCentroid( o.x(:, 1:end-1) );
            x1 = o.x(:, 1);
            xn = o.x(:, end-1);
            xnPlus1 = o.x(:, end);
            xr = o.reflectPointAboutPoint(xnPlus1, xo);
            if(o.J(xr) < o.J(x1))
                if o.printInfo, disp("Expansion"), end
                xe = xo + o.chi * (xr - xo);
                if (o.J(xe) < o.J(xr))
                    o.x(:,end) = xe;
                else
                    o.x(:,end) = xr;
                end
                return;
            elseif(o.J(xr) <= o.J(xn))
                if o.printInfo, disp("Replace"), end
                o.x(:, end) = xr;
                return;
            else
                if(o.J(xr) > o.J(xnPlus1))
                    xc = o.gamma*xnPlus1 + (1 - o.gamma)*xo;
                    if(o.J(xc) < o.J(xnPlus1))
                        if o.printInfo, disp("Inside Contraction"), end
                        o.x(:,end) = xc;
                        return;
                    end
                else
                    xc = xo + o.gamma * (xr - xo);
                    if(o.J(xc) < o.J(xnPlus1))
                        if o.printInfo, disp("Outside Contraction"), end
                        o.x(:,end) = xc;
                        return;
                    end
                end
                if o.printInfo, disp("Shrink"), end
                x1 = o.x(:,1);
                o.x = x1 + o.sigma * (o.x - x1);
                return;
            end
        end
        function [minPos] = getCurrentMinimum(o)
            p = o.sortPoints(o.x);
            minPos = p(:,1);
        end
        % Find the minimum of the function
        function [minPos] = minimize(o)
            while(true)
                if(o.shouldTerminate())
                    break;
                end
                o.step();
            end
            minPos = o.getCurrentMinimum();
            return
        end
        function [vals] = getValues(o, p)
            sz = size(p);
            vals = zeros([1,sz(2)]);
            for i=1:sz(2)
                vals(i) = o.J(p(:,i));
            end
        end
        function [ordered] = sortPoints(o, p)
            sz = size(p);
            vals = o.getValues(p);
            combined = [p;vals]';
            combinedSorted = sortrows(combined, sz(1)+1)';
            ordered = combinedSorted(1:end-1, :);
        end
        function [terminate] = shouldTerminate(o)
            std_dev = std(o.getValues(o.x)); % get the standard deviation
            if(std_dev < o.threshold)
                terminate = true;
                if(o.printInfo)
                    disp("Std_dev=" + std_dev + ". Should Terminate");
                end
            else
                terminate = false;
            end
        end
        function [centroid] = calculateCentroid(~,p)
            sz = size(p);
            dims = sz(2);
            op = ones([1, dims]) / dims;
            centroid = (op * p')';
        end
        function [reflected] = reflectPointAboutPoint(o, p1, p2)
            reflected = p2 + o.rho*(p2 - p1);
        end
    end
end


