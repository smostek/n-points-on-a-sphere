classdef thomson < shape
methods

function s = thomson(n,tol,cutoff)
    % Solve the Thomson problem for a given value of n
    %
    % The Thomson problem is a subset of my larger problem that looks for
    % the lowest-energy distribution of electrons on the sphere. (Or at
    % least, it is assumed to be a subset because it is assumed that every
    % Thomson solution is in equilibrium)
    % 
    % This program works because of the assumption that the Thomson
    % distribution is the only STABLE equilibrium (there are exceptions, notably 
    % at 16, but it holds way more than not) which means the hilariously simple 
    % minimization technique of moving each point a small step in the direction 
    % of the force acting on it (almost) always leads to the Thomson solution.
    %
    % There does exist the complication of how large of a step to take, but
    % after much data-combing I was able to split that problem into two
    % regimes with step sizes that depend only on n.
    arguments
        n double {mustBeInteger}
        % Number of vertices of solution

        tol = 1e-8 
        % Value of CS at which minimization stops

        cutoff {mustBeInteger} = 10
        % Number of attempts to make before giving up
    end
    if nargin >0
        if n > 100
            warning('Really high values are unlikely to work')
        end
        tries = 0;
        ohno = true;
        fprintf('Searching for Solutions:\n')
        while ohno && tries <cutoff
            % generate random initial position:
            C=coorConvert([ones(1,n); asin(2*rand(1,n)-1); pi*(2*rand(1,n)-1)],'sph');

            %initial coarse minimization:
            [C,ohno] = simpleMin(C,0.5,10^(-0.971*n^(0.214)));
            if ~ohno
            %secondary refining minimization:
                if n > 5
                    pow = -0.45*acosh(n-2.6);
                else
                    pow = -0.5;
                end
                [C,ohno] = simpleMin(C,tol,10^pow);
            end
            tries = tries +1;
            if ohno
                fprintf('Attempt %d failed :(\n',tries)
            end
        end

        if ~ohno
            fprintf('Attempt %d succeeded!\n',tries)
        end
        s.cartCoor = C;
        s.minOrd = -log10(CS(C));
    end

        function [C,flag] = simpleMin(C,tol,dt)
            [cs,T]=CS(C);
            mincs = cs;
            while cs >tol && cs < 2*mincs
                C = C+T*dt;
                [cs,T]=CS(C);
                if cs <mincs
                    mincs = cs;
                end
            end
            flag = cs>tol; % 0 means successful minimization
        end
end

end
end