function p = polyfix(x,y,n,xfix,yfix,xder,dydx)
%POLYFIX Fit polynomial p to data, but specify value at specific points
%
%   p = polyfix(x,y,n,xfix,yfix) 
% finds the coefficients of the polynomial of degree n that fits the data 
% in a least-squares sense, with the constraint that polyval(p,xfix) = yfix
%
%   p = polyfix(x,y,n,xfix,yfix,xder,dydx)
% uses the additional constraint that the derivative at xder = dydx
% 
% The polynomial order n must be high enough to match all specified values
% NOTE: For the lowest order allowed, p will fit the constraints, but 
%       may disregard x and y.
%
% Example 1:
% x = linspace(0,2,100)';y = sin(pi*x)+ 0.1*randn(100,1);
% p = polyfix(x,y,3,[0,2],[0,0]);plot(x,y,'.',x,polyval(p,x));
% 
% Example 2:
% x = linspace(0,1,100)';y = sin(x*pi/2) + 0.1*randn(100,1);
% p = polyfix(x,y,4,[],[],[0 1],[1 0]);plot(x,y,'.',x,polyval(p,x))
% 
% See also: polyfit, polyval

% Are Mjaavatten, Telemark University College, Norway, November 2015

% Revision history
% 2015-11-28: Version 1.0
% 2015-12-07: Version 1.1:
%             Added option for specifying derivatives
%             The output is now a row vector, for consistency with polyfit
%             Added test for polynomial degree
% 2017-08-23: Version 1.2
%             Fixed trivial errors in the help section

    %% Make sure all input arrays are column vectors of compatible sizes:
    x = x(:);
    y = y(:);
    nfit = length(x);
    if ~(length(y)== nfit)
        error('x and y must have the same size');
    end
    xfix = xfix(:);
    yfix = yfix(:);
    nfix = length(xfix);
    if ~(length(yfix)== nfix)
        error('xfit and yfit must have the same size');
    end
    if nargin > 5 % Derivatives specified
        xder = xder(:);
        dydx = dydx(:);
    else
        xder = [];
        dydx = [];
    end
    nder = length(xder); 
    if ~(length(dydx) == nder)
        error('xder and dydx must have the same size');
    end
    nspec = nfix + nder;
    specval = [yfix;dydx];
 
    %% First find A and pc such that A*pc = specval
    A = zeros(nspec,n+1); 
    % Specified y values
    for i = 1:n+1
        A(1:nfix,i) = ones(nfix,1).*xfix.^(n+1-i);
    end
    % Specified values of dydx
    if nder > 0
        for i = 1:n
            A(nfix +(1:nder),i) = (n-i+1)*ones(nder,1).*xder.^(n-i);
        end
    end

    if nfix > 0
        lastcol = n+1;
        nmin = nspec - 1;
    else
        lastcol = n;   % If only derivatives, p(n+1) is arbitrary
        nmin = nspec;
    end
    if n < nmin
        error(['Polynomial degree too low. Cannot match all constraints']);
    end    
    %% Find the unique polynomial of degree nmin that fits the constraints. 
    firstcol = n-nmin+1;   % A(:,firstcol_lastcol) detrmines p0   
    pc0 = A(:,firstcol:lastcol)\specval;  % Satifies A*pc = specval
    % Now extend to degree n and pad with zeros:
    pc = zeros(n+1,1);
    pc(firstcol:lastcol) = pc0;    % Satisfies A*pcfull = yfix
    
    % Column i in matrix X is the (n-i+1)'th power of x values 
    X = zeros(nfit,n+1);
    for i = 1:n+1
        X(:,i) = ones(nfit,1).*x.^(n+1-i);
    end
        
    % Subtract constraints polynomial values from y. 
    yfit = y-polyval(pc,x);
    
    %% We now find the p0 that mimimises (X*p0-yfit)'*(X*p0-yfit)
    %  given that A*p0 = 0
    B = null(A);    % For any (n+1-nspc by 1) vector z, A*B*z = 0     
    z = X*B\yfit;   % Least squares solution of X*B*z = yfit
    p0 = B*z;       % Satisfies A*p0 = 0;
    p = p0'+pc';    % Satisfies A*p = b; 
end