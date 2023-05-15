function [d, p, Xfdh, Yfdh] = NCA(X,Y)
% Necessary Condition Analysis (NCA) for identifying necessary (but not
% sufficient) conditions in data sets with discrete (with many values) 
% and continuous necessary conditions.
% The ceiling line is estimated using ceiling regression with Free Disposal
% Hull (FDH) to obtain the "empty space" of the left upper corner.
% If a negative slope is found for Y is found then the "empty space" of left 
% lower corner is estimated by mirroring Y.
%
% For a detailed introduction into NCA, see the comprehensive documentation at:
% https://www.erim.eur.nl/necessary-condition-analysis/
%
% Reference:
% Dul, J. (2016). Necessary Condition Analysis (NCA): Logic and methodology 
% of "necessary but not sufficient" causality. Organizational Research 
% Methods, 19(1), 10-52.
%
% Format: [d, beta] = NCA(X,Y)
% X    - vector of determinant of size n x 1
% Y    - matrix of outcome of size n x m
%
% Output:
% d    - effect size of the necessary condition with size m x 1
% p    - coefficient s of the linear fit: offset and slope of the ceiling line
%        that defines the "empty space" with size m x 2
% Xfdh - estimated X-values of the Free Disposal Hull (FDH)
% Yfdh - estimated Y-values of the Free Disposal Hull (FDH)
% ______________________________________________________________________
%
% Christian Gaser
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________

% transpose if necessary
if size(X,1) == 1
  X = X';
end

  % transpose if necessary
if size(X,1) ~= size(Y,1)
  Y = Y';
end

% check size
if size(X,1) ~= size(Y,1)
  error('Size mismatch between X and Y.');
end

% check whether slope is negative
beta = linear_fit(X,Y);

% if slope is negative we have to use -Y and rather estimate the left
% lower traingle
if beta(:,2) < 0
  flipy = true;
  Y = -Y;
else
  flipy = false;
end

Xfdh = [];
Yfdh = [];
[Xsort, Xind] = sort(X);

% initialize maximum of Yi
YimX = -1e10*ones(1,size(Y,2));

% go through each X-value
for j = 1:numel(X)
  
  % find entrY in sorted X-values
  i = find(X == Xsort(j));
  
  % if we have multiple Y-values for the same X-value we need the largest
  % Y-value
  if numel(i) > 1
    Yi = max(Y(i,:));
    i = i(1);
  else
    Yi = Y(i,:);
  end
  
  % if Yi exceeds previous maximum we add the new values to the fdh-list
  ind = Yi > YimX;
  if any(ind)
    YimX = Yi;
    Xfdh = [Xfdh; X(Xind(i))];
    Yfdh = [Yfdh; YimX];
  end
end

if flipy
  Yfdh = -Yfdh;
end

p = linear_fit(Xfdh,Yfdh);

Xmin = min(X);
Xmax = max(X);
Ymin = min(Y);
Ymax = max(Y);

% estimate effect size
S = (Xmax - Xmin) * (Ymax - Ymin); % whole area
C = (Xmax - Xmin) * (beta(:,2)*Xmax - beta(:,2)*Xmin) * 0.5; % left upper triangle
d = abs(C./S');


% ______________________________________________________________________
function beta = linear_fit(X,Y)
% Estimate offset and slope for the data Y = f(X) using linear fit
%
% Format: beta = linear_fit(X,Y)
% X    - vector of determinant of size n x 1
% Y    - matrix of outcome of size n x m
%
% Output:
% beta - offset and slope of y with size m x 2


% transpose if necessary
if size(X,1) ~= 1
  X = X';
end

% check size
if size(X,1) ~= size(Y,1)
  Y = Y';
end

% add constant to X
G = [ones(size(X)); X];

% estimate beta
beta = Y*pinv(G);