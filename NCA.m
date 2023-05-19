function [d, p] = NCA(X,Y)
% Necessary Condition Analysis (NCA) for identifying necessary (but not
% sufficient) conditions in data sets with discrete (with many values) 
% and continuous necessary conditions.
% The ceiling line is estimated using ceiling regression with Free Disposal
% Hull (FDH) to obtain the "empty space" of the left upper corner.
% If a negative slope is found for Y is found then the "empty space" of right 
% upper corner is estimated by mirroring X.
%
% For a detailed introduction into NCA, see the comprehensive documentation at:
% https://www.erim.eur.nl/necessary-condition-analysis/
% https://bookdown.org/ncabook/advanced_nca2/
%
% Reference:
% Dul, J. (2016). Necessary Condition Analysis (NCA): Logic and methodology 
% of "necessary but not sufficient" causality. Organizational Research 
% Methods, 19(1), 10-52.
%
% Format: [d, beta] = NCA(X,Y)
% X    - vector of determinant of size n x 1
% Y    - matrix of outcome of size m x n
%
% Output:
% d    - effect size of the necessary condition with size m x 1
% p    - coefficient s of the linear fit: offset and slope of the ceiling line
%        that defines the "empty space" with size m x 2
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
  
% initialize output
d = zeros(size(Y,2),1);
p = zeros(size(Y,2),2);

% check whether slope is negative
beta = linear_fit(X,Y);

% if slope is negative we have to mirror Y and estimate the left lower corner
flipx = beta(:,2) < 0;

Xfdh = nan(size(X));
Yfdh = nan(size(Y,2),size(Y,1),'single');

% go through all entries with a positive slope and estimate
% the empty space of the left upper triangle
if any(~flipx)
  [Xsort, Xind] = sort(X);
  
  Ypos = Y(:,~flipx);
  
  % initialize maximum of Yi
  Yimx = -1e10*ones(1,size(Ypos,2));
  
  % go through each X-value
  for j = 1:numel(X)
    
    % find entry in sorted X-values
    i = find(X == Xsort(j));
    
    % if we have multiple Y-values for the same X-value we need the largest
    % Y-value
    if numel(i) > 1
      Yi = max(Ypos(i,:));
      i = i(1);
    else
      Yi = Ypos(i,:);
    end
    
    % if Yi exceeds previous maximum we add the new values to the fdh-peers
    ind = Yi > Yimx;
    if any(ind)
      % only set new maximum for entries where values exceeded
      Yimx(ind) = Yi(ind);
      Yi(~ind) = NaN;
      
      % build list of FDH peers
      Xfdh(Xind(i)) = X(i);
      Yfdh(~flipx,Xind(i)) = Yi;
    end
  end
  
  % get offset and slope of the FDH peers
  [d0, p0] = NCA_effect_size(Xfdh,Yfdh(~flipx,:),X,Ypos);

  d(~flipx) = d0;
  p(~flipx,:) = p0;
end

% and go through all entries with a negative slope and mirror X to estimate
% the empty space of the right upper triangle
if any(flipx)
  X = flipud(X);
  [Xsort, Xind] = sort(X,'descend');
  
  Yneg = flipud(Y(:,flipx));
  
  % initialize maximum of Yi
  Yimx = -1e10*ones(1,size(Yneg,2));
  
  % go through each X-value
  for j = 1:numel(X)
    
    % find entry in sorted X-values
    i = find(X == Xsort(j));
    
    % if we have multiple Y-values for the same X-value we need the largest
    % Y-value
    if numel(i) > 1
      Yi = max(Yneg(i,:));
      i = i(1);
    else
      Yi = Yneg(i,:);
    end
    
    % if Yi exceeds previous maximum we add the new values to the fdh-peers
    ind = Yi > Yimx;
    if any(ind)
      % only set new maximum for entries where values exceeded
      Yimx(ind) = Yi(ind);
      Yi(~ind) = NaN;
      
      % build list of FDH peers
      Xfdh(Xind(i)) = X(i);
      Yfdh(flipx,Xind(i)) = Yi;
    end
  end
  
  % get offset and slope of the FDH peers
  [d0, p0] = NCA_effect_size(Xfdh,Yfdh(flipx,:),X,Yneg);

  d(flipx) = d0;
  p(flipx,:) = p0;
end

% ______________________________________________________________________
function [d, beta] = NCA_effect_size(Xfdh,Yfdh,X,Y)
% Estimate effects size for NCA and offset and slope for the data Y = f(X) 
% using linear fit
%
% Area of the "empty space" is estimated using the ceiling which is explained 
% here:
% https://bookdown.org/ncabook/advanced_nca2/mathematical.html#fig:math
%
% Format: [d, beta] = NCA_effect_size(X,Y)
% Xfdh - estimated X-values (peers) of the Free Disposal Hull (FDH)
%        with size n x 1
% Yfdh - estimated Y-values (peers) of the Free Disposal Hull (FDH)
%        with size m x n
% X    - vector of original X-values of size n x 1
% Y    - matrix of Y outcome of size m x n
%
% Output:
% d    - effect size of NCA with size m x 1
% beta - offset and slope of y with size m x 2

n = numel(Xfdh);

[my, ny] = size(Yfdh);
[mx, nx] = size(Xfdh);

% transpose if necessary
if n ~= ny
  Yfdh = Yfdh';
end

% transpose if necessary
if n ~= nx
  Xfdh = Xfdh';
end

G = [ones(size(Xfdh)); Xfdh];

Xmin = min(X);
Xmax = max(X);
Ymin = min(Y);
Ymax = max(Y);

% estimate whole area
S = (Xmax - Xmin) * (Ymax - Ymin);

% if Y has no NaNs we can estimate beta as matrix
if all(isfinite(Yfdh(:)))
  beta = Yfdh*pinv(G);
  ind_pos = beta(:,2) > 0;
  
  % estimate area of left (or right) upper triangle
  Xcmax = (Ymax - beta(:,1))/beta(:,2);
  Xcmax = max([Xmin Xcmax]);
  
  Ycmin(ind_pos)  = beta(:,1) + beta(:,2)*min(Xfdh);
  Ycmin(~ind_pos) = beta(:,1) + beta(:,2)*max(Xfdh);
  C(ind_pos)  = (Xcmax - Xmin) * (Ymax - Ycmin) * 0.5;
  C(~ind_pos) = (Xcmax - Xmax) * (Ymax - Ycmin) * 0.5;

  d = C./S';
else
  % otherwise we have to find finite entries and limit the values and go 
  % through each row
  beta = zeros(size(Yfdh,1),2);
  d = zeros(size(Yfdh,1),1);
  for i=1:size(Yfdh,1)
    ind = isfinite(Yfdh(i,:));
    if any(ind)
      beta(i,:) = Yfdh(i,ind)*pinv(G(:,ind));

      % estimate area of left (or right) upper triangle
      Xcmax = (Ymax(i) - beta(i,1))/beta(i,2);
      Xcmax = max([Xmin Xcmax]);
      if (beta(i,2)) > 0
        Ycmin = beta(i,1) + beta(i,2)*min(Xfdh(ind));
        C = (Xcmax - Xmin) * (Ymax(i) - Ycmin) * 0.5;
      else
        Ycmin = beta(i,1) + beta(i,2)*max(Xfdh(ind));
        C = (Xcmax - Xmax) * (Ymax(i) - Ycmin) * 0.5;
      end
      d(i) = C/S(i);
    end
  end
end

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

G = [ones(size(X)) X];
beta = (pinv(G)*Y)';
