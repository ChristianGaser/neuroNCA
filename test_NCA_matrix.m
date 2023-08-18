function test_NCA_matriX

%clf
rng(2)
sc = 20;
m = 5;

X = (1:20*sc)';
Xm = repmat(X,1,m);

X = X(randperm(numel(X)));
Y = 2*sign(randn(1,m)).*Xm + 5*sc*randn(size(Xm));
%Y = 2*Xm + 10*sc*randn(size(Xm));

[d, p, Xfdh, Yfdh] = NCA(X,Y);
beta = linear_fit(X,Y);

i = 5;

d = d(i);
Yfdh = Yfdh(:,i);
ind = isfinite(Yfdh);
Xfdh = Xfdh(ind);
Yfdh = Yfdh(ind);
p = p(i,:);

x = X;
y = Y(:,i);

Xmin = min(Xfdh);
Xmax = max(Xfdh);

figure(11)
scatter(x,y,'g')
hold on

Xfit = linspace(Xmin,Xmax,100);
Yfit = p(1) + p(2)*Xfit;

fprintf('Effect size of necessary condition: %g\n',d);
plot(Xfdh,Yfdh,'r.')
plot(Xfit,Yfit);
hold off
legend('Data','CE-FDH','CR-FDH')

% ______________________________________________________________________
function beta = linear_fit(X,Y)
% Estimate offset and slope for the data Y = f(X) using linear fit
%
% Format: beta = linear_fit(X,Y)
% X    - vector of determinant of size n x 1
% Y    - matrix of outcome of size m x n
%
% Output:
% beta - offset and slope of y with size m x 2

n = numel(X);

[my, ny] = size(Y);
[mx, nx] = size(X);

% transpose if necessary
if n ~= ny
  Y = Y';
end

% transpose if necessary
if n ~= nx
  X = X';
end

G = [ones(size(X)); X];

% if Y has no NaNs we can estimate beta as matrix
if all(isfinite(Y(:)))
  beta = Y*pinv(G);
else
  % otherwise we have to find finite entries and limit the values and go 
  % through each row
  beta = zeros(size(Y,1),2);
  for i=1:size(Y,1)
    ind = isfinite(Y(i,:));
    beta(i,:) = Y(i,ind)*pinv(G(:,ind));
  end
end