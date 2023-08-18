function demo_NCA
% demo tool for univariate NCA
% ______________________________________________________________________
%
% Christian Gaser
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________

rng(2)

% simulate y with a slope of 2
x = 1:100;
y =  10 + 2*x + 50*randn(size(x));

y =  [10 + 2*x + 50*randn(size(x)); 5 + 2.5*x + 50*randn(size(x))];

[d, p, Xfdh, Yfdh] = NCA_switched(x,y);

d
return
Xmin = min(x);
Xmax = max(x);

figure(1)
clf
subplot(2,1,1)
scatter(x,y,'g')
hold on

Xfit = linspace(Xmin,Xmax,100);
Yfit = p(1) + p(2)*Xfit;

Acc = 100*(1 - sum(y > (p(1) + p(2)*x))/numel(x));

str = 'Example 1 with large effect size';
fprintf('%s\n',str);
fprintf('Effect size of necessary condition: %g\n',d);
fprintf('Accuracy of the ceiling line: %g\n',Acc);

plot(Xfdh,Yfdh,'r.')
plot(Xfit,Yfit);
hold off
legend('Data','CE-FDH Data','CR-FDH Fit','Location','NorthWest')
title(str)

fprintf('\n');
% 2nd example with more added noise and smaller effect size
%y = 10 + 2*x + 50*randn(size(x));

subplot(2,1,2)

y = repmat(y,2,1);

if 0
x0 = y;
y0 = x;
x = x0; y = y0;
[d, p, Xfdh, Yfdh] = NCA(x,y);
Xmin = min(x);
Xmax = max(x);
scatter(x,y,'g')
else

n = size(y,1);
d = zeros(n,1);
p = zeros(n,2);
Xfdh = nan(size(x));
Yfdh = nan(size(x,2),size(x,1),'single');
tic
for i=1:n
  [d(i), p(i,:), Xfdh, Yfdh(:,i)] = NCA_vector(x,y(i,:));
%  [d(i), p(i,:),Yfdh(:,i), Xfdh] = NCA_vector(y(i,:),x);
end
toc
tic
[d, p, Xfdh, Yfdh] = NCA(x,y);
toc

d
Xmin = min(y(1,:));
Xmax = max(y(1,:));
scatter(y(1,:),x,'g')
end

hold on

Xfit = linspace(Xmin,Xmax,100);
Yfit = p(1,1) + p(1,2)*Xfit;

Acc = 100*(1 - sum(y > (p(1) + p(2)*x))/numel(x));

str = 'Example 2 with small effect size';
fprintf('%s\n',str);
fprintf('Effect size of necessary condition: %g\n',d);
%fprintf('Accuracy of the ceiling line: %g\n',Acc);

plot(Yfdh,Xfdh,'r.')
plot(Yfit,Xfit);
hold off
legend('Data','CE-FDH Data','CR-FDH Fit','Location','NorthWest')
title(str)
