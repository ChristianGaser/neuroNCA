function demo_NCA

rng(2)

% simulate y with a slope of 2
x = 1:100;
y =  10 + -2*x + 50*randn(size(x));

[d, p, Xfdh, Yfdh] = NCA(x,y);

Xmin = min(x);
Xmax = max(x);

figure
clf
subplot(2,1,1)
scatter(x,y,'g')
hold on

Xfit = linspace(Xmin,Xmax,100);
Yfit = p(1) + p(2)*Xfit;

Acc = 100 - sum(y > (p(1) + p(2)*x))/numel(x);

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
y = 10 + 2*x + 200*randn(size(x));
[d, p, Xfdh, Yfdh] = NCA(x,y);

subplot(2,1,2)
scatter(x,y,'g')
hold on

Xfit = linspace(Xmin,Xmax,100);
Yfit = p(1) + p(2)*Xfit;

Acc = 100 - sum(y > (p(1) + p(2)*x))/numel(x);

str = 'Example 2 with small effect size';
fprintf('%s\n',str);
fprintf('Effect size of necessary condition: %g\n',d);
fprintf('Accuracy of the ceiling line: %g\n',Acc);

plot(Xfdh,Yfdh,'r.')
plot(Xfit,Yfit);
hold off
legend('Data','CE-FDH Data','CR-FDH Fit','Location','NorthWest')
title(str)
