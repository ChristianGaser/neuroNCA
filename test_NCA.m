function test_NCA

%clf
rng(2)

x = 1:100;
x = x(randperm(numel(x)));

y = 2*x + 50*randn(size(x));

[d, p, Xfdh, Yfdh] = NCA(x,y);

Xmin = min(x);
Xmax = max(x);

figure(11)
subplot(2,1,1)
scatter(x,y,'g')
hold on

Xfit = linspace(Xmin,Xmax,100);
Yfit = p(1) + p(2)*Xfit;

Acc = 100 - sum(y > (p(1) + p(2)*x))/numel(x);

fprintf('Effect size of necessary condition: %g\n',d);
fprintf('Accuracy of the ceiling line: %g\n',Acc);

plot(Xfdh,Yfdh,'r.')
plot(Xfit,Yfit);
hold off
legend('Data','CE-FDH','CR-FDH')

% add more noise
%y = 2*x + 150*randn(size(x));
y = -y;
[d, p, Xfdh, Yfdh] = NCA(x,y);

subplot(2,1,2)
scatter(x,y,'g')
hold on

Xfit = linspace(Xmin,Xmax,100);
Yfit = p(1) + p(2)*Xfit;

Acc = 100 - sum(y > (p(1) + p(2)*x))/numel(x);

fprintf('Effect size of necessary condition: %g\n',d);
fprintf('Accuracy of the ceiling line: %g\n',Acc);

plot(Xfdh,Yfdh,'r.')
plot(Xfit,Yfit);
hold off
legend('Data','CE-FDH','CR-FDH')
