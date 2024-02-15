clc
x = linspace(-1.5, 1.5, 51);
y = linspace(-2.5, 2.5, 51*2+1);
[X, Y] = meshgrid(x, y);
points = [X(:), Y(:)];
scatter(points(:, 1), points(:, 2), 'o');
axis equal;
hold on;
selectedPoints = [];
setappdata(gcf, 'pointsData', points);
set(gcf, 'WindowButtonDownFcn', @mouseClickCallback);