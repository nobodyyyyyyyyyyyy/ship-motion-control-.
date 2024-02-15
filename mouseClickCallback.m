function mouseClickCallback(~, ~)
    points = getappdata(gcf, 'pointsData');
    selectedPoints = getappdata(gcf, 'selectedPoints');
    currentPoint = get(gca, 'CurrentPoint');
    mouseX = currentPoint(1, 1);
    mouseY = currentPoint(1, 2);
    distances = sqrt((points(:, 1) - mouseX).^2 + (points(:, 2) - mouseY).^2);
    [~, index] = min(distances);
    mouseX = points(index, 1);
    mouseY = points(index, 2);
    isPointSelected = false;
    for i = 1:size(selectedPoints, 1)
        if isequal(selectedPoints(i, :), [mouseX, mouseY])
            selectedPoints(i, :) = [];
            isPointSelected = true;
            break;
        end
    end
    if ~isPointSelected
        selectedPoints = [selectedPoints; mouseX, mouseY];
    end
    clf;
    scatter(points(:, 1), points(:, 2), 'o');
    hold on;
    scatter(selectedPoints(:, 1), selectedPoints(:, 2), 50, 'filled', 's', 'MarkerEdgeColor', 'r', 'LineWidth', 2);
    axis equal;
    setappdata(gcf, 'selectedPoints', selectedPoints);
    disp('选中点坐标:');
    disp(selectedPoints);
end