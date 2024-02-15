function [centroid,distence,angles]=ship_control()
datapoints=[ 
    0    2.5000;
   -0.2400    2.3529;
   -0.4200    2.2059;
   -0.5400    2.0588;
   -0.7200    1.9118;
   -0.9000    1.7157;
   -1.0800    1.5196;
   -1.2600    1.3235;
   -1.3800    1.1275;
   -1.5000    0.8824;
   -1.5000    0.6373;
   -1.5000    0.2941;
   -1.5000   -0.0490;
   -1.5000   -0.5392;
   -1.5000   -0.9314;
   -1.5000   -1.3235;
   -1.5000   -1.7157;
   -1.5000   -2.0098;
   -1.5000   -2.3039;
   -1.5000   -2.5000;
   -1.2000   -2.5000;
   -0.9600   -2.5000;
   -0.7200   -2.5000;
   -0.4800   -2.5000;
   -0.2400   -2.5000;
         0   -2.5000;
         ]
symmetry = [-datapoints(:,1),datapoints(:,2)];
symmetry = flip(symmetry);
datapoints = [datapoints;symmetry];
figure();
scatter(datapoints(:,1),datapoints(:,2),'r.',SizeData=50);
axis equal
hold on;
centroid = mean(datapoints);
plot(centroid(1),centroid(2),'blue','Marker','o');
distence=pdist2(centroid,datapoints);
angles = atan2(datapoints(:, 1) -centroid(1), datapoints(:, 2) - centroid(2)) * 180 / pi;
% close;