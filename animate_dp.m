function [] = animate_dp(theta1, theta2, filename)
% animate and record the movement of the musculoskeletal model
global ground humerus urh hand ...

gh = [ground; humerus];   %  hu = [humerus; urh];

figure('Position', [500 500 1020 1020])
v = VideoWriter('results/'+filename);
v.FrameRate = 100;       %v.Quality = 100;
open(v)

[x, y, z] = sphere;
alpha 0.15; 


for i = 1:length(theta1)

    J1 = comp_rot3(theta1(i))';
    J2 = comp_rot3(theta2(i))';

    urh_g = urh*J1;
    hu = [humerus; urh_g];
    hand_g = (hand - urh)*J2 + urh_g;
    uh = [urh_g; hand_g];
    
    surf(x*0.03,                        y*0.03,                         z*0.03, 'FaceColor', 'k')
    hold on, xlim([-0.6 0.6]); ylim([-0.6 0.6])
    xlabel('X'), ylabel('Y'), zlabel('Z')
    view([0 90]);
    % joints
    surf(x*0.03+humerus(1),             y*0.03+humerus(2),              z*0.03, 'FaceColor', 'b', 'EdgeColor', 'none')
    surf(x*0.015+urh_g(1),              y*0.015+urh_g(2),               z*0.015,'FaceColor', 'r', 'EdgeColor', 'none')
    surf(x*0.02+hand_g(1),              y*0.02+hand_g(2),               z*0.03, 'FaceColor', 'g', 'EdgeColor', 'none')
    
    % bones
    plot3(gh(:, 1), gh(:, 2), zeros(2, 1), 'k', 'LineWidth', 10)
    plot3(hu(:, 1), hu(:, 2), zeros(2, 1), 'k', 'LineWidth', 10)
    plot3(uh(:, 1), uh(:, 2), zeros(2, 1), 'k', 'LineWidth', 10)
    hold off
%     
    % record frame
    frame = getframe(gcf);
    set(gcf, 'Units', 'pixels', 'Position', [500 500 1020 1020])             % to avoid reshaping of first frame
    frame = getframe(gcf);
    writeVideo(v, frame)
    i
end
close(v)
end