% fig=open('No6trials/PFC_correct_trajectories.fig');
fig=open('No6trials/PPC_alldistpoints.fig');

% % patches
% % Step 1: Create a 3D figure
% subplot(1,2,1);
% hold on;
% pc1_points=[-10 13]; pc2_points=[-12 9];
% plot_patches(pc1_points,pc2_points, [0 0.5],'r');
% plot_patches(pc1_points,pc2_points,[2 2.5],'r');
% plot_patches(pc1_points,pc2_points,[4 5], 'b');
% hold on
% xlim(pc1_points);
% ylim(pc2_points);
% zlim([-1 4.5]);
% view([-5 5]);
% zlabel('Time(s)')
% % xlim([-4 4]);
% % ylim([-4 4]);
% view ([-190 5]);
% 
% subplot(1,2,2);
% hold on;
% 
% pc1_points=[-8 12]; pc2_points=[-9 8];
% plot_patches(pc1_points,pc2_points, [0 0.5],'r');
% plot_patches(pc1_points,pc2_points,[2 2.5],'r');
% plot_patches(pc1_points,pc2_points,[4 5], 'b');
% 
% xlim(pc1_points);
% ylim(pc2_points);
% zlim([-1 4.5]);
% view([-5 5]);
% zlabel('Time(s)')
% % xlim([-4 4]);
% % ylim([-4 4]);
% view ([-190 5]);
set(fig, 'Renderer', 'painters');
exportgraphics(fig, 'No6trials/PPC_alldistpoints.emf', 'ContentType', 'vector');
print(fig, 'No6trials/PPC_alldistpoints.pdf', '-vector', '-bestfit', '-dwinc');
% print('-depsc', 'No6trials/PPC_avg_corr_err_trajectories.eps');


% function plot_patches(pc1_points,pc2_points,time, Color)
% vertices = [...
%     pc1_points(1), pc2_points(1), time(1);   % Vertex 1 (bottom-left, front)
%     pc1_points(2), pc2_points(1), time(1);   % Vertex 2 (bottom-right, front)
%     pc1_points(2), pc2_points(1), time(2);   % Vertex 3 (top-right, front)
%     pc1_points(1), pc2_points(1), time(2);   % Vertex 4 (top-left, front)
%     pc1_points(1), pc2_points(2), time(1);   % Vertex 5 (bottom-left, back)
%     pc1_points(2), pc2_points(2), time(1);   % Vertex 6 (bottom-right, back)
%     pc1_points(2), pc2_points(2), time(2);   % Vertex 7 (top-right, back)
%     pc1_points(1), pc2_points(2), time(2);   % Vertex 8 (top-left, back)
% ];
% 
% faces = [...
%     1, 2, 3, 4;  % Face 1 (front)
%     2, 6, 7, 3;  % Face 2 (right)
%     6, 5, 8, 7;  % Face 3 (back)
%     5, 1, 4, 8;  % Face 4 (left)
%     4, 3, 7, 8;  % Face 5 (top)
%     1, 5, 6, 2;  % Face 6 (bottom)
% ];
% 
% % Step 3: Colorize the rectangular patch
% 
% patch1=patch('Vertices', vertices, 'Faces', faces, 'FaceColor', Color,'EdgeColor', Color);
% set(patch1, 'FaceAlpha', 0.05, 'EdgeAlpha', 0.05);
% hold on
% % grid on;
% end