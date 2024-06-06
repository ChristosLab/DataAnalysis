fig = open('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\Figure3\NS_BS_histogram20231201.fig');
% ylim([0.4 0.6])
set(fig, 'Renderer', 'painters');
exportgraphics(fig, 'C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\Figure3\NS_BS_histogram20231201.emf', 'ContentType', 'vector');
print(fig, 'C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\Figure3\NS_BS_histogram20231201.pdf', '-vector', '-bestfit', '-dwinc');


