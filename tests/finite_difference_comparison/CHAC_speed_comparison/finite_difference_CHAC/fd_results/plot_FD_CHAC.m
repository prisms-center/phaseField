function [z,c_0,c_1] = plot_FD_CHAC()

% --------------------------------------------------------------------------------------------------------
% Plot the concentration
% --------------------------------------------------------------------------------------------------------

run_num = '001';

domain_info = importdata(['domain_info' run_num '.frt']);

xg = domain_info(1);
yg = domain_info(2);
zg = domain_info(3);

output_num = '00';
c_0 = assemble_MPI_output(['con_r' run_num '_' output_num '_'],'',run_num);
eta_0 = assemble_MPI_output(['eta_r' run_num '_' output_num '_'],'',run_num);

x = assemble_MPI_output(['xgrid' run_num '_' output_num '_'],'',run_num);
y = assemble_MPI_output(['ygrid' run_num '_' output_num '_'],'',run_num);
z = assemble_MPI_output(['zgrid' run_num '_' output_num '_'],'',run_num);

output_num = '01';
c_1 = assemble_MPI_output(['con_r' run_num '_' output_num '_'],'',run_num);
eta_1 = assemble_MPI_output(['eta_r' run_num '_' output_num '_'],'',run_num);

threeD_precipitate_plotter(c_1,xg,yg,zg,1);
drawnow;

vtkwrite('fd_results.vtk','structured_grid',x,y,z,...
    'scalars','c_0',c_0,'scalars','c_1',c_1);


end

function threeD_precipitate_plotter(c,xg,yg,zg,fig_num)

plot_start_x = 1;
plot_end_x = xg;
plot_start_y = 1;
plot_end_y = yg;
plot_start_z = 1;
plot_end_z = zg;

figure(fig_num);
whitebg('w')

ite = 1;
isovL = 0.06;
isovC = 0.06;

[NewY,NewX,NewZ] = meshgrid(plot_start_y:plot_end_y,plot_start_x:plot_end_x,plot_start_z:plot_end_z);

p1 = patch(isosurface(NewY,NewX,NewZ,c(plot_start_x:plot_end_x,plot_start_y:plot_end_y,plot_start_z:plot_end_z),isovL));
p2 = patch(isocaps(NewY,NewX,NewZ,c(plot_start_x:plot_end_x,plot_start_y:plot_end_y,plot_start_z:plot_end_z),isovL));

isocolors(c(plot_start_x:plot_end_x,plot_start_y:plot_end_y,plot_start_z:plot_end_z),p1)
set(p1,'FaceColor','interp','EdgeColor','none')
isonormals(c(plot_start_x:plot_end_x,plot_start_y:plot_end_y,plot_start_z:plot_end_z),p1)
isocolors(c(plot_start_x:plot_end_x,plot_start_y:plot_end_y,plot_start_z:plot_end_z),p2)
set(p2,'FaceColor','interp','EdgeColor','none')

set(gca,'Projection','perspective')
view(3); daspect([1,1,1]);
axis([1 yg 1 xg 1 zg]);box on;
lightangle(45,60)
lighting phong

end
