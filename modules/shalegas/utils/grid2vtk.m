function grid2vtk(fname,G,dataName,dataArray)
%GRID2VTK Summary of this function goes here
%Detailed explanation goes here

unit_x=ft;
unit_y=psia;

nx=G.cartDims(1);ny=G.cartDims(2);
if(numel(G.cartDims)==3)
    nz=G.cartDims(3);
else
    nz=1;
end
coord_x=G.nodes.coords(1:nx+1,1);
coord_y=G.nodes.coords(1:nx+1:(nx+2)*ny,2);

fid = fopen(fname, 'w'); 
fprintf(fid, '# vtk DataFile Version 2.0\n');
fprintf(fid, 'Rectilinear grid\n');
fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET RECTILINEAR_GRID\n');
fprintf(fid, 'DIMENSIONS %d %d %d\n',nx+1,ny+1,nz+1);
fprintf(fid, 'X_COORDINATES %d float\n', nx+1);
fprintf(fid, '%0.10f ', coord_x(:)');fprintf(fid, '\n');
fprintf(fid, 'Y_COORDINATES %d float\n', ny+1);
fprintf(fid, '%0.10f ', coord_y(:)'); fprintf(fid, '\n');
fprintf(fid, 'Z_COORDINATES %d float\n', nz+1);
fprintf(fid, '%0.10f ', [0 coord_x(end)/10.0]);fprintf(fid, '\n');
fprintf(fid, 'CELL_DATA %d\n', nx*ny*nz);
fprintf(fid, ['SCALARS ', dataName, ' float\n']);
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid, '%0.10f ', dataArray(1:nx*ny*nz)');
fclose(fid);

end

