function Frac2vtk(fname,physdim,G,fl)
%GRID2VTK Summary of this function goes here
%Detailed explanation goes here

if(numel(physdim)==3)
    DZ=physdim(3);
else
    DZ=G.fractureHeight;
end

NumFracs=size(fl,1);
pts=[fl(:,1:2) 0.*ones(NumFracs,1) fl(:,3:4) 0.*ones(NumFracs,1)];
pts=[pts; [fl(:,3:4) DZ.*ones(NumFracs,1) fl(:,1:2) DZ.*ones(NumFracs,1) ]];
NumPts=size(pts,1)*2;

FracPolyId=[1:NumPts/2];
FracPolyId=reshape(FracPolyId,[2,NumFracs])';
FracPolyId=[FracPolyId FracPolyId+FracPolyId(end)]-1;

%Hydraulic Fractures tag=0, Natural Fracs=1
FracTags=ones(NumFracs,1);
FracTags(1:G.NumHFs)=0;

fid = fopen(fname, 'w'); 
fprintf(fid, '# vtk DataFile Version 2.0\n');
fprintf(fid, 'Fracture Geometry\n');
fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET POLYDATA\n');
fprintf(fid, 'POINTS %d float\n', NumPts);
fprintf(fid, '%0.10f %0.10f %0.10f %0.10f %0.10f %0.10f\n', pts');fprintf(fid, '\n');
fprintf(fid, 'POLYGONS %d %d\n', NumFracs,(4+1)*NumFracs);
fprintf(fid, '4 %d %d %d %d\n', FracPolyId');fprintf(fid, '\n');
fprintf(fid, 'CELL_DATA %d\n', NumFracs);
fprintf(fid, ['SCALARS ', 'FractureTag', ' int\n']);
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid, '%d %d %d %d %d %d\n', FracTags');
fclose(fid);

end

