function [ dipoleMomentMatrix ] = constructDipoleMomentMatrix( latticeHeight,latticeWidth,...
    unitCellHeight,unitCellWidth,dipoleUnitCell)
% constructDipoleMomentMatrix does as the name implies
%   The function constructs a matrix describing the dipole moments of
%   the individual molecules of a crystal given a set of input parameters

% Rework crystal dimensions (if necessary) so that it can be filled by
% tiling the unit cell
ratioHeight=0;
ratioWidth=0;
ratioHeight=latticeHeight/unitCellHeight;
ratioWidth=latticeWidth/unitCellWidth;

if mod(ratioHeight,1)~=0
	newHeight=ceil(ratioHeight)*unitCellHeight;
else
    newHeight=latticeHeight;
end

if mod(ratioWidth,1)~=0
	newWidth=ceil(ratioWidth)*unitCellWidth;
else
    newWidth=latticeWidth;
end

% Create dipole moment matrix for crystal by tiling the unit cell
dipoleMomentMatrix=zeros(newHeight,newWidth,2);
dipoleMomentMatrix=repmat(dipoleUnitCell,ceil(ratioHeight),ceil(ratioWidth));

% Truncate matrix (if need be) to original dimensions
dipoleMomentMatrix=dipoleMomentMatrix(1:latticeHeight,1:latticeWidth,:);

end

