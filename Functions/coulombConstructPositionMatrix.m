function [ coulombPositionMatrix, dipoleChargeMatrix ] = coulombConstructPositionMatrix( ...
    latticeHeight,latticeWidth,basisVector1,basisVector2,unitCellHeight,...
    unitCellWidth,dipoleUnitCell,dipoleLengthUnitCell)
%constructPositionMatrix does as the name implies
%   The function creates a matrix of the positions in physical
%   space of notional dipole charges given a set of input parameters

% Construct a matrix representing the physical location of lattice points

latticePositionMatrix=zeros(latticeHeight,latticeWidth,2);
k=0;
for k=1:latticeHeight
	f=0;
	for f=1:latticeWidth
		latticePositionMatrix(k,f,1)=((f-1)*basisVector1(1,1) + (k-1)*basisVector2(1,1));
		latticePositionMatrix(k,f,2)=((f-1)*basisVector1(2,1) + (k-1)*basisVector2(2,1));
	end
end

% Construct a matrix of dipole moments by calling custom function
% constructDipoleMomentMatrix

[ dipoleMomentMatrix ] = constructDipoleMomentMatrix( latticeHeight,latticeWidth,...
    unitCellHeight,unitCellWidth,dipoleUnitCell);

% Construct a matrix containing unit vectors pointing in the dipole
% directions
dipoleMagnitudeIntermed=zeros(latticeHeight,latticeWidth,2);
dipoleMagnitudeMatrix=zeros(latticeHeight,latticeWidth);
normedDipoleMatrix=zeros(latticeHeight,latticeWidth,2);

dipoleMagnitudeIntermed=dipoleMomentMatrix.^2;
dipoleMagnitudeMatrix=sqrt(sum(dipoleMagnitudeIntermed,3));
normedDipoleMatrix(:,:,1)=dipoleMomentMatrix(:,:,1)./dipoleMagnitudeMatrix;
normedDipoleMatrix(:,:,2)=dipoleMomentMatrix(:,:,2)./dipoleMagnitudeMatrix;

% Construct latticePositionMatrix, which contains the positions of the notional
% point charges constituting the dipole. latticePositionMatrix(k,f,1,:) contains
% the vector indicating the position of the positive charge of dipole k,f
% while latticePositionMatrix(k,f,2,:) contains the vector indicating the position
% of the negative charge of dipole k,f.

% Construct a matrix containing the dipole lengths

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

% Create dipole length matrix by tiling the unit cell
dipoleLengthMatrix=zeros(newHeight,newWidth);
dipoleLengthMatrix=repmat(dipoleLengthUnitCell,ceil(ratioHeight),ceil(ratioWidth));

% Truncate matrix (if need be) to original dimensions
dipoleLengthMatrix=dipoleLengthMatrix(1:latticeHeight,1:latticeWidth);

% Compute coulombPositionMatrix
coulombPositionMatrix=zeros(latticeHeight,latticeWidth,2,2);

coulombPositionMatrix(:,:,1,1)=latticePositionMatrix(:,:,1)+0.5*...
    dipoleLengthMatrix.*normedDipoleMatrix(:,:,1);
coulombPositionMatrix(:,:,1,2)=latticePositionMatrix(:,:,2)+0.5*...
    dipoleLengthMatrix.*normedDipoleMatrix(:,:,2);
coulombPositionMatrix(:,:,2,1)=latticePositionMatrix(:,:,1)-0.5*...
    dipoleLengthMatrix.*normedDipoleMatrix(:,:,1);
coulombPositionMatrix(:,:,2,2)=latticePositionMatrix(:,:,2)-0.5*...
    dipoleLengthMatrix.*normedDipoleMatrix(:,:,2);

% Calculate dipoleChargeMatrix for energy calculation
dipoleChargeMatrix=dipoleMagnitudeMatrix./dipoleLengthMatrix;

end