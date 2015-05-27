function [ totalEnergy ] = calculateDipoleEnergy( latticeHeight,...
    latticeWidth,positionMatrix,dipoleMomentMatrix )
%calculateDipoleEnergy does as the name implies
%   This function calculates the energy associated with a crystal of
%   dipoles.  It does this by repeatedly calling individualDipoleEnergy

individualEnergyMatrix=zeros(latticeHeight,latticeWidth);
totalEnergy=0;

% Construct matrix of total energies associated with each dipole.

k=0;
for k=1:latticeHeight
    f=0;
    for f=1:latticeWidth
        individualEnergyMatrix(k,f)=individualDipoleEnergy( latticeHeight,...
    latticeWidth,positionMatrix,dipoleMomentMatrix,k,f);
    end
end

% Sum over all dipoles.

totalEnergy=sum(individualEnergyMatrix(:));

end

function [ individualEnergy ] = individualDipoleEnergy( latticeHeight,...
    latticeWidth,positionMatrix,dipoleMomentMatrix,k,f)
% individualDipoleEnergy calculates its namesake
%   The function calculates the energy associated with a single dipole's
%   (the one at row k, column f) interaction with the rest of the dipoles
%   in the crystal.

% Initialize variables

epsilon=8.854187817e-12;

rVectorMatrix=zeros(latticeHeight,latticeWidth,2);
r1Matrix=zeros(latticeHeight,latticeWidth,2);
rScalarMatrix=zeros(latticeHeight,latticeWidth);
normedRVectorMatrix=zeros(latticeHeight,latticeWidth,2);
rScalarIntermed=zeros(latticeHeight,latticeWidth,2);
energyMatrix=zeros(latticeHeight,latticeWidth);
p1Matrix=zeros(latticeHeight,latticeWidth,2);

% Construct a matrix containing the displacement vector from every dipole
% to the target dipole

r1Matrix(:,:,1)=positionMatrix(k,f,1);
r1Matrix(:,:,2)=positionMatrix(k,f,2);

rVectorMatrix=r1Matrix-positionMatrix;

% Construct matrices containing the magnitude of the displacement vector
% and a normed displacement vector matrix

rScalarIntermed=rVectorMatrix.^2;
rScalarMatrix(:,:)=sqrt(sum((rScalarIntermed),3));

normedRVectorMatrix(:,:,1)=rVectorMatrix(:,:,1)./rScalarMatrix;
normedRVectorMatrix(:,:,2)=rVectorMatrix(:,:,2)./rScalarMatrix;

% Construct the matrix of interaction energies

p1Matrix(:,:,1)=dipoleMomentMatrix(k,f,1);
p1Matrix(:,:,2)=dipoleMomentMatrix(k,f,2);

energyMatrix=(1/(8*pi*epsilon))*(sum((p1Matrix.*dipoleMomentMatrix),3)-...
    3*(sum((p1Matrix.*normedRVectorMatrix),3)).*...
    (sum((dipoleMomentMatrix.*normedRVectorMatrix),3)))./(rScalarMatrix.^3);

energyMatrix(k,f)=0;

% Sum over the matrix to find the energy associated with the k,fth dipole

individualEnergy=sum(energyMatrix(:));

end

