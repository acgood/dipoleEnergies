function [ totalEnergy ] = coulombCalculateDipoleEnergy( latticeHeight,...
    latticeWidth,coulombPositionMatrix,dipoleChargeMatrix )
%calculateDipoleEnergy does as the name implies
%   This function calculates the energy associated with a crystal of
%   dipoles.  It does this by repeatedly calling individualDipoleEnergy

coulombIndividualEnergyMatrix=zeros(latticeHeight,latticeWidth,2);
totalEnergy=0;
wonkCount=zeros(2,1);

% Construct matrix of total energies associated with each dipole.

k=0;
for k=1:latticeHeight
    f=0;
    for f=1:latticeWidth
        g=0;
        for g=1:2
            [coulombIndividualEnergyMatrix(k,f,g),wonkCount]=coulombIndividualDipoleEnergy( ...
                latticeHeight,latticeWidth,coulombPositionMatrix,dipoleChargeMatrix,k,f,g...
                ,wonkCount);
        end
    end
end

% Sum over all dipoles.

totalEnergy=sum(coulombIndividualEnergyMatrix(:));

% Give error messages

if wonkCount(1,1)>0
    disp('WARNING: Notional charges occupy the same physical space. ')
    disp('Results may be wonky.')
end

if wonkCount(2,1)>0
    disp('WARNING: Charges of the same sign occupy the same physical space. ')
    disp('The error could be quite serious.')
end

end

function [ individualEnergy,wonkCount ] = coulombIndividualDipoleEnergy( latticeHeight,...
    latticeWidth,coulombPositionMatrix,dipoleChargeMatrix,k,f,g,wonkCount)
% individualDipoleEnergy calculates its namesake
%   The function calculates the energy associated with a single dipole's
%   (the one at row k, column f) interaction with the rest of the dipoles
%   in the crystal.

% Initialize variables

epsilon=8.854187817e-12;

rVectorMatrix=zeros(latticeHeight,latticeWidth,2,2);
r1Matrix=zeros(latticeHeight,latticeWidth,2,2);
rScalarMatrix=zeros(latticeHeight,latticeWidth,2);
rScalarIntermed=zeros(latticeHeight,latticeWidth,2,2);
energyMatrix=zeros(latticeHeight,latticeWidth,2);

% Construct a matrix containing the displacement vector from every dipole
% to the target dipole

r1Matrix(:,:,:,1)=coulombPositionMatrix(k,f,g,1);
r1Matrix(:,:,:,2)=coulombPositionMatrix(k,f,g,2);

rVectorMatrix=r1Matrix-coulombPositionMatrix;

% Construct matrices containing the magnitude of the displacement vector
% and a normed displacement vector matrix

rScalarIntermed=rVectorMatrix.^2;
rScalarMatrix(:,:,:)=sqrt(sum((rScalarIntermed),4));

% Test for wonkiness related to notional charges occupying the same
% physical space
wonkTest=zeros(latticeWidth,2);
wonkTest=any(rScalarMatrix==0);

if sum(wonkTest(:))>=2
    wonkCount(1,1)=wonkCount(1,1)+1;
end

% If charges of the same sign are on top of one another, warn that
% the error could be quite serious.

if sum(wonkTest(:,:,g))>=2
    wonkCount(2,1)=wonkCount(2,1)+1;
end

% Create charge matrices that make the computation more efficient
charge1Matrix=zeros(latticeHeight,latticeWidth,2);
otherChargesMatrix=zeros(latticeHeight,latticeWidth,2);
sign=0;

if g==1
    sign=1;
elseif g==2
    sign=-1;
end

charge1Matrix(:,:,:)=sign*dipoleChargeMatrix(k,f);
otherChargesMatrix(:,:,1)=dipoleChargeMatrix;
otherChargesMatrix(:,:,2)=-dipoleChargeMatrix;

energyMatrix=(charge1Matrix.*otherChargesMatrix./rScalarMatrix)/(8*pi*epsilon);

% No self-interaction or charges on top of one another
energyMatrix(k,f,g)=0;
energyMatrix(energyMatrix==inf)=0;
energyMatrix(energyMatrix==-inf)=0;

% Sum over the matrix to find the energy associated with the k,fth dipole

individualEnergy=sum(energyMatrix(:));

end

