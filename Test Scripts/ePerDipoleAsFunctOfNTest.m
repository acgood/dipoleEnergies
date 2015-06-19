% This script is designed to calculate the energy per dipole as a function
% of lattice size.

% Define lattice parameters
basisVector1=[1;0];
basisVector2=[0;1];
unitCellHeight=1;
unitCellWidth=1;
dipoleUnitCell(1,1,:)=[0;1];

energiesForSquareCrystals=zeros(200,2);

k=0;
for k=1:200
    latticeHeight=0;
    latticeWidth=0;
    
    latticeHeight=k;
    latticeWidth=k;
    
    positionMatrix=zeros(latticeHeight,latticeWidth,2);
    dipoleMomentMatrix=zeros(latticeHeight,latticeWidth,2);
    totalEnergy=0;
    
    [ positionMatrix ] = constructPositionMatrix( latticeHeight,latticeWidth,...
    basisVector1,basisVector2);

    [ dipoleMomentMatrix ] = constructDipoleMomentMatrix( latticeHeight,latticeWidth,...
    unitCellHeight,unitCellWidth,dipoleUnitCell);

    [ totalEnergy ] = calculateDipoleEnergy( latticeHeight,...
    latticeWidth,positionMatrix,dipoleMomentMatrix );

    energiesForSquareCrystals(k,1)=totalEnergy;
    energiesForSquareCrystals(k,2)=totalEnergy/(latticeHeight*latticeWidth);
end

save ePerDipoleAsFuncOfN.dat energiesForSquareCrystals -ascii