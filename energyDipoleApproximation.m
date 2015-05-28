% This script is designed to calculate the interaction energy of a lattice
% of dipoles using the dipole approximation.

% Query user for lattice parameters

[ latticeHeight,latticeWidth,basisVector1,basisVector2,...
    unitCellHeight,unitCellWidth,dipoleUnitCell] = getCrystalParameters;

% Create position matrix

[ positionMatrix ] = constructPositionMatrix( latticeHeight,latticeWidth,...
    basisVector1,basisVector2);

% Create dipole moment matrix

[ dipoleMomentMatrix ] = constructDipoleMomentMatrix( latticeHeight,latticeWidth,...
    unitCellHeight,unitCellWidth,dipoleUnitCell);

% Calculate the crystal's electrostatic dipole energy

[ totalEnergy ] = calculateDipoleEnergy( latticeHeight,...
    latticeWidth,positionMatrix,dipoleMomentMatrix );

energyPerDipole=totalEnergy/(latticeHeight*latticeWidth);