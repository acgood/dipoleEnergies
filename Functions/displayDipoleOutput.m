function displayDipoleOutput( totalEnergy,latticeHeight,latticeWidth)
%displayDipoleOutput displays the output of the calculation
%   Given a total energy and the dimensions of the crystal, this function
%   displays the total energy and energy per dipole in Joules and eV.

A=sprintf('Total energy for the lattice: %e J',totalEnergy);
B=sprintf('Total energy for the lattice: %e eV',totalEnergy/(1.602177e-19));
C=sprintf('Energy per dipole: %e J',totalEnergy/(latticeHeight*latticeWidth));
D=sprintf('Energy per dipole: %e eV',totalEnergy/(latticeHeight*latticeWidth*1.602177e-19));

disp(A)
disp(B)
disp(C)
disp(D)

end

