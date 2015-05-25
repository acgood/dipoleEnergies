% This script is designed to calculate the interaction energy of a lattice
% of dipoles using the dipole approximation.

% Initialize all variables and arrays
totalEnergy=0;
energyPerDipole=0;
latticeWidth=0;
latticeHeight=0;
basisVector1=[0;0];
basisVector2=[0;0];

% Input physical constants
epsilon=8.85418782e-12;

% Query user for lattice parameters
% Note: the basis vectors and position matrix are implicitly expressed in a 
% rectangular Cartesian system
latticeHeightInitial=input('Please enter lattice height (in molecules): ');
latticeWidthInitial=input('Please enter lattice width (in molecules): ');
disp('Note: the following vectors must be entered with respect to the X-Y basis.')
disp('Also note that they should not necessarily be unit vectors.')
basisVector1=input('Please enter the first basis vector, in [x;y] format (in meters): ');
basisVector2=input('Please enter the second basis vector, in [x;y] format (in meters): ');

% Create dipole moment matrix
% Query user for structure of basic repeating unit
unitCellHeight=input('Please enter the height (in molecules) of the basic structure: ');
unitCellWidth=input('Please enter the width (in molecules) of the basic structure: ');
k=0;
for k=1:unitCellHeight
	f=0;
	for f=1:unitCellWidth
		dipoleMomentInputString=sprintf('Please input dipole moment vector for molecule in row %d, column %d, in [x;y] format: ', k, f);
		dipoleUnitCell(k,f,:)=input(dipoleMomentInputString);
	end
end
% Compare dimensions of lattice to dimensions of unit cell; expand lattice if 
% necessary to make dimensions an integer multiple of unit cell dimensions
ratioHeight=0;
ratioWidth=0;
ratioHeight=latticeHeightInitial/unitCellHeight;
ratioWidth=latticeWidthInitial/unitCellWidth;
isHeightChanged=0;
isWidthChanged=0;
if mod(ratioHeight,1)~=0
	latticeHeight=ceil(ratioHeight)*unitCellHeight;
	isHeightChanged=1;
else
	latticeHeight=latticeHeightInitial;
end

if mod(ratioHeight,1)~=0
	latticeWidth=ceil(ratioWidth)*unitCellWidth;
	isWidthChanged=1;
else
	latticeWidth=latticeWidthInitial;
end

if isHeightChanged==1 | isWidthChanged==1
	changeMessage=sprintf('Crystal dimensions have been resized to %d rows by %d columns.',latticeHeight,latticeWidth)
end
% Create dipole moment matrix for crystal by tiling the unit cell
dipoleMomentMatrix=zeros(latticeHeight,latticeWidth,2);
dipoleMomentMatrix=repmat(dipoleUnitCell,ceil(ratioHeight),ceil(ratioWidth));


% Create position matrix
positionMatrix=zeros(latticeHeight,latticeWidth,2);
k=0;
for k=1:latticeHeight
	f=0;
	for f=1:latticeWidth
		positionMatrix(k,f,1)=((f-1)*basisVector1(1,1) + (k-1)*basisVector2(1,1));
		positionMatrix(k,f,2)=((f-1)*basisVector1(2,1) + (k-1)*basisVector2(2,1));
	end
end

% Sum up the energy contribution of each interaction term
% Note: term1, term2, etc refer to the terms in the dipole energy equation
energyRaw=0;
k=0;
for k=1:latticeHeight
	f=0;
	for f=1:latticeWidth
		g=0;
		for g=1:latticeHeight
			h=0;
			for h=1:latticeWidth

				if k==g && f==h

					increment=0;

				else
					
					increment=0;
					p1=[0;0];
					p2=[0;0];
					r1=[0;0];
					r2=[0;0];
					rVector=[0;0];
					rScalar=0;
					rVectorNormalized=[0;0];

					p1(1,1)=dipoleMomentMatrix(k,f,1);
					p1(2,1)=dipoleMomentMatrix(k,f,2);
					p2(1,1)=dipoleMomentMatrix(g,h,1);
					p2(2,1)=dipoleMomentMatrix(g,h,2);
					r1(1,1)=positionMatrix(k,f,1);
					r1(2,1)=positionMatrix(k,f,2);
					r2(1,1)=positionMatrix(g,h,1);
					r2(2,1)=positionMatrix(g,h,2);

					rVector=r1-r2;
					rScalar=norm(rVector);
					rVectorNormalized=rVector/rScalar;

					term1=0;
					term2=0;
					term3=0;
					term1=p1(1,1)*p2(1,1)+p1(2,1)*p2(2,1);
					term2=p1(1,1)*rVectorNormalized(1,1)+p1(2,1)*rVectorNormalized(2,1);
					term3=p2(1,1)*rVectorNormalized(1,1)+p2(2,1)*rVectorNormalized(2,1);
					increment=(term1-3*term2*term3)/(rScalar^3);
				end

				energyRaw=energyRaw+increment;
			end
		end
	end
end

% Multiply by constants and calculate relevant quantities

totalEnergy=energyRaw/(8*pi*epsilon);
energyPerDipole=totalEnergy/(latticeWidth*latticeHeight);



