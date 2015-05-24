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
latticeHeight=input('Please enter lattice height (in molecules): ');
latticeWidth=input('Please enter lattice width (in molecules): ');
disp('Note: the following vectors must be entered with respect to the X-Y basis.')
basisVector1=input('Please enter the first basis vector, in [x;y] format: ');
basisVector2=input('Please enter the second basis vector, in [x;y] format: ');

% Create dipole moment matrix
dipoleMomentMatrix=zeros(latticeHeight,latticeWidth,2);
k=0;
for k=1:latticeHeight
	f=0;
	for f=1:latticeWidth
		dipoleMomentMatrix(k,f,1)=0;
		dipoleMomentMatrix(k,f,2)=1;
	end
end

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



