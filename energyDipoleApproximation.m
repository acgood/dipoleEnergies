% This script is designed to calculate the interaction energy of a lattice
% of dipoles using the dipole approximation.

% Initialize all variables and arrays
TotalEnergy=0;
EnergyPerDipole=0;
LattWidth=0;
LattHeight=0;
basis1=[0;0];
basis2=[0;0];

% Input physical constants
epsilon=8.85418782e-12;

% Input lattice and dipole parameters
% Note: the basis vectors and position matrix are implicitly expressed in a 
% rectangular Cartesian system
LattHeight=1;
LattWidth=1;
basis1=[1;0];
basis2=[0;1];

% Create dipole moment matrix
dipm=zeros(LattHeight,LattWidth,2);
k=0;
for k=1:LattHeight
	f=0;
	for f=1:LattWidth
		dipm(k,f,1)=0;
		dipm(k,f,2)=1;
	end
end

% Create position matrix
posn=zeros(LattHeight,LattWidth,2);
k=0;
for k=1:LattHeight
	f=0;
	for f=1:LattWidth
		posn(k,f,1)=((f-1)*basis1(1,1) + (k-1)*basis2(1,1));
		posn(k,f,2)=((f-1)*basis1(2,1) + (k-1)*basis2(2,1));
	end
end

% Sum up the energy contribution of each interaction term
% Note: Creating an energy matrix expressing the interaction energy of each
% dipole with every other would be too memory intensive
eraw=0;
k=0;
for k=1:LattHeight
	f=0;
	for f=1:LattWidth
		g=0;
		for g=1:LattHeight
			h=0;
			for h=1:LattWidth

				if k==g && f==h

					inc=0;

				else
					
					inc=0;
					p1=[0;0];
					p2=[0;0];
					r1=[0;0];
					r2=[0;0];
					rvector=[0;0];
					rscalar=0;
					rvecnorm=[0;0];

					p1(1,1)=dipm(k,f,1);
					p1(2,1)=dipm(k,f,2);
					p2(1,1)=dipm(g,h,1);
					p2(2,1)=dipm(g,h,2);
					r1(1,1)=posn(k,f,1);
					r1(2,1)=posn(k,f,2);
					r2(1,1)=posn(g,h,1);
					r2(2,1)=posn(g,h,2);

					rvector=r1-r2;
					rscalar=norm(rvector);
					rvecnorm=rvector/rscalar;

					term1=0;
					term2=0;
					term3=0;
					term1=p1(1,1)*p2(1,1)+p1(2,1)*p2(2,1);
					term2=p1(1,1)*rvecnorm(1,1)+p1(2,1)*rvecnorm(2,1);
					term3=p2(1,1)*rvecnorm(1,1)+p2(2,1)*rvecnorm(2,1);
					inc=(term1-3*term2*term3)/(rscalar^3);
				end

				eraw=eraw+inc;
			end
		end
	end
end

% Multiply by constants and calculate relevant quantities

TotalEnergy=eraw/(8*pi*epsilon);
EnergyPerDipole=TotalEnergy/(LattWidth*LattHeight);



