function [ latticeHeight,latticeWidth,basisVector1,basisVector2,...
    unitCellHeight,unitCellWidth,dipoleUnitCell,dipoleLengthUnitCell] = coulombGetCrystalParameters
%getCrystalParameters prompts the user for information on a dipole crystal
%   More detailed explanation to follow?

% Initialize variables
latticeWidth=0;
latticeHeight=0;
basisVector1=[0;0];
basisVector2=[0;0];
unitCellHeight=0;
unitCellWidth=0;

% Prompt user for all inputs but dipole moment
% Note: all vectors are implicitly expressed in a 
% rectangular Cartesian system
latticeHeight=input('Please enter lattice height (in molecules): ');
latticeWidth=input('Please enter lattice width (in molecules): ');
disp('Note: the following vectors must be entered with respect to the X-Y basis.')
disp('Also note that they should not necessarily be unit vectors.')
basisVector1=input('Please enter the first basis vector, in [x;y] format (in meters): ');
basisVector2=input('Please enter the second basis vector, in [x;y] format (in meters): ');
unitCellHeight=input('Please enter the height (in molecules) of the basic structure: ');
unitCellWidth=input('Please enter the width (in molecules) of the basic structure: ');

% Prompt user for information about dipole moments of basic repeating unit
disp('Please enter all dipole moments in Debye.')
dipoleUnitCell=zeros(unitCellHeight,unitCellWidth,2);
dipoleLengthUnitCell=zeros(unitCellHeight,unitCellWidth);

dipoleQuery=input('Do all dipoles have the same dipole length? (y/n): ','s');

if strcmp(dipoleQuery,'y')
    
    dipoleLengthUnitCell(:,:)=input('Please enter the dipole length (in meters): ');
    k=0;
    for k=1:unitCellHeight
        f=0;
        for f=1:unitCellWidth
            dipoleMomentInputString=sprintf('Please input dipole moment vector for molecule in row %d, column %d, in [x;y] format: ', k, f);
            dipoleUnitCell(k,f,:)=3.335641e-30*input(dipoleMomentInputString);
        end
    end
    
elseif strcmp(dipoleQuery,'n')
    
    k=0;
    for k=1:unitCellHeight
        f=0;
        for f=1:unitCellWidth
            dipoleMomentInputString=sprintf('Please input dipole moment vector for molecule in row %d, column %d, in [x;y] format: ', k, f);
            dipoleUnitCell(k,f,:)=3.335641e-30*input(dipoleMomentInputString);
            dipoleLengthInputString=sprintf('Please enter the length of the dipole in row %d, column %d (in meters): ', k, f);
            dipoleLengthUnitCell(k,f)=input(dipoleLengthInputString);
        end
    end
    
end

end