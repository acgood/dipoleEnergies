latticeWidth=10;
latticeHeight=10;
basisVector1=0.5e-9*[1;0];
basisVector2=0.5e-9*[0;1];
unitCellHeight=1;
unitCellWidth=1;
azimuthalResolution=361;
polarResolution=91;

squareEPerDipoleOutputMatrix=zeros(azimuthalResolution,polarResolution);

k=0;
for k=1:azimuthalResolution
    f=0;
    for f=1:polarResolution
        
        dipoleUnitCell=zeros(unitCellHeight,unitCellWidth,3);
        totalEnergy=0;
        
        dipoleUnitCell(1,1,:)=[sin((f-1)*2*pi/360)*cos((k-1)*2*pi/360);sin((f-1)*2*pi/360)*sin((k-1)*2*pi/360);cos((f-1)*2*pi/360)];
        
        [ positionMatrix ] = constructPositionMatrix( latticeHeight,latticeWidth,...
            basisVector1,basisVector2);
        
        [ dipoleMomentMatrix ] = constructDipoleMomentMatrix( latticeHeight,latticeWidth,...
            unitCellHeight,unitCellWidth,dipoleUnitCell);
        
        [ totalEnergy ] = calculateDipoleEnergy( latticeHeight,...
            latticeWidth,positionMatrix,dipoleMomentMatrix );
        
        squareEPerDipoleOutputMatrix(k,f)=totalEnergy/(latticeWidth*latticeHeight);
    end
end

save SquareEPerDipole3DDipoleMoments.dat squareEPerDipoleOutputMatrix -ascii        