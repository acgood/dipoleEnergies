function [ positionMatrix ] = constructPositionMatrix( latticeHeight,latticeWidth,...
    basisVector1,basisVector2)
%constructPositionMatrix does as the name implies
%   The function creates a matrix of the positions in physical
%   space of crystal molecules given a set of input parameters

positionMatrix=zeros(latticeHeight,latticeWidth,2);
k=0;
for k=1:latticeHeight
	f=0;
	for f=1:latticeWidth
		positionMatrix(k,f,1)=((f-1)*basisVector1(1,1) + (k-1)*basisVector2(1,1));
		positionMatrix(k,f,2)=((f-1)*basisVector1(2,1) + (k-1)*basisVector2(2,1));
	end
end


end

