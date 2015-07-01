theta=[0:360];
phi=[0:90];
phi=90-phi;

theta=2*pi*theta/360;
phi=2*pi*phi/360;
[P,T]=meshgrid(phi,theta);
[X,Y,Z]=sph2cart(T,P,1);
surf(X,Y,Z,squareEPerDipoleOutputMatrix);
colorbar; shading interp; daspect([1 1 1]); axis tight;