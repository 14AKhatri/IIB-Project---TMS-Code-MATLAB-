function [my_sup,my_dir,lookup,lookupinv,vecscale] = sphericalmode(lmax, coilradius, n)
%% function to generate spherical harmonic modes 
% inputs: lmax - max. degree of spherical harmonics
% coilradius - radius of sphere on which harmonics will be calculated
% n - depth parameter for triangulating sphere - used in generating mesh 

[mesh.Vertices, mesh.Faces] = spheretribydepth(n); 
% spherical mesh with "n" defining depth of triangulation
% mesh.Vertices - sotres corrdinates of mesh vertices
% mesh.Faces - stores triangular face indices

mesh.centerpoints=(mesh.Vertices(mesh.Faces(:,1),:) + mesh.Vertices(mesh.Faces(:,2),:) + mesh.Vertices(mesh.Faces(:,3),:))./3; 
%compute centre point of each triangle
    % ..average of triangle's vertices


xs=coilradius*mesh.centerpoints(:,1)/max(mesh.centerpoints(:,1));
ys=coilradius*mesh.centerpoints(:,2)/max(mesh.centerpoints(:,2));
zs=coilradius*mesh.centerpoints(:,3)/max(mesh.centerpoints(:,3));
%rescale centrepoints so they are on sphere 


nt=size(mesh.Faces,1); %returns length of first dimension
edges=reshape(mesh.Vertices(reshape(mesh.Faces(:,2:3)',nt*2,1),:)',3,2,nt)-repmat(reshape(mesh.Vertices(mesh.Faces(:,1),:)',3,1,nt),[1,2,1]);
mesh.size=squeeze(1/2*sqrt(sum(cross(edges(:,1,:),edges(:,2,:),1).^2,1)));
% mesh.size - area of each triangular face on mesh 
    %found using cross-product

clear edges nt

%% create lookup table - maps each mode combination to index i, to be used for spherical harmonics
%%lookupinv - does the reverse mapping to retrieve index from modes
i=1;
for mode1=1:lmax*2+1 %M
    for mode2=abs(-lmax-1+mode1)+(1)*((mode1-lmax-1)==0):lmax % def. of harmonics - M upto +- L

        lookup(i,1)=mode1;
        lookup(i,2)=mode2;
        
        lookupinv(mode1,mode2)=i;
        
        i=i+1;
        
    end
end
%lookup - two col. matrix - [mode index for L degree of spherical harmonics, index for M order of sph. harm. ]
%   the table maps index i to each valid pair of L & M upto lmax. - easy
%   way to look up L & M based on i

%lookupinv - 2D matrix - Row = mode for L values; Col = mode for M values
%   easy way to retrieve i based on specific values of L & M

vecscale=(mesh.size/max(mesh.size));
%scaling factor based on ratio between each triangle, and max area in mesh

%% calculate real spherical harmonics
for L = 1:lmax

    for M = -L:L
        
        f = curl(spherefun.sphharm(L,M)); %the curl generates tangent to the surface harmonics
        mytempsphh = feval(f,mesh.centerpoints(:,1),mesh.centerpoints(:,2),mesh.centerpoints(:,3))';

        Mxs(:,lookupinv(lmax+1+M,L)) = mytempsphh(:,1);
        Mys(:,lookupinv(lmax+1+M,L)) = mytempsphh(:,2);
        Mzs(:,lookupinv(lmax+1+M,L)) = mytempsphh(:,3);
        
    end
    
end
% for each degree L & M (upto lmax), the spherical harmonics calculated at
% each triangle's center point (spherefun.sphharm)
% curl applied to create vector field - evaluated at each center point
% stored in matrices Mxs "" - each row corresponds to centre point& column
% ... corresponds to specific harmonic mode

clear mytempsphh

%% format for processing

my_sup = [ys xs zs]; %rescaled corrdinates of center points of each tria.

my_dir(:,1,:) = Mys;
my_dir(:,2,:) = Mxs;
my_dir(:,3,:) = Mzs;
