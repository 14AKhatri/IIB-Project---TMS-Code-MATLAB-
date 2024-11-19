function [my_sup,my_dir,lookup,lookupinv,vecscale] = sphericalmode(lmax, coilradius, n)

[mesh.Vertices, mesh.Faces] = spheretribydepth(n);

mesh.centerpoints=(mesh.Vertices(mesh.Faces(:,1),:) + mesh.Vertices(mesh.Faces(:,2),:) + mesh.Vertices(mesh.Faces(:,3),:))./3;

xs=coilradius*mesh.centerpoints(:,1)/max(mesh.centerpoints(:,1));
ys=coilradius*mesh.centerpoints(:,2)/max(mesh.centerpoints(:,2));
zs=coilradius*mesh.centerpoints(:,3)/max(mesh.centerpoints(:,3));


nt=size(mesh.Faces,1);
edges=reshape(mesh.Vertices(reshape(mesh.Faces(:,2:3)',nt*2,1),:)',3,2,nt)-repmat(reshape(mesh.Vertices(mesh.Faces(:,1),:)',3,1,nt),[1,2,1]);
mesh.size=squeeze(1/2*sqrt(sum(cross(edges(:,1,:),edges(:,2,:),1).^2,1)));

clear edges nt


i=1;
for mode1=1:lmax*2+1
    for mode2=abs(-lmax-1+mode1)+(1)*((mode1-lmax-1)==0):lmax

        lookup(i,1)=mode1;
        lookup(i,2)=mode2;
        
        lookupinv(mode1,mode2)=i;
        
        i=i+1;
        
    end
end

vecscale=(mesh.size/max(mesh.size));

%% calculate real spherical harmonics
for L = 1:lmax

    for M = -L:L
        
        f = curl(spherefun.sphharm(L,M));
        mytempsphh = feval(f,mesh.centerpoints(:,1),mesh.centerpoints(:,2),mesh.centerpoints(:,3))';

        Mxs(:,lookupinv(lmax+1+M,L)) = mytempsphh(:,1);
        Mys(:,lookupinv(lmax+1+M,L)) = mytempsphh(:,2);
        Mzs(:,lookupinv(lmax+1+M,L)) = mytempsphh(:,3);
        
    end
    
end

clear mytempsphh

%% format for processing

my_sup = [ys xs zs];

my_dir(:,1,:) = Mys;
my_dir(:,2,:) = Mxs;
my_dir(:,3,:) = Mzs;
