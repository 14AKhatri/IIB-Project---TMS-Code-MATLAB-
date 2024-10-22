function [my_sup,my_dir,lookup] = planarmode(gridsize,planarsize,m_max,n_max,coillift,my_rot)

d=[0; 0; 1]; %Normalenvektor zur Ebene

[x,y] = meshgrid(linspace(-planarsize/2,planarsize/2,gridsize), linspace(-planarsize/2,planarsize/2,gridsize));

W(1:gridsize,1:gridsize,1:m_max*n_max)=zeros;
z(1:gridsize,1:gridsize)=coillift;


ii=1;
for m=1:m_max
    for n=1:n_max
        
        Y(:,:)=(sin(n*(x+planarsize/2)*pi/planarsize)).*(sin(m*(y+planarsize/2)*pi/planarsize));
        [U1,V1]=gradient(Y(:,:));
        
        lookup(ii,1)= m;
        lookup(ii,2)= n;
        ii=ii+1;
        
        for i=1:gridsize
            for j=1:gridsize
                b=[U1(i,j);V1(i,j);0];
                c=cross(d,b);
                U(i,j,(m-1)*n_max+n)=c(1,:);
                V(i,j,(m-1)*n_max+n)=c(2,:);
            end
        end
        
    end
end
clear U1 V1

my_sup(:,1) = reshape(x,[],1);
my_sup(:,2) = reshape(y,[],1);
my_sup(:,3) = reshape(z,[],1);

my_sup = my_sup * my_rot;

for i=1:size(lookup,1)

    my_dir(:,1,i) = reshape(U(:,:,i),[],1);
    my_dir(:,2,i) = reshape(V(:,:,i),[],1);
    my_dir(:,3,i) = reshape(W(:,:,i),[],1);

    my_dir(:,:,i) = my_dir(:,:,i) * my_rot;

end

end