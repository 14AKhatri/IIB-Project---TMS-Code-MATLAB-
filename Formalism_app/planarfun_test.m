syms x y l m ps

Y(x,y,l,m,ps) = (sin(l*(x+ps/2)*pi/ps)).*(sin(m*(y+ps/2)*pi/ps));

Y_grad(x,y,l,m,ps) = [Y(x-0.5,y,l,m,ps)-Y(x+0.5,y,l,m,ps), Y(x,y-0.5,l,m,ps)-Y(x,y+0.5,l,m,ps), 0];

Y_curl(x,y,l,m,ps) = cross(Y_grad(x,y,l,m,ps), [0,0,1]);




k=1;
l=1;

for i=-100:5:100

    for j=-100:5:100

        cur(k,l,1:3) = double(Y_curl(i,j,1,1,200));
        k=k+1;

    end
    l=l+1;
    k=1;
end

quiver(cur(:,:,1),cur(:,:,2))
axis equal

f(x,y) = Y_curl(x,y,1,1,200)+0.3*Y_grad(x,j,1,1,200)