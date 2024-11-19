
f = 0;
my_max = 10;

for L = 1:my_max
    for M = -L:L
        
        f = f + CSUMM(lookupinv(l_max+1+M,L),my_max) * spherefun.sphharm(L,M);
        
    end
end

my_curl = curl(f);
my_grad = grad(f);

%%
clear coil
my_k = -0.002;
% my_mass = -0.45;
my_mass = -0.35;
my_dir(2,1:3) = zeros;

start = [0.66 0 0.75];
% start = [0	-30.3521463954584	85.5888994896230];
% start = [0	-25.3521463954584	87.5888994896230];
start = start * 1/norm(start);
i=1;
coil(i,:) = start;

my_mul = 10/coilradius; % sets resolution of the interpolation

for i=2:20000
    
    my_del_curl = feval(my_curl,coil(i-1,1),coil(i-1,2),coil(i-1,3));
    my_del_grad = feval(my_grad,coil(i-1,1),coil(i-1,2),coil(i-1,3));
    my_del_grad = (my_del_grad/norm(my_del_grad))*my_k;
    
    my_dir(1,:) = (my_del_curl+my_del_grad)'/(coilradius*.1);

    coil(i,:) = coil(i-1,:) + my_dir(1,:) + my_dir(2,:)*my_mass;
    
    coil(i,:) = coil(i,:) * 1/norm(coil(i,:));

    my_dir(2,:) = my_dir(1,:);
    
end

% coil = coil.*coilradius;

%% inductance

L = InductanceInt(coil.*1e-3.*coilradius)

%%

figure

% s.EdgeColor = 'none';
% s.FaceColor = [1,1,1];
% surf(xs,ys,zs,s);
quiver(my_curl)
hold on
% contour(f,-.5:0.02:.5), axis off
line(coil(:,1),coil(:,2),coil(:,3))
hold off
axis tight
axis equal
set(gca,'visible','off')


%%

figure

% s.EdgeColor = 'none';
% s.FaceColor = [1,1,1];
% surf(xs,ys,zs,s);
quiver(my_curl)
hold on
% contour(f,-.5:0.02:.5), axis off
line(coil(:,1),coil(:,2),coil(:,3))
hold off
axis tight
axis equal
set(gca,'visible','off')
%% CLOSE THE COIL

addit(:,1)=linspace(coil(1,1),coil(end,1),1000);
addit(:,2)=linspace(coil(1,2),coil(end,2),1000);
addit(:,3)=linspace(coil(1,3),coil(end,3),1000);

coil2=[coil; addit];

%% inductance

L = InductanceInt(coil2.*1e-3)

%% COIL SHORTENING

coil3=coil;

coil3(11501:end,:)=[];

coil4=resample(coil3,2,10);

[th,ph,Rr]=cart2sph(coil4(:,1),coil4(:,2),coil4(:,3));
Rr(:)=90;
[coil4(:,1),coil4(:,2),coil4(:,3)]=sph2cart(th,ph,Rr);


coil4=coil4.*1e-3;
L=InductanceInt(coil4)
coil4=coil4.*1e+3;

%%
figure

[xss,yss,zss]=sphere;
xss=xss*coilradius;
yss=yss*coilradius;
zss=zss*coilradius;

s.EdgeColor = 'none';
s.FaceColor = [1,1,1];
% surf(xss,yss,zss,s);
mesh(xss,yss,zss,'EdgeColor','k','FaceAlpha',1,'EdgeAlpha',.1)
% hold on
line(-coil4(:,1),coil4(:,2),coil4(:,3))
% plot3(-start(1),start(2),start(3),'+','Color','r','MarkerSize',10);
% hold off
axis tight
axis equal
set(gca,'visible','off')