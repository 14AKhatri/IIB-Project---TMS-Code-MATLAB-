
my_point1=[ 100  100  000];
my_point2=[-100  100  000];
my_point3=[-100 -100  000];
my_point4=[ 100 -100  000];


[xs,ys,zs] = sphere(50);

xs = xs*coilrad(3);
ys = ys*coilrad(3);
zs = zs*coilrad(3);

%%

eul=[pi/8 0 pi/4];

coilrad = [0 0 105];

my_rot = eul2rotm(eul,'XYZ');

corr_vec = [0 0 1];

corr_vec_trans = corr_vec * my_rot;

my_coillift = coilrad * my_rot;

rot_p1 = my_point1 * my_rot + my_coillift;
rot_p2 = my_point2 * my_rot + my_coillift;
rot_p3 = my_point3 * my_rot + my_coillift;
rot_p4 = my_point4 * my_rot + my_coillift;

%%

figure
surf(xs,ys,zs)
axis equal
hold on

patch([rot_p1(1) rot_p2(1) rot_p3(1) rot_p4(1)],[rot_p1(2) rot_p2(2) rot_p3(2) rot_p4(2)],[rot_p1(3) rot_p2(3) rot_p3(3) rot_p4(3)],'blue','FaceAlpha',0.2)
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

%%

hotspot = [2 4 1];

hotspot = hotspot./norm(hotspot);

el_ang = acos(dot(corr_vec,hotspot));

R_1 = eul2rotm([el_ang 0 0],'XYZ');

corr_vec_R1 =corr_vec * R_1;

az_ang = acos(dot(corr_vec_R1,hotspot));

R_2 = eul2rotm([0 0 az_ang],'XYZ');

corr_vec_R2 = corr_vec_R1 * R_2;

new_vec = corr_vec * R_1 * R_2;

%% get rotation from 2 vectors

% two random 3D vectors
% p0 = randi(10,3,1)
% p1 = randi(10,3,1)

p0 = corr_vec';
p1 = hotspot';
% calculate cross and dot products
C = cross(p0, p1) ; 
D = dot(p0, p1) ;
NP0 = norm(p0) ; % used for scaling
if ~all(C==0) % check for colinearity    
    Z = [0 -C(3) C(2); C(3) 0 -C(1); -C(2) C(1) 0] ; 
    R = (eye(3) + Z + Z^2 * (1-D)/(norm(C)^2)) / NP0^2 ; % rotation matrix
else
    R = sign(D) * (norm(p1) / NP0) ; % orientation and scaling
end
% R is the rotation matrix from p0 to p1, so that (except for round-off errors) ...
R * p0      % ... equals p1 
inv(R) * p1 % ... equals p0
