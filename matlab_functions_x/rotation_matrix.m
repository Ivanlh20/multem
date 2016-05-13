function[Rm]=rotation_matrix(u0, theta)
u0 = u0/norm(u0);
theta = theta*pi/180;
alpha = 1-cos(theta);
beta = sin(theta);

Rm = zeros(3,3);
Rm(1,1) = 1.0 + alpha*(u0(1)*u0(1)-1);
Rm(2,1) = u0(2)*u0(1)*alpha + u0(3)*beta;
Rm(3,1) = u0(3)*u0(1)*alpha - u0(2)*beta;

Rm(1,2) = u0(1)*u0(2)*alpha - u0(3)*beta;
Rm(2,2) = 1.0 + alpha*(u0(2)*u0(2)-1);
Rm(3,2) = u0(3)*u0(2)*alpha + u0(1)*beta;

Rm(1,3) = u0(1)*u0(3)*alpha + u0(2)*beta;
Rm(2,3) = u0(2)*u0(3)*alpha - u0(1)*beta;
Rm(3,3) = 1.0 + alpha*(u0(3)*u0(3)-1);