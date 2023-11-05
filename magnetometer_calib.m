clear all
%
% Original algorithm
%
% https://arxiv.org/ftp/arxiv/papers/1506/1506.02776.pdf
%
function [xc, yc, zc, R] = sphere_fit(x, y, z)
	[N,dum] = size(x);
	Sx = sum(x); Sy = sum(y); Sz = sum(z);
	Sxx = sum(x.*x); Syy = sum(y.*y);
	Szz = sum(z.*z); Sxy = sum(x.*y);
	Sxz = sum(x.*z); Syz = sum(y.*z);
	Sxxx = sum(x.*x.*x); Syyy = sum(y.*y.*y);
	Szzz = sum(z.*z.*z); Sxyy = sum(x.*y.*y);
	Sxzz = sum(x.*z.*z); Sxxy = sum(x.*x.*y);
	Sxxz = sum(x.*x.*z); Syyz =sum(y.*y.*z);
	Syzz = sum(y.*z.*z);
	A1 = Sxx +Syy +Szz;
	a = 2*Sx*Sx-2*N*Sxx;
	b = 2*Sx*Sy-2*N*Sxy;
	c = 2*Sx*Sz-2*N*Sxz;
	d = -N*(Sxxx +Sxyy +Sxzz)+A1*Sx;

	e = 2*Sx*Sy-2*N*Sxy;
	f = 2*Sy*Sy-2*N*Syy;
	g = 2*Sy*Sz-2*N*Syz;
	h = -N*(Sxxy +Syyy +Syzz)+A1*Sy;
	j = 2*Sx*Sz-2*N*Sxz;
	k = 2*Sy*Sz-2*N*Syz;
	l = 2*Sz*Sz-2*N*Szz;
	m = -N*(Sxxz +Syyz + Szzz)+A1*Sz;
	delta = a*(f*l - g*k)-e*(b*l-c*k) + j*(b*g-c*f);

	xc = (d*(f*l-g*k) -h*(b*l-c*k) +m*(b*g-c*f))/delta;
	yc = (a*(h*l-m*g) -e*(d*l-m*c) +j*(d*g-h*c))/delta;
	zc = (a*(f*m-h*k) -e*(b*m-d*k) +j*(b*h-d*f))/delta;
	R = sqrt(xc^2+yc^2+zc^2+(A1-2*(xc*Sx+yc*Sy+zc*Sz))/N);
endfunction

%
% Accumulating version which doesn't need to store multiple samples
%
global Sx = 0;   global Sy = 0; 	global Sz = 0
global Sxx = 0;  global Syy = 0;
global Szz = 0;  global Sxy = 0;
global Sxz = 0;  global Syz = 0;
global Sxxx = 0; global Syyy = 0;
global Szzz = 0; global Sxyy = 0;
global Sxzz = 0; global Sxxy = 0;
global Sxxz = 0; global Syyz = 0;
global Syzz = 0;

function sphere_fit_accumulate(x, y, z)
	global Sx  ;   global Sy; 	global Sz;
	global Sxx ;  global Syy;
	global Szz ;  global Sxy;
	global Sxz ;  global Syz;
	global Sxxx; global Syyy;
	global Szzz; global Sxyy;
	global Sxzz; global Sxxy;
	global Sxxz; global Syyz;
	global Syzz;
	Sx   = Sx+x;         Sy   = Sy+y; Sz = Sz+z;
	Sxx  = Sxx+x.*x;     Syy  = Syy+y.*y;
	Szz  = Szz+z.*z;     Sxy  = Sxy+x.*y;
	Sxz  = Sxz+x.*z;     Syz  = Syz+y.*z;
	Sxxx = Sxxx+x.*x.*x; Syyy = Syyy+y.*y.*y;
	Szzz = Szzz+z.*z.*z; Sxyy = Sxyy+x.*y.*y;
	Sxzz = Sxzz+x.*z.*z; Sxxy = Sxxy+x.*x.*y;
	Sxxz = Sxxz+x.*x.*z; Syyz = Syyz+y.*y.*z;
	Syzz = Syzz+y.*z.*z;
endfunction

function [xc, yc, zc, R] = sphere_fit_finalize(N)
	global Sx  ;   global Sy; 	global Sz;
	global Sxx ;  global Syy;
	global Szz ;  global Sxy;
	global Sxz ;  global Syz;
	global Sxxx; global Syyy;
	global Szzz; global Sxyy;
	global Sxzz; global Sxxy;
	global Sxxz; global Syyz;
	global Syzz;
	A1 = Sxx +Syy +Szz;
	a = 2*Sx*Sx-2*N*Sxx;
	b = 2*Sx*Sy-2*N*Sxy;
	c = 2*Sx*Sz-2*N*Sxz;
	d = -N*(Sxxx +Sxyy +Sxzz)+A1*Sx;

	e = 2*Sx*Sy-2*N*Sxy;
	f = 2*Sy*Sy-2*N*Syy;
	g = 2*Sy*Sz-2*N*Syz;
	h = -N*(Sxxy +Syyy +Syzz)+A1*Sy;
	j = 2*Sx*Sz-2*N*Sxz;
	k = 2*Sy*Sz-2*N*Syz;
	l = 2*Sz*Sz-2*N*Szz;
	m = -N*(Sxxz +Syyz + Szzz)+A1*Sz;
	delta = a*(f*l - g*k)-e*(b*l-c*k) + j*(b*g-c*f);

	xc = (d*(f*l-g*k) -h*(b*l-c*k) +m*(b*g-c*f))/delta;
	yc = (a*(h*l-m*g) -e*(d*l-m*c) +j*(d*g-h*c))/delta;
	zc = (a*(f*m-h*k) -e*(b*m-d*k) +j*(b*h-d*f))/delta;
	R = sqrt(xc^2+yc^2+zc^2+(A1-2*(xc*Sx+yc*Sy+zc*Sz))/N);
endfunction

%
% Main code
%

xyz=csvread('data.csv');

x=xyz(:,1);
y=xyz(:,2);
z=xyz(:,3);

% Original version
[xc, yc, zc, r] = sphere_fit(x,y,z)

% Accumulating version
[N,dum] = size(x);
num=0;
for i = 1:N
	num=num+1;
	sphere_fit_accumulate(x(i),y(i),z(i));
end
[xc, yc, zc, r] = sphere_fit_finalize(num)


%scatter3(x,y,z)
scatter3((x-xc)/r, (y-yc)/r, (z-zc)/r)
axis equal