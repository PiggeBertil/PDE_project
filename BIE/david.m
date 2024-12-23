%% Problem 1
N= 500; % number of dicretization points on the curve
tvec= linspace(-pi + 2*pi/N, pi, N);

rvec= 3+cos(4.* tvec + pi); % analytic expression for r(t)
rprimvec= -4*sin(4.*tvec+pi); % analytic expression for r'(t)
rbisvec= -16*cos(4.*tvec+pi); % analytic expression for r''(t)

y1 = rvec .* cos(tvec);
y2 = rvec .* sin(tvec);

%Compute normal vectors 
nu1= rvec .* cos(tvec) + rprimvec .* sin(tvec);
nu2= rvec .* sin(tvec) - rprimvec .* cos(tvec);
nu1= nu1 ./ sqrt( rvec.^2+ rprimvec.^2 );
nu2= nu2 ./ sqrt( rvec.^2+ rprimvec.^2 );


%We create the A matrix
k = 1;
omega = k;
A = zeros(N);
dsdtvec = sqrt(rprimvec.^2+rvec.^2);
for i = 1:N
    for j = 1:N
        nu_i = [nu1(i), nu2(i)];
    
        r_j = [y1(j), y2(j)];
        r_i = [y1(i), y2(i)];
        
        diff = r_i - r_j;
        hankel_diff = besselh(1,1,k*norm(diff));
        auxillary = (1i*k/4)*hankel_diff/(norm(diff)); 
    
        A(i,j) = dot(nu_i, diff)*auxillary;
    end
end

%We need to fix the diagonals
for i = 1:N
numerator = rprimvec(i)^2 - 0.5*rbisvec(i)*rvec(i)+0.5*rvec(i)^2;
denominator = 2*pi*(rprimvec(i)^2 + rvec(i)^2)^(3/2);
dsdt = dsdtvec(i);
A(i,i) = numerator/denominator;
end

subplot(1,3,1)
imagesc(tvec,tvec,real(A).')
axis xy
colorbar
title('re kmat')

subplot(1,3,2)
imagesc(tvec,tvec,imag(A).')
axis xy
colorbar
title('im kmat')

subplot(1,3,3)
imagesc(tvec,tvec,abs(A).')
axis xy
colorbar
title('abs kmat')

%% Problem 2

p = [0; -3];
uAn = @(x,y) -(1i/4)*besselh(0,1,k*norm([x;y]-p));

%{
dx_u = @(x,y) 1i*k./(4*vecnorm([x;y]-p)) .* besselh(1,1,k*vecnorm([x;y]-p)).*(x-p(1));
dy_u = @(x,y) 1i*k./(4*vecnorm([x;y]-p)) .* besselh(1,1,k*vecnorm([x;y]-p)).*(y-p(2));

g1 = dx_u(y1,y2).*nu1;
g2 = dy_u(y1,y2).*nu2;

gvec = g1.'+g2.';
%}

gvec = 1i*k/4 * besselh(1,1,k*vecnorm([y1;y2]-p))./vecnorm([y1;y2]-p) .* (dot(([y1;y2]-p),[nu1;nu2]));
kmat = (-eye(N)/2+ 2*pi/N* A *diag(dsdtvec));

hvec = kmat \ gvec.';
hvec = hvec.';

M= 300; % pick something here
x1field = linspace(-4, 4, M);
x2field = linspace(-4, 4, M);
vfield = zeros(M,M);
uAnField = zeros(M,M);


for ix1=1:M
for ix2=1:M
    x1=x1field(ix1);
    x2=x2field(ix2);
    t=angle(complex(x1,x2));
    radius= 3 + cos(4 * t + pi); %formula for rvec
    if x1^2+ x2^2 < radius^2
        phivec= (-1i / 4) * besselh(0,1, k * vecnorm([y1;y2] - [x1;x2])); %kernel expression in terms of xi, yi, nui, i=1,2
        vfield(ix1,ix2) = (phivec *(hvec.* dsdtvec).')*2*pi/N; %trapezoidal rule
        uAnField(ix1,ix2) = uAn(x1,x2);
    end
end
end

t = 0;
ufield_real = real(vfield*exp(-1i*omega*t));
ufield_im = imag(vfield*exp(-1i*omega*t));

uAnField_real = real(uAnField);
uAnField_im = imag(uAnField);

subplot(2,3,1)
imagesc(x1field, x2field, ufield_real.')
axis xy
colormap turbo
pbaspect([1 1 1])
colorbar
title('Numerical real solution')

subplot(2,3,4)
imagesc(x1field, x2field, ufield_im.')
axis xy
colormap turbo
pbaspect([1 1 1])
colorbar
title('Numerical im solution')

subplot(2,3,2)
imagesc(x1field, x2field, uAnField_real.')
axis xy
colormap turbo
pbaspect([1 1 1])
colorbar
title('Analytic solution real')


subplot(2,3,5)
imagesc(x1field, x2field, uAnField_im.')
axis xy
colormap turbo
pbaspect([1 1 1])
colorbar
title('Analytic solution im')
%% Problem 2 part 2
%Find the average in a square inside
average_numerical_real = mean(mean(ufield_real(M/3:2*M/3, M/3:2*M/3)));
average_analytic_real = mean(mean(uAnField_real(M/3:2*M/3, M/3:2*M/3)));

average_numerical_im = mean(mean(ufield_im(M/3:2*M/3, M/3:2*M/3)));
average_analytic_im = mean(mean(uAnField_im(M/3:2*M/3, M/3:2*M/3)));

subplot(2,3,3)
imagesc(x1field, x2field, log10(abs(ufield_real-average_numerical_real+average_analytic_real-uAnField_real)).')
axis xy
colormap turbo
pbaspect([1 1 1])
colorbar
title('Log10 error')

subplot(2,3,6)
imagesc(x1field, x2field, log10(abs(ufield_im-average_numerical_im+average_analytic_im-uAnField_im)).')
axis xy
colormap turbo
pbaspect([1 1 1])
colorbar