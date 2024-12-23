%% Problem 1
N= 1000; % number of dicretization points on the curve
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
vecdsdt = sqrt(rprimvec.^2+rvec.^2);
for i = 1:N
    for j = 1:N
        nu_j = [nu1(j), nu2(j)];

        r_j = [y1(j), y2(j)];
        r_i = [y1(i), y2(i)];
        
        diff = r_j - r_i;
        hankel_diff = besselh(1,k*norm(diff));
        auxillary = (1i*k/4)*hankel_diff/(norm(diff)); 

   
        
        A(i,j) = dot(nu_j, diff)*auxillary;
    end
end

%We need to fix the diagonals
for i = 1:N
    numerator = rprimvec(i)^2 - 0.5*rbisvec(i)*rvec(i)+0.5*rvec(i)^2;
    denominator = 2*pi*(rprimvec(i)^2 + rvec(i)^2)^(3/2);
    
    A(i,i) = numerator/denominator;
end



p = [0; -3];
uAn = @(x,y) -(1i/4)*besselh(0,k*norm([x;y]-p));
dx_u = @(x,y) 1i*k/(4*norm([x;y]-p)) * besselh(1,k*norm([x;y]-p))*(x-p(1));
dy_u = @(x,y) 1i*k/(4*norm([x;y]-p)) * besselh(1,k*norm([x;y]-p))*(x-p(2));

g1 = zeros(1,N);
g2 = zeros(1,N);
for i = 1:N
    g1(i) = dx_u(y1(i), y2(i))*nu1(i);
    g2(i) = dy_u(y1(i), y2(i))*nu2(i);
end
gvec = g1'+g2';

kmat = (-eye(N)/2+ 2*pi/N* A *diag(vecdsdt));

hvec = kmat \ gvec;
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
            phivec= (-1i / 4) * besselh(0, k * vecnorm([y1;y2] - [x1;x2])); %kernel expression in terms of xi, yi, nui, i=1,2
            vfield(ix1,ix2) = (phivec *(hvec.* vecdsdt).')*2*pi/N; %trapezoidal rule
            uAnField(ix1,ix2) = real(uAn(x1,x2));
        end
    end
end
ufield = real(vfield);
%% 
imagesc(x1field, x2field, ufield.')
axis xy
colormap turbo
pbaspect([1 1 1])
colorbar
title('Analytical solution of u(x,y)')



