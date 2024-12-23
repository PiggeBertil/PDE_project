N=300;
tvec = linspace(-pi+2*pi/N, pi, N);
rvec = 3+cos(4.*tvec+pi);
rprimvec = -4*sin(4.*tvec+pi);
rbisvec = -16*cos(4.*tvec+pi);
y1 = rvec .*cos(tvec);
y2 = rvec .*sin(tvec);
nu1 = rvec .* cos(tvec) + rprimvec .* sin(tvec);
nu2 = rvec .* sin(tvec) - rprimvec .* cos(tvec);
nu1 = nu1 ./ sqrt( rvec.^2+ rprimvec.^2 );
nu2 = nu2 ./ sqrt( rvec.^2+ rprimvec.^2 );

A_k = zeros(N,N);
vecdsdt = sqrt(rprimvec.^2+rvec.^2);
p = [0,4];
nu = [nu1; nu2];
k=1;
for i = 1:N
    for j = 1:N
        nu_i = [nu1(i), nu2(i)];
    
        r_j = [y1(j), y2(j)];
        r_i = [y1(i), y2(i)];
        
        diff = r_i - r_j;
        hankel_diff = besselh(1,1,k*norm(diff));
        auxillary = (1i*k/4)*hankel_diff/(norm(diff)); 
    
        A_k(i,j) = dot(nu_i, diff)*auxillary;
    end
end

for i=1:N
    taljare = rprimvec(i)^2-0.5*rbisvec(i)*rvec(i) + 0.5*rvec(i)^2;
    namnare = 2*pi*(rprimvec(i)^2+rvec(i)^2)^(3/2);
    A_k(i,i) = taljare/namnare;    
end

%%
M = 300; % pick something here
x1field = linspace(-4, 4, M);
x2field = linspace(-4, 4, M);
ufield = zeros(M,M);
exactfield = zeros(M,M);

gvec = (nu1.*y1+nu2.*(y2-4))./(sqrt((y1).^2+(y2-4).^2).*(1i*k/4*besselh(1,k.*sqrt(y1.^2+(y2-4).^2))));
ineverse_matrix = (-eye(N)/2+ 2*pi/N * A_k * diag(vecdsdt));
hvec = ineverse_matrix\gvec.';
hvec = hvec.';

for ix1=1:M
    for ix2=1:M
        x1=x1field(ix1);
        x2=x2field(ix2);
        t=angle(complex(x1,x2));
        radius= 3+cos(4*t+pi); %formula for rvec
        if x1^2+ x2^2< radius^2
           phivec = -1i/4 * besselh(0,k*sqrt((x1-y1).^2 + (x2-(y2)).^2)); %kernel expression in terms of xi, yi, nui, i=1,2
           exactfield(ix1,ix2) = real(1i/4*besselh(0,k*sqrt(x1.^2+(x2-4).^2)));
           ufield(ix1,ix2) = (phivec*((hvec.* vecdsdt).')*2*pi/N)-4;
        end
    end
end

ufield = real(ufield);
%% 
imagesc(x1field, x2field, exactfield.')
axis xy
colormap turbo
pbaspect([1 1 1])
colorbar
title('Analytical solution of u(x,y)')

%%
errorfield = log10(abs(ufield-exactfield));
imagesc(x1field, x2field, ufield.')
k=1;
axis xy
colormap turbo
pbaspect([1 1 1])
colorbar
title('Numerical solution of u(x,y)')

%%
imagesc(x1field, x2field, errorfield.')

axis xy
colormap turbo
pbaspect([1 1 1])
colorbar
title('Error plot')
