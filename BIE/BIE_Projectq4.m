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
p = [0;4];
%

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
           phivec = (-1i / 4) * besselh(0, k * vecnorm([y1;y2] - [x1;x2])); %kernel expression in terms of xi, yi, nui, i=1,2
           exactfield(ix1,ix2) = real(1i/4*besselh(0,k*sqrt(x1.^2+(x2-4).^2)));
           ufield(ix1,ix2) = (phivec*((hvec.* vecdsdt).')*2*pi/N);
        end
    end
end

ufield = real(ufield);
errorfield = log10(abs(ufield-exactfield));

%% 
imagesc(x1field, x2field, exactfield.')
axis xy
colormap turbo
pbaspect([1 1 1])
colorbar
title('Analytical solution of u(x,y)')

%%
imagesc(x1field, x2field, ufield.')
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
