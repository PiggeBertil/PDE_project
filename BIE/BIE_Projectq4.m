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
k = 1;

for i = 1:N
    for j = 1:N
        nu_i = [nu1(i), nu2(i)];
    
        r_j = [y1(j), y2(j)];
        r_i = [y1(i), y2(i)];
        
        difference = r_i - r_j;
        hankel_func = besselh(1,1,k*norm(difference));
        taljare = (1i*k/4)*hankel_func;
        namnare = (norm(difference)); 
    
        A_k(i,j) = dot(nu_i, difference)*taljare/namnare;
    end
end

for i = 1:N
    taljare = rprimvec(i)^2 - 0.5*rbisvec(i)*rvec(i)+0.5*rvec(i)^2;
    namnare = 2*pi*(rprimvec(i)^2 + rvec(i)^2)^(3/2);
    A_k(i,i) = taljare/namnare;
end

%%
M = 300; % pick something here
x1field = linspace(-4, 4, M);
x2field = linspace(-4, 4, M);
ufield = zeros(M,M);
exactfield = zeros(M,M);

gvec = (1i*k/4) * besselh(1,1,k*vecnorm([y1;y2]-p))./vecnorm([y1;y2]-p) .* (dot(([y1;y2]-p),[nu1;nu2]));
%gvec = (nu1.*y1+nu2.*(y2-4))./(sqrt((y1).^2+(y2-4).^2).*(1i*k/4*besselh(1,k.*sqrt(y1.^2+(y2-4).^2))));
hvec = (-eye(N)/2+ 2*pi/N * A_k * diag(vecdsdt))\gvec.';
hvec = hvec.';

for ix1=1:M
   for ix2=1:M
        x1=x1field(ix1);
        x2=x2field(ix2);
        t=angle(complex(x1,x2));
        radius= 3+cos(4*t+pi); %formula for rvec
        if x1^2+ x2^2< radius^2
           phivec = (-1i / 4) * besselh(0, k * vecnorm([y1;y2] - [x1;x2])); %kernel expression in terms of xi, yi, nui, i=1,2
           exactfield(ix1,ix2) = real(-1i/4*besselh(0,k*sqrt(x1.^2+(x2-4).^2)));
           ufield(ix1,ix2) = (phivec*((hvec.* vecdsdt).')*2*pi/N);
        end
    end
end
num_Sub_matrix = ufield(M/3:M*2/3,M/3:M*2/3);
mean_sub_matrix = mean(num_Sub_matrix, 'all');
real_sub_matrix = exactfield(M/3:M*2/3,M/3:M*2/3);
mean_real_sub = mean(real_sub_matrix,'all');

t=0;
re_ufield = real(ufield*exp(-1i*k*t));
im_ufield = imag(ufield*exp(-1i*k*t));
errorfield = log10(abs(re_ufield-exactfield));

%% 
imagesc(x1field, x2field, exactfield.')
axis xy
colormap turbo
pbaspect([1 1 1])
colorbar
title('Analytical solution of u(x,y)')

%%
imagesc(x1field, x2field, re_ufield.')
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
