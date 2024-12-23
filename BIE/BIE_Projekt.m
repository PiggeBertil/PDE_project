N=1000;
tvec = linspace(-pi+2*pi/N, pi, N);
rvec = 3+cos(4*tvec+pi);
rprimvec = -4*sin(4.*tvec+pi);
rbisvec = -16*cos(4.*tvec+pi);
y1 = rvec .*cos(tvec);
y2 = rvec .*sin(tvec);
nu1 = rvec .* cos(tvec) + rprimvec .* sin(tvec);
nu2 = rvec .* sin(tvec) - rprimvec .* cos(tvec);
nu1 = nu1 ./ sqrt( rvec.^2+ rprimvec.^2 );
nu2 = nu2 ./ sqrt( rvec.^2+ rprimvec.^2 );

A = zeros(N,N);
vecdsdt = sqrt(rprimvec.^2+rvec.^2);
for i=1:N
    for j=1:N
        nu = [nu1(i), nu2(i)];
        
        r_i = [y1(i), y2(i)];
        r_j = [y1(j), y2(j)];
        
        param = norm(r_i-r_j)^2;
        namnare = 2*pi .* param;

        A(i,j) = dot(nu, (r_i-r_j)/namnare);
    end
end

for i=1:N
    taljare = rprimvec(i)^2-0.5*rbisvec(i)*rvec(i) + 0.5*rvec(i)^2;
    namnare = 2*pi*(rprimvec(i)^2+rvec(i)^2)^(3/2);
    A(i,i) = taljare/namnare;    
end

imagesc(tvec,tvec,A.')
axis xy
colorbar