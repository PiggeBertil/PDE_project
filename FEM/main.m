
%Create A and B
N = size(p, 2);
A = zeros(N, N);
B = zeros(N, N);

%iterate over all nodes
for element = 1 : size(t, 2)
    %Calculate A_nn and B_nn for adjecent node.
    nn = t(1:3, element);
    A(nn, nn) = A(nn,nn) + IntMatrix( p(:, nn) );
    B(nn, nn) = B(nn, nn) + IntMatrix2( p(:, nn) );
end


nn= size(p, 2);
intnodes= setdiff(1:nn, e(1, :)); % list of non-boundary nodes 
ADir= A(intnodes, intnodes); % removes boundary rows/columns
A=ADir;
BDir= B(intnodes , intnodes);
B=BDir;


%Minimize det(A-lamda*B) to recieve our eigenvalues approximations
% Compute the generalized eigenvalues
[eigenvectors, eigenvalues] = eig(A, B);

% Extract eigenvalues from the diagonal matrix
eigenvalues = diag(eigenvalues);

% Sort the eigenvalues by their absolute value
[sorted_eigenvalues, indices] = sort(abs(eigenvalues));

% Select the 50 smallest eigenvalues
smallest_eigenvalues = sorted_eigenvalues(1:50);

% Display results
disp('The 50 smallest eigenvalues are:');
%plot(smallest_eigenvalues)

%plot u for all f:s
%pdeplot(p,e,t,XYData=eigenvectors(:,10))

%% b) calculate errrors

analytic_eigenvalues = [];

for n = 1:50
    for m = 1:50
        analytic_eigenvalues = [analytic_eigenvalues; pi^2*(n^2/4 + m^2)];
    end
end

analytic_eigenvalues = sort(analytic_eigenvalues);

errors = abs(analytic_eigenvalues(1:50)-smallest_eigenvalues);

max_error = max(errors);

disp(max_error(1));


%% Plot eigenvalues and weyl's approximation

AREA = 2*pi;
CIRCOMFRENCE = 2^(3/2)*pi;

js = 1:50;

weyls_eigenvalues_linear = js*4*pi/AREA;

weyls_eigenvalues_curve = weyls_eigenvalues_linear + (CIRCOMFRENCE/AREA^(3/2))*(4*pi*js).^0.5;

hold on
plot(smallest_eigenvalues);
plot(weyls_eigenvalues_curve);
plot(weyls_eigenvalues_linear);

legend('smallest eigenvalues','weyls curvelinear aproximation', 'weyls linear aproximation');

title('Eigenvalues for D_4');