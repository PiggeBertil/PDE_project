function A0 = IntMatrix(nodes)
  % Input: 2 x 3 matrix, node coords as columns
  % Output: 3 x 3 matrix of integrals for stiffness matrix
  e1= nodes(:, 1)- nodes(:, 3);  % choose 3rd node as origin
  e2= nodes(:, 2)- nodes(:, 3);
  basis= [e1, e2];
  dualbasis= inv(basis');       % computes the dual basis
  grads= [dualbasis(:, 1), dualbasis(:, 2), -dualbasis(:, 1)- dualbasis(:, 2)];
  area= det(basis)/2;          % computes the area of the triangle
  A0= grads' * grads * area;    % returns the 9 inner products
end
