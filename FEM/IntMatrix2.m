function B0 = IntMatrix2(nodes)
  % Input: 2 x 3 matrix, node coords as columns
  % Output: 3 x 3 matrix of integrals for stiffness matrix
  e1= nodes(:, 1)- nodes(:, 3);  % choose 3rd node as origin
  e2= nodes(:, 2)- nodes(:, 3);
  basis= [e1, e2];
  multiple = det(basis)/24;          % computes the area of the triangle
  B0 = (ones(3,3) + eye(3, 3))*multiple;  % returns the 9 inner products
end
