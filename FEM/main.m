
%Create A
N = size(p, 2);
A = zeros(N, N);


for element = 1 : size(t, 2)
    nn = t(1:3, element);
    disp("nn:")
    disp(nn)
    A(nn, nn) = A(nn, nn) + 1;
end
for boundryelement = 1 : size(e, 2)
    nn = e(1:2, boundryelement);
    A(nn, nn) = A(nn, nn) + 1;
end