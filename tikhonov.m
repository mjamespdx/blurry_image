function Xapp = tikhonov(lambda,Dapp)

% This function uses Tikhonov regularization to approximate an original
% matrix from a given noisy matrix.

% lambda - regularization parameter
% Dapp - noisy data matrix

[Dm,Dn] = size(Dapp); % size constraints for function
L = 0.45; % given parameter

% create B, A, I matrices
B = zeros(Dm,Dm);
for i = 1:Dm
    B(i,i) = 1 - 2*L;
    if (i <= Dm-1)
	B(i,i+1) = L;
	B(i+1,i) = L;
    end
end

A = B^25; % non-singular transformation matrix
I = eye(size(A)); % Dm x Dm identity matrix

% solve approximate X original matrix
for j = 1:Dn
    d = Dapp(:,j);
    x = (A’*A + lambda^2*I)\A’*d;
    Xapp(:,j) = x;
end

end