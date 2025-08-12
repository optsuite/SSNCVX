%% generate a random positive definte matrix

function A = rand_psd(n)
    A = rand(n);
    A = A'*A + 0.1*eye(n);
end