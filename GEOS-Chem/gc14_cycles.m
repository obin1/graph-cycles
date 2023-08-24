Badj = readmatrix("GC14ReciprocityMatrixSparse.csv")
S = sparse(Badj(:,1),Badj(:,2),Badj(:,3))

cycleMax = 12
elapsedTime = zeros(cycleMax,1);
for i = 1:6
    [primes, elapsedTime(i)] = CycleCount(S,i)
end