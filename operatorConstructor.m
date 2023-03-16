function operator = operatorConstructor(gates,registers,nbRegisters,qc)

%is built for registers of equal length so far
%each gate acts on a register

nbGates = length(gates);

parfor r = 1:1:nbRegisters
    paddedGates{r} = eye(2^qc);
end

for n = 1:1:nbGates
    targetRegister = registers(n);
    paddedGates{targetRegister} = gates{n};
end

operator = paddedGates{1};

for r = 2:1:nbRegisters
    operator = kron(paddedGates{r},operator);
end

end
% operators need to be tensored for rightmost qubit and onwards