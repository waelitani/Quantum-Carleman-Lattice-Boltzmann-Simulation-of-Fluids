function sol = classicalNumericalSolution(funcToSolve,initialValues,nbVars,tend,Dt)

tarray = 0:Dt:tend;

f = sym('f',[nbVars,1]);
for n = 1:1:nbVars
    var{n} = f(n);
end

%initialize solution array
sol = zeros(nbVars,length(tarray));
sol(:,1) = initialValues;
drivingFunction = zeros(nbVars,1);

%evolve the solution
for tt = 2:1:length(tarray)
    for n = 1:1:nbVars
        prev{n} = sol(n,tt-1);
    end
    
    for n = 1:1:nbVars
        drivingFunction(n) = funcToSolve{n}(prev);
    end
    
    sol(:,tt) = sol(:,tt-1)+Dt*drivingFunction;
end

end