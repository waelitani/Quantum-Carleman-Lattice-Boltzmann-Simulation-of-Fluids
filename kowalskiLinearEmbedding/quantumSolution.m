function [sol] = quantumSolution(op,opa,opad,initialState,nbVars,tend,Dt,tau,funcNorm,options)

state = initialState;

tarray = 0:Dt:tend;

sol = zeros(nbVars,length(tarray));

sol(:,1) = expectation(opa,opad,state,nbVars,options);

for tind = 2:1:length(tarray)
    state = op((tind-1)*Dt)*state;
    sol(:,tind) = expectation(opa,opad,state,nbVars,options)/(funcNorm((tind-1)*Dt))^2;
end

end
