function [evolutionOp,opa,opad] = linearembedding(a,ad,q,p,nn,funcToSolve,Dt,tau,qc,options,nbVars,funcNorm)

if nargin < 12
  funcNorm = @(t) 1+0*t;
end

opa = @(x) operatorConstructor({a},[x],nbVars,qc);
opad = @(x) operatorConstructor({ad+(options.IncludeTruncationCorrection)*2^qc*nn},[x],nbVars,qc);
opq = @(x) operatorConstructor({q},[x],nbVars,qc);
% checks only the first differential equation to determine whether we are
% initializing or evolving
%this needs to be generalized
opp = @(x) operatorConstructor({p+(options.IncludeTruncationCorrection)*(isa(funcToSolve{1},'function_handle'))*2^qc*nn},[x],nbVars,qc);

opX = @(x) operatorConstructor({eye(2^qc)-(options.IncludeTruncationCorrection)*2^(-qc)*nn},[x],nbVars,qc);

sts = 2^(qc*nbVars);

zeroOp = zeros(sts,sts);
oneOp = eye(sts);

operator = zeroOp;

if options.UsePositionMomentumEmbedding == true
    pushingOp = @(x) opp(x);
    variableOp = @(x) opq(x);
else
    pushingOp = @(x) opad(x);
    variableOp = @(x) opa(x);
end

for n = 1:1:nbVars
    passVars{n} = feval(variableOp,n);
end

% maxSparsity = 0;

for n = 1:1:nbVars
    if isa(funcToSolve{n},'function_handle')
        toadd = (feval(pushingOp,n)*feval(funcToSolve{n},passVars))*opX(n);
%         maxSparsity = max(max(sum(toadd~=0,2)),maxSparsity);
        operator = operator + toadd;
    else
        toadd =  feval(pushingOp,n)*funcToSolve{n};
%         maxSparsity = max(max(sum(toadd~=0,2)),maxSparsity);
        operator = operator + toadd;
        
    end
end

% maxSparsity

doperator = conj(transpose(operator));

%dynamic time rescaling is currently hard-coded
H = @(t) 0.5*(operator+doperator)%*funcNorm(t);
% H = @(t) operator;%*funcNorm(t);
% H = @(t) (operator)%*funcNorm(t);

% sparsityH = max(sum(H(0)~=0,2))
% maxNormH = norm(H(0))

if options.TaylorExpandEvolutionOperator == false
    if options.UsePositionMomentumEmbedding == true
        evolutionOp = @(t) expm(-1i*Dt*H(t));
    end
else
    if options.UsePositionMomentumEmbedding == true
        evolutionOp = @(t) t*zeroOp;
        for oo = 0:1:options.TaylorExpansionOrder
            evolutionOp = @(t) evolutionOp(t)+(-1i*Dt*H(t))^oo;
        end   
    end
    
end

end
