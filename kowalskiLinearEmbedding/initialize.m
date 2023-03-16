function [state] = initialize(a,ad,q,p,nn,initialValues,tau,qc,nbVars,options)

nbVars = length(initialValues);

% overallNorm = 1;
% 
% for v = 1:1:nbVars
%     overallNorm = overallNorm*sqrt(truncatedNorm(initialValues(v),qc));
% end

for v = 1:1:nbVars
   funcToInitialize{v} = initialValues(v);
end

options.normalizeInitialState = false;
options.DoNotUseDisplacementOperator = false;
options.TaylorExpandEvolutionOperator = false;
options.TaylorExpansionOrder = 0;

initOp = linearembedding(a,ad,q,p,nn,funcToInitialize,1,tau,qc,options,nbVars);

state = zeros(2^(qc*nbVars),1);

% if options(3) == true
%     state(end) = 1;
% else
    state(1) = 1;
% end

state = initOp(0)*state;

end
