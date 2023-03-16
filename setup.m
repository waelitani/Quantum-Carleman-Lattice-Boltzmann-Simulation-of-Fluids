function [evolutionOp,initialState,opa,opad,norm] = setup(funcToSolve,Dt,initialValues,qc,options);

%this setup is limited to real-valued functions with real initial values

[a,ad,q,p,nn] = operators(qc);

[initialState,norm] = initialize(a,ad,q,p,nn,initialValues,qc,options);

[evolutionOp,opa,opad] = linearembedding(a,ad,q,p,nn,funcToSolve,Dt,qc,options,norm);

end