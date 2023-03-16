function value = expectation(opa,opad,state,nbVars,options)

value = zeros(nbVars,1);

if options.UsePositionMomentumEmbedding == true

    q = @(x)(1/sqrt(2))*(opa(x)+opad(x));

    for v = 1:1:nbVars
            value(v) = conj(state')*q(v)*state;
            value(v) = norm(value(v));
    end
    
else
    
    for v = 1:1:nbVars
%             value(v) = conj(state')*opad(v)*opa(v)*state;
            value(v) = conj(state')*opa(v)*state;
%             value(v) = state(2);
            value(v) = norm(value(v));

    end
    
end

end
