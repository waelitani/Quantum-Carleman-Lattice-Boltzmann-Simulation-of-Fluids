function options = quantumOptions(positionmomentum,lastlevelpush,taylorEvolution,taylorOrder)

if isempty(positionmomentum)
    positionmomentum = false;
end

if isempty(lastlevelpush)    
    lastlevelpush = false;
end

if isempty(taylorEvolution)
    taylorEvolution = false;
end

if isempty(taylorOrder)    
    taylorOrder = 0;
end

options = struct();    

options.UsePositionMomentumEmbedding = positionmomentum;
options.IncludeTruncationCorrection = lastlevelpush;
options.TaylorExpandEvolutionOperator = taylorEvolution;
options.TaylorExpansionOrder = round(max(taylorOrder,0));
%add arguments block
%include repeat normalization, return state, return norm
%normalized is not implemented

end
