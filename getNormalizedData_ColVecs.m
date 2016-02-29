function normTraces = getNormalizedData_ColVecs(traces)
% return normalized column-vectors

len = size(traces,1);
normVals = (sum(traces.*traces)).^0.5;

normTraces = traces ./ repmat(normVals,len,1);

end