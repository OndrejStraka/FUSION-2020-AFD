function innerEdges = generateinneredges(states)

% The dimension of the state
dimStates = length(states);

% Preallocate inner edges
innerEdges = cell(1,dimStates);

% Compute inner edges for each element of the vector state
for n = 1:dimStates
    innerEdges{n} = states{n}(1:end-1) +...
        0.5*(states{n}(2:end) - states{n}(1:end-1));
end