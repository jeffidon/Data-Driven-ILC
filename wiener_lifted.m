function lifted_Y = wiener_lifted(G, Y, ny, f)

% %% error checking
% if G(1) ~= 1
%     error('the first element of G must be 1!!!');
% end

%% process by G
lifted_Y = [];
for ii = 1:f
    y = Y( (ii-1)*ny+1 : ii*ny );
    z = [];
    for jj = 1:length(G)
        z = [ z; y.^G(jj) ];
    end
    lifted_Y = [lifted_Y; z];
end