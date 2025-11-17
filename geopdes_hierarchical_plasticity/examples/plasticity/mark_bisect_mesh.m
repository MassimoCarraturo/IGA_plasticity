function marked = mark_bisect_mesh(hmsh_scalar)

marked = cell(hmsh_scalar.nlevels,1); % bisect every element
for ilev = 1:hmsh_scalar.nlevels
    marked{ilev} = (hmsh_scalar.active{ilev})';
end

end