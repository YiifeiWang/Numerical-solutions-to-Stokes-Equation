function out=plot_mat(A)

[m,n] = size(A);

x=1:m;
y=1:n;

figure('Color', 'w')
Handle=pcolor(x, y, full(A));
% For more control options, type: 
% get(Handle); select the option and use set(Handle, 'Option', 'Value')
caxis([min(min(A)) max(max(A))]')
colormap(jet)
CBar_Handle=colorbar('EastOutside');
set(get(CBar_Handle, 'YLabel'), 'String', 'Data Value Scale')
title('Title')
xlabel('x-axis')
ylabel('y-axis')
out=0;

end