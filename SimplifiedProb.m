maxN=20;
y=zeros(maxN,1);

figure('Position',[50 250 550 750])
tiledlayout(2,1)
a1=nexttile;
hold on
a2 = nexttile;
hold on
for n = 2:maxN
    x = linearEquilibrium(n);
    scatter(a1,x,n*ones(n,1),50, 'b','filled')
    y=linspace(-1,1,n)';
    plot(a2,y,x-y,'Color',rand(3,1))
    % y(n)=x(2);
end


function x = linearEquilibrium(n)
x = linspace(-1,1, n)';
F = force(x);
cs = sum(abs(F));
while cs > 1e-6
    x = x+F/n^3;
    F = force(x);
    cs = sum(abs(F));
end

    function F = force(x)
        d = x - x';
        F = sum(sign(d)./d.^2,2,'omitmissing');
        F([1 n]) = 0;
    end
end