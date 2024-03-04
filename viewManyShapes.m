function viewManyShapes(s, Opts)
% Options: vistype, detailviews, ratings, title
arguments
    s
    Opts.visType = 'normal'
    Opts.detailviews (1,2) cell = {'solid','area'}
    Opts.ratings (1,2) double = [1 2];
    Opts.title = 'name'
end
Opts.title = validatestring(Opts.title, {'name','fingerprint','none'});

% Figure
f = gcf;
clf
set(f,'Units', 'normalized', 'Position', [0 0.1 0.8 0.8]);
% 'Next' Button:
uicontrol('Style','pushbutton','String','>', 'Units','normalized', ...
    'Position', [0.95 0.45 0.03 0.1], 'Callback', @nextcallback);
% 'Last' Button:
uicontrol('Style','pushbutton','String','<', 'Units','normalized', ...
    'Position', [0.02 0.45 0.03 0.1], 'Callback', @lastcallback);
% Main Panel
mp = uipanel('Parent', f, 'Position', [0.06 0 0.84 1]);

ratingInstructions = {'Force Score', 'ascend';
                      'Area Score',  'descend';
                      'Uniqueness Score','ascend';
                      'Distance Score', 'descend'};
n=1;
if isa(s, 'cell')
    alls = s;
    s = s{n};
else
    alls = {s};
end
S = length(s);
DetailMode=false;
currentI=1;
TL = viewAll;

function axs = viewAll
    N = length(s);
    if N <=12
        axs = tiledlayout(mp, ceil(N/4), 4);
    else
        rows = ceil(N/4);
        % sub panel
        axs(3) = uipanel('Parent',mp, 'Position', [0 1-rows/3 1 rows/3]);
        %scroll bar
        axs(2) = uicontrol('Parent',f, 'Style','slider', 'units', 'normalized', ...
            'Position', [0.9 0.1 0.04 0.8], 'Min', 1-rows/3, 'Max',0, 'Value', 0, ...
            'Callback', {@scbrcallback, axs(3), rows});
        axs(1)=tiledlayout(axs(3), rows, 4);
    end
    %input all the axes
    for j=1:N
        ax = nexttile(axs(1));
        ax.UserData = j;
        ax.ButtonDownFcn = @clikcallback;
        ax.HitTest = 'on'; ax.PickableParts = 'all';
        s(j).see(Opts.visType, 'axes',ax)
        if strcmp(Opts.title, 'name')
            title(ax, s(j).name(:,1))
        elseif strcmp(Opts.title, 'fingerprint')
            title(ax, cellfun(@num2str, s(j).fingerprint,'un',0))
        end
    end
end

function viewDetails(I)
    trg = s(I);
    delete(TL)
    TL = tiledlayout(mp, 2, 2);
    % 'Back' Button
    TL(2)=uicontrol('Style','pushbutton', 'String', 'Return', 'Units','normalized', ...
        'Position', [0.07 0.92 0.08 0.05], 'Callback', @returnCallback);
    % Label
    [cs,T]=CS(trg.cartCoor);
    TL(3)=uicontrol('Style','text','String',{['cs: ' num2str(cs)], ...
        ['ts: ' num2str(mag(sum(T,1)))], ['com: ' num2str(mag(mean(trg.cartCoor,1)))]}, ...
        'Units', 'normalized', 'Position',[0.45,0.7,0.1,0.05]);
    ax = nexttile(TL(1));
    trg.see(Opts.detailviews{1}, 'axes',ax)
    ax = nexttile(TL(1));
    trg.see(Opts.detailviews{2},'axes',ax)

    if any(arrayfun(@(a) isempty(s(a).ratings), 1:S))
        s.rate;
    end
    ax = nexttile(TL(1));
    scorePlot(Opts.ratings(1))

    ax = nexttile(TL(1));
    scorePlot(Opts.ratings(2))
    
        function scorePlot(scoreI)
            scores = sort(arrayfun(@(a) s(a).ratings(scoreI), 1:S), ratingInstructions{scoreI,2});
            colr = zeros(S,3);
            colr(find(scores==trg.ratings(scoreI),1),:)= [1 0 1];
            scatter(ax, 1:S, scores, 50, colr, 'filled');
            title(ax, [ratingInstructions{scoreI,1} ' = ' num2str(trg.ratings(scoreI))])
        end

    end

function scbrcallback(src, ~, arg1, arg2)
    %Scroll bar
    val = src.Value;
    set(arg1, 'Position', [0 1-arg2/3-val 1 arg2/3])
end

function clikcallback(src,~)
    DetailMode = true;
    currentI = src.UserData;
    viewDetails(currentI)
end

function nextcallback(~,~)
    if DetailMode
    currentI = rem(currentI,length(s))+1;
    viewDetails(currentI);
    elseif n<length(alls)
        n = n+1;
        delete(TL)
        drawnow
        s = alls{n};
        S = length(s);
        TL = viewAll;
    end
end

function lastcallback(~,~)
    if DetailMode
        currentI = mod(currentI-2,length(s))+1;
        viewDetails(currentI);
    elseif n>1
        n = n-1;
        delete(TL)
        drawnow
        s = alls{n};
        S = length(s);
        TL = viewAll;
    end
end

function returnCallback(~,~)
    delete(TL)
    DetailMode=false;
    TL = viewAll;
end

end