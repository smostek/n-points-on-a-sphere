function viewManyShapes(s, Opts)
% Visualize any number of shapes in one figure
% Click on the image of any shape to go to 'detail view' and see more info about that
% shape and how it compares to the others
arguments
    s 
    % Shape vector or a cell array with each cell containing shape vectors

    Opts.visType = 'flat'
    % Visualization type to be used for the main display window
    % see shape/see for options

    Opts.detailviews (1,2) cell = {'solid','symmetry'}
    % Visualization types (2) to be used within the detail view
    % see shape/see for options

    Opts.ratings (1,2) double = [1 2];
    % Rating systems to rank the shapes by in the detail view
    % see shape/rate for options
end
if iscell(s)
    isACell=true;
    n=1;
    alls = s;
    s = s{n};
else
    isACell=false;
    alls = {s};
end
if ~isa(s,'shape')
    error('s must be either a shape vector or a cell vector of shape vectors')
end

% Figure
f = gcf;
if strcmp(gcf().UserData,'many')
    clf
elseif ~isempty(gcf().Children)
    f = figure;
end
set(f,'Units', 'normalized', 'Position', [0 0.1 0.8 0.8],'Color','w','UserData','many');
% 'Next' Button:
but(1) = uicontrol('Style','pushbutton','String','>', 'Units','normalized', ...
    'Position', [0.95 0.45 0.03 0.1], 'Callback', @nextcallback,'BackgroundColor','w');
% 'Last' Button:
but(2) = uicontrol('Style','pushbutton','String','<', 'Units','normalized', ...
    'Position', [0.02 0.45 0.03 0.1], 'Callback', @lastcallback,'BackgroundColor','w');
% Main Panel
mp = uipanel('Parent', f, 'Position', [0.06 0 0.84 1],'BackgroundColor','w', ...
    'FontName','Liberation Serif','FontSize',24,'TitlePosition','leftbottom');
% Selector
if isACell
    sel = uicontrol('Parent',f,'Style','popupmenu','Units','normalized', ...
        'Position',[0.01,0.95,0.04,0.03],'String',cellstr(num2str((1:length(alls))')), ...
        'Callback',@selCallback);
end

ratingInstructions = {'Force Score', 'ascend';
                      'Distance Score', 'descend'};

S = length(s);
DetailMode=false;
currentI=1;
TL = viewAll;

function axs = viewAll
    N = length(s);
    set(f,'Tag','on')
    if N <=12
        axs = tiledlayout(mp, ceil(N/4), 4);
    else
        rows = ceil(N/4);
        % sub panel
        axs(3) = uipanel('Parent',mp, 'Position', [0 1-rows/3 1 rows/3],'BackgroundColor','w');
        %scroll bar
        axs(2) = uicontrol('Parent',f, 'Style','slider', 'units', 'normalized', ...
            'Position', [0.9 0.1 0.04 0.8], 'Min', 1-rows/3, 'Max',0, 'Value', 0, ...
            'Callback', {@scbrcallback, axs(3), rows},'BackgroundColor','w');
        axs(1)=tiledlayout(axs(3), rows, 4);
    end
    if ~isACell
        set(but,'Visible','off')
    else 
        set(mp,'Title',['n = ' num2str(n)])
    end
    %input all the axes
    for j=1:N
        ax = nexttile(axs(1));
        set(ax,'UserData',j,'ButtonDownFcn',@clikcallback,'HitTest','on','PickableParts','all')
        s(j).see(Opts.visType, 'axes',ax, 'clear',false)
        set(ax.Children,'Hittest','off')
        title(ax,symmetryName(s(j).symmetry))
    end
    set(f,'Tag','')
end

function viewDetails(I)
    trg = s(I);
    delete(TL)
    set(f,'Tag','on')
    TL = tiledlayout(mp, 2, 2);
    % 'Back' Button
    TL(2)=uicontrol('Style','pushbutton', 'String', 'Return', 'Units','normalized', ...
        'Position', [0.07 0.92 0.08 0.05], 'Callback', @returnCallback);
    set(but,'Visible','on')
    % Label
    cs=CS(trg.cartCoor);
    TL(3)=uicontrol('Style','text','String',{['cs: ' num2str(cs)], ...
        ['com: ' num2str(vecnorm(mean(trg.cartCoor,2)))]}, ...
        'Units', 'normalized', 'Position',[0.45,0.7,0.1,0.05],'BackgroundColor','w');
    ax = nexttile(TL(1));
    trg.see(Opts.detailviews{1}, 'axes',ax)
    ax = nexttile(TL(1));
    trg.see(Opts.detailviews{2},'axes',ax, 'clear',false)

    if any(arrayfun(@(a) isempty(s(a).ratings), 1:S))
        s.rate;
    end
    ax = nexttile(TL(1));
    scorePlot(Opts.ratings(1))

    ax = nexttile(TL(1));
    scorePlot(Opts.ratings(2))
    set(f,'Tag','')
    
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
    elseif isACell && n<length(alls)
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
    elseif isACell && n>1
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

function selCallback(~,~)
    n = sel.Value;
    DetailMode = false;
    delete(TL)
    drawnow
    s = alls{n};
    S = length(s);
    TL = viewAll;
end

end