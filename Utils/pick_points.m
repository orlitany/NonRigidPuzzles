function idx = pick_points(M,hh)


function [txt] = myupdatefcn(~,event_obj)
    pos = get(event_obj,'Position');
    [~,idx] = ismember(pos,M.VERT,'rows');
    txt = {['idx: ',num2str(idx)]};
end

    figure, plot_mesh(M,hh)
    
    h = datacursormode;
    set(h,'SnapToDataVertex','on');
    set(h,'UpdateFcn',@myupdatefcn);
        
    % Use ALT to select multiple points, and the code below to have the
    % points saved into an array.
    
%     f = getCursorInfo(h);
%     
%     a = struct2cell(f);
%     n = size(a,3);
%     
%     idx = zeros(n,3);
%     
%     for i=1:n
%         p = a(2,:,i);
%         idx(i,:) = p{1};
%     end
%     
end
