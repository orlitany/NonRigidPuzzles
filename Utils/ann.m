function varargout = ann(method, varargin)

%error(nargchk(3, inf, nargin));

% some predicates
is_normal_matrix = @(x) isnumeric(x) && ndims(x) == 2 && isreal(x) && ~issparse(x);
is_posint_scalar = @(x) isnumeric(x) && isscalar(x) && x == fix(x) && x > 0;
is_switch = @(x) islogical(x) && isscalar(x);
is_float_scalar = @(x) isfloat(x) && isscalar(x); 

% % Xr and Xq
% require_arg(is_normal_matrix(Xr), 'Xr should be a full numeric real matrix');
% require_arg(is_normal_matrix(Xq), 'Xq should be a full numeric real matrix');
%         
% [d, n] = size(Xr);
% require_arg(size(Xq, 1) == d, 'The point dimensions in Xr and Xq are inconsistent.')
% 
% % k
% require_arg(is_posint_scalar(k), 'k should be a positive integer scalar');
% require_arg(k <= n, 'The value k exceeds the number of reference points');         

k = 1;

switch method,
    
    case 'init'
        
        ref_pts = varargin{1};
        varargin = varargin(2:end);
        
    case 'search'
        
        tree = varargin{1};
        query_pts = varargin{2};
        if length(varargin) > 2,
            k = varargin{3};
            varargin = varargin(4:end);
        else
            varargin = varargin(3:end);
        end
       
    case 'deinit'
        
        tree = varargin{1};
        varargin = varargin(2:end);

    case 'close'
        %tree = varargin{1};
        
    otherwise
        
        error('INVALID COMMAND');
        
end


% options
opts = struct( ...
    'use_bdtree', false, ...
    'bucket_size', 1, ...
    'split', 'suggest', ...
    'shrink', 'suggest', ...
    'search_sch', 'std', ...
    'eps', 0, ...
    'radius', 0);

if ~isempty(varargin)
    opts = setopts(opts, varargin{:});
end

%require_opt(is_switch(opts.use_bdtree), 'The option use_bdtree should be a logical scalar.');
require_opt(is_posint_scalar(opts.bucket_size), 'The option bucket_size should be a positive integer.');

split_c = get_name_code('splitting rule', opts.split, ...
                        {'std', 'midpt', 'sl_midpt', 'fair', 'sl_fair', 'suggest'});

if opts.use_bdtree                    
    shrink_c = get_name_code('shrinking rule', opts.shrink, ...
                             {'none', 'simple', 'centroid', 'suggest'});
else
    shrink_c = int32(0);
end

ssch_c = get_name_code('search scheme', opts.search_sch, ...
                        {'std', 'pri', 'fr'});

require_opt(is_float_scalar(opts.eps) && opts.eps >= 0, ...
            'The option eps should be a non-negative float scalar.');

use_fix_rad = strcmp(opts.search_sch, 'fr');
if use_fix_rad
    require_opt(is_float_scalar(opts.radius) && opts.radius > 0, ...
                'The option radius should be a positive float scalar in fixed-radius search');
    rad2 = opts.radius  * opts.radius;
else
    rad2 = 0;
end
        
% main (invoking mexann)

internal_opts = struct( ...
    'use_bdtree', opts.use_bdtree, ...
    'bucket_size', int32(opts.bucket_size), ...
    'split', split_c, ...
    'shrink', shrink_c, ...
    'search_sch', ssch_c, ...
    'knn', int32(k), ...
    'err_bound', opts.eps, ...
    'search_radius', rad2);

%[nnidx, dists] = mexann(Xr, Xq, internal_opts);


switch method,
    
    case 'init'
        
        tree = mexann('createKdTree', ref_pts, internal_opts);
        varargout(1) = { tree };
        
    case 'search'
        
        [nnidx, dists] = mexann('performAnnkSearch', tree, query_pts, internal_opts);
        nnidx = nnidx + 1;        % from zero-based to one-based        
        if nargout >= 2
            dists = sqrt(dists);  % from squared distance to euclidean
            if use_fix_rad
                dists(nnidx == 0) = inf;
            end
        end
        varargout(1) = { nnidx };
        varargout(2) = { dists };
        
    case 'deinit'
        
        mexann('deleteKdTree', tree);
        
    case 'close'        
        mexann('annClose');
        clear mexann;


end








% Auxiliary function

function c = get_name_code(optname, name, names)

require_opt(ischar(name), ['The option ' optname ' should be a string indicating a name.']);

cidx = find(strcmp(name, names));
require_opt(~isempty(cidx), ['The option ' optname ' cannot be assigned to be ' name]);

c = int32(cidx - 1);


function require_arg(cond, msg)

if ~cond
    error('ann_mwrapper:annquery:invalidarg', msg);
end

function require_opt(cond, msg)

if ~cond
    error('ann_mwrapper:annquery:invalidopt', msg);
end




function opts = setopts(opts0, varargin)

if isempty(opts0)
    opts = [];
elseif isstruct(opts0) && isscalar(opts0)
    opts = opts0;        
else
    error('dmtoolbox:setopts:invalidarg', ...
        'opts0 should be either a struct scalar or empty.');
end

if nargin > 1
    fparam = varargin{1};
    if isstruct(fparam)
        if nargin > 2
            error('dmtoolbox:setopts:invalidarg', ...
                'No input arguments are allowed to follow the struct parameter');
        end
        params = fparam;
    elseif iscell(fparam)
        if nargin > 2
            error('dmtoolbox:setopts:invalidarg', ...
                'No input arguments are allowed to follow the cell parameter');
        end
        params = fparam;
    elseif ischar(fparam)
        params = varargin;
    else
        error('dmtoolbox:setopts:invalidarg', 'The input argument list is illegal.');
    end    
else
    return;
end

%% main delegate

if iscell(params)
    opts = setopts_with_cell(opts, params);
else
    opts = setopts_with_struct(opts, params);
end


%% core functions

function opts = setopts_with_cell(opts, params)

names = params(1:2:end);
values = params(2:2:end);
n = length(names);

if length(values) ~= n
    error('dmtoolbox:setopts:invalidarg', 'The names and values should form pairs');
end

for i = 1 : n
    opts.(names{i}) = values{i};
end

function opts = setopts_with_struct(opts, params)

fns = fieldnames(params);
n = length(fns);

for i = 1 : n
    fn = fns{i};
    opts.(fn) = params.(fn);
end



