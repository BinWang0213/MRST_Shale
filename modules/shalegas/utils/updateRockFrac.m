function G = updateRockFrac(G, K_frac, frac_type,varargin)
%{
    Set up fracture permeability
    All of input parameters are distributed into each sub-function

    Author: Bin Wang(binwang.0213@gmail.com)
    Date: Dec.2018
%}


opt = struct('permtype', 'homogeneous');
opt = merge_options(opt, varargin{:});

if ~isempty(varargin) && strcmp(opt.permtype,'homogeneous')==0 && strcmp(opt.permtype,'heterogeneous')==0
    fprintf(['Wrong argument string. Valid Arguments: ''Homogeneous'' or ''heterogeneous''.\n',...
        'Setting homogeneous permeability in fractures...\n']);
end

for i = 1:numel(fieldnames(G.FracGrid))
    Gf = G.FracGrid.(['Frac',num2str(i)]);
    if strcmp(frac_type,'HydraulicFracs')
        if(i>G.NumHFs),continue,end
    end
    if strcmp(frac_type,'NaturalFracs')
        if(i<=G.NumHFs),continue,end
    end
    G.FracGrid.(['Frac',num2str(i)]).rock.perm = ones(Gf.cells.num, 1)*darcy*K_frac;
end

end