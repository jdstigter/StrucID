function attributes = StrucIDset(varargin)

if ~isempty(varargin) && isstruct(varargin{1})
    attributes = varargin{1};
    varargin(1)=[];
else
    IntOptions = odeset('AbsTol',1e-12);
    attributes = struct('Obs',true,'Integrator',@ode45,'IntegratorOptions',IntOptions,'Tf',3,...
        'JacobiComplex',true,'TypeOfAnalysis','1','nTraject',1,'kRed',0,'nyMissing',0,...
        'Tsim',9,'IncludeAdjointSens',false);
end

narg=numel(varargin);
for k=1:2:narg
    if isfield(attributes,varargin{k})
        attributes.(varargin{k})=varargin{k+1}; % no check on values is performed here!
    else
        disp('Wrong field name. Please, check!');
    end
end
