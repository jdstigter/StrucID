function StrucIDmakedirectedgraph(model,options)
% Creates a directed graph of the model in the App
np=numel(model.parNames);
nu=numel(model.inputNames);

Observability = options.Obs;
OutputSelect = model.OutputSelect;
ny=numel(OutputSelect);

varValues=cat(1,model.parValues,model.stateValues);
% parValues=varValues(1:np); % not needed because included in varValues
uInput=model.U;
% uInput=model.U; %inputValues=InputTable.Data.inputValues;
% compute values for algebraic equations
vAlgebra=model.vAlg;
% vAlgebraValues=vAlgebra(stateValues,parValues,inputValues);

% compute adjacency matrix for directed graph construction
fDyn=model.fDyn;
% metafDyn=str2func('xdotTemplate');
hObs=model.Y;
% metahObs=str2func('yTemplate');

% NaN's are substituted with random values
index=find(isnan(varValues));
varValues(index)=0.5+rand(size(index)); stateValues=varValues(np+1:end);
A = admDiffComplex(@(x) xdotTemplate(1,x,uInput,varValues,fDyn,vAlgebra),1,stateValues);
Aadj=double(A~=0);

% Find observable/controllable nodes in the model
if Observability
    C = admDiffComplex(@(x) yTemplate(1,x,uInput,varValues,hObs,vAlgebra,OutputSelect),1,stateValues);
    Cadj=(C~=0); if ny>1, Cadj=any(Cadj); end
    ObsNodes=find(Cadj);
else % controllability
    u0=rand(nu,1);
    C = admDiffComplex(@(u) xdotTemplate(1,stateValues,u,varValues,fDyn,vAlgebra),1,u0)';
    Cadj=(C~=0); if nu>1, Cadj=any(Cadj); end
    ObsNodes=find(Cadj);
end

w=digraph(Aadj'); fig1=figure(1); fig1.Resize = 'off'; H=plot(w);
highlight(H,ObsNodes,'NodeColor','r');