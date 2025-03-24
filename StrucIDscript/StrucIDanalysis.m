function [Iden,NumberZeroSingularValues,GapSize,CorParIndex,thetaReal,Message] = StrucIDanalysis(model,options)
% main analysis for all cases (Identifiability, Observability,
% Controllability, Model Reduction, Minimial Output Set)

% Set defaults
ProbeNetwork=false; % default is no probing for a MOS
nck=1; % 1 is the default value for number of possible sensor combinations.
StateReduction=false; % default is no state reduction

nTraject=options.nTraject; % number of trajectories to be computed
Concatenate = (nTraject>1); % concatenation of trajectories yes/no

ChoiceAnalysis=options.TypeOfAnalysis;
switch ChoiceAnalysis
    case '2' % state reduction
        StateReduction=true;
        kRed=options.kRed; % order reduced model kRed<=nx
    case '3' % MOS
        ProbeNetwork=true;
        nyMissing=options.nyMissing; % number of missing sensors
end

% Pick up variables from the model and option structures
OutputIndex = model.OutputSelect;
% return if no sensors have been selected
if (isempty(OutputIndex) && not(ProbeNetwork))
    Iden=4; NumberZeroSingularValues=0; GapSize=0; CorParIndex=[];
    Message = [];
    return;
end

Observability = options.Obs;
Controllability = ~Observability;

nx=numel(model.stateNames);
np=numel(model.parNames);
nu=numel(model.inputNames);
nth = np + nx;

if ProbeNetwork
    ny=nx; % only directly measured states are considered when probing a network for a MOS
else
    ny=numel(OutputIndex);
end

% check if computation Adjoint Sensitivities is actually needed
ComputeForwardSens = (ny>0 && Observability);
ComputeAdjointSens = (nu>0 && (Controllability || options.IncludeAdjointSens));

% return if no inputs are present
if (nu==0 && Controllability)
    Iden=0; NumberZeroSingularValues=0; GapSize=0; CorParIndex=[];
    Message = sprintf('%s\n\n','Cannot analyse controllability if no inputs are present!');
    return;
end

% ParIndex is a ROW vector of indices that indicate the unknown parameters
% in the total set of parameters
ParIndex=model.UnknownVarIndex; nthA=numel(ParIndex);
varNames=cat(1,model.parNames,model.icNames);
varValues=cat(1,model.parValues,model.stateValues);

JacobiComplex=options.JacobiComplex;
if JacobiComplex
    odefnc=@meta;
    odefncAdj=@metaAdjoint;
else
    odefnc=@metaFD;
    odefncAdj=@metaAdjointFD;
end

% 6.2 INTEGRATION OF ODEs + SENSITIVITIES
% pickup tolerance integration
opts = options.IntegratorOptions;
% pickup integrator choice
integrator = options.Integrator;

% function handles for dynamics, IC, input, and algebraic relations
fDyn=model.fDyn; % function handle to xdot=f(x,u,th,v)
modelx0=model.x0; % function handle to initial conditions
vAlgebra=model.vAlg; % function handle to algebraic eqns
uInput=model.U; % function handle to input signal/function

% storage variables for sensitivities in forward and adjoint eqns
dxdthRows=cell(nTraject,1);
dxdthAdjointRows=cell(nTraject,1);
StateTraj=cell(nTraject,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.1 BEGIN ANALYSIS
thetaNom=varValues;
% Randomly draw theta from a uniform distribution
% FIX the parameters that have been specified numerically in model file
ind=find(isnan(thetaNom)); nNaN=numel(ind);
if (nNaN>0)
    thetaLow=0.5*ones(nth,1); thetaHigh=1.5*ones(nth,1);
    thetaNom(ind)=thetaLow(ind)+rand(nNaN,1).*(thetaHigh(ind)-thetaLow(ind));
end

if Concatenate
    thetaReal=zeros(nth,nTraject);
    FixedVarIndex=model.FixedVarIndex;
    for j=1:nTraject
        randomMultipliers=0.8+0.4*rand(nth,1);
        randomMultipliers(FixedVarIndex)=1;
        thetaReal(:,j)=randomMultipliers.*thetaNom;
    end
else
    thetaReal=thetaNom;
end

disp('Clock started');
tic;

Tf=options.Tf;
flagTf=true;
while flagTf
    lastwarn('',''); % reset and clean-up last warning message

    try
        % integration times
        T0=0;
        Nt=max(30,nth+1);
        timeForward=linspace(T0,Tf,Nt);
        timeBackward=linspace(Tf,T0,Nt);

        for j=1:nTraject
                theta=thetaReal(:,j);
                x0=modelx0(theta);
                x0Sens=admDiffComplex(@(th) modelx0(th),1,theta);
                x0Aug=cat(1,x0,x0Sens(:));
                if ComputeForwardSens
                    [~,X] = integrator(@(t,x) odefnc(t,x,uInput,theta,fDyn,...
                        vAlgebra,[nx nth]),timeForward,x0Aug,opts);
                else % trajectory x(t) is needed for Adjoint sensitivities computation
                    [~,X] = integrator(@(t,x) xdotTemplate(t,x,uInput,theta,fDyn,...
                        vAlgebra),timeForward,x0,opts);
                end

                dxdthRows{j}=X(:,(nx+1):end);
                StateTraj{j}=X(:,1:nx);

            if ComputeAdjointSens
                xf=X(end,1:nx)'; % start backwards integration with final state xf
                theta(np+1:end)=xf;
                xfSens=eye(nx);
                xfAug=cat(1,xf,xfSens(:));

                [~,X] = integrator(@(t,x) odefncAdj(t,x,uInput,theta,fDyn,...
                    vAlgebra,[nx nx]),timeBackward,xfAug,opts);

                dxdthAdjointRows{j}=X(:,(nx+1):end);
            end
        end

        [~,warnID]=lastwarn();
        if ~isempty(warnID)
            error('Integration Tolerance not met with reasonable step-size');
        end

        flagTf = false; % successfull integration performed!
    catch
        Tf = Tf/2;
        disp(['Final integration time was halved (because of accuracy) to a value Tf = ',num2str(Tf)]);
    end
end

% Analysis of sensitivities (complete the mapping from dxdth to dydth)
switch ChoiceAnalysis
    % make observation function and (if not probing)
    % draw directed graph for sensor combination chosen
    case {'1','2'}
        h=model.Y;
    case '3' % case network probing
        allSensors=1:nx;
        SelectedSensors=OutputIndex';
        availableSensors=setdiff(allSensors,SelectedSensors);
        nSensors=numel(availableSensors);
        if (nSensors>nyMissing)
            Trials=nchoosek(availableSensors,nSensors-nyMissing);
        else
            Iden=3; NumberZeroSingularValues=0; GapSize=0; CorParIndex=[];
            Message =[];
            return;
        end
        [nck,~]=size(Trials);
end

Symmetries=[];
% Build the output sensitivity matrix
for comb=1:nck
    dydthAnalysisAllTraj=[];
    dudx0AnalysisAllTraj=[];
    for k=1:nTraject
        dydthAnalysis=[];
        dudx0Analysis=[];
        if ProbeNetwork
            sensorIndex=sort([Trials(comb,:),SelectedSensors]);
            ny=numel(sensorIndex);
        else
            sensorIndex=model.OutputSelect;
        end
        theta=thetaReal(:,k);

        for j=1:Nt
            if ComputeForwardSens
                if ProbeNetwork
                    dxdth=reshape(dxdthRows{k}(j,:),nx,nth);
                    dydth=dxdth(sensorIndex,:); % since we measure states directly!
                else
                    X=StateTraj{k}(j,:);
                    dhdx=admDiffComplex(@(x) yTemplate(timeForward(j),x,uInput,theta,h,vAlgebra,sensorIndex),1,X);
                    dxdth=reshape(dxdthRows{k}(j,:),nx,nth);
                    dhdth = admDiffComplex(@(th) yTemplate(timeForward(j),X,uInput,th,h,vAlgebra,sensorIndex),1,theta);
                    dydth = dhdx*dxdth + dhdth;
                end
                dydthAnalysis=cat(1,dydthAnalysis,dydth(:,ParIndex));
            end
            % compute Adjoint sensitivities if input signals are available and if they are needed
            if ComputeAdjointSens
                X=StateTraj{k}(j,:); u0=uInput(timeForward(j),theta);
                dfdu=admDiffComplex(@(u) xdotTemplate(timeForward(j),X,u,theta,fDyn,vAlgebra),1,u0);
                dxdx0=reshape(dxdthAdjointRows{k}(j,:),nx,nx);
                dudx0=dfdu'*dxdx0';
                dudx0Analysis=cat(1,dudx0Analysis,dudx0);
            end
        end
        % append sensitivities for each trajectory
        if ComputeForwardSens
            dydthAnalysisAllTraj=cat(1,dydthAnalysisAllTraj,dydthAnalysis);
        end
        if ComputeAdjointSens
            dudx0AnalysisAllTraj=cat(1,dudx0AnalysisAllTraj,dudx0Analysis);
        end
    end

    switch ChoiceAnalysis
        case {'1','3'}
            if Observability
                [~,S,V]=svd(dydthAnalysisAllTraj,'econ');
            elseif Controllability
                [~,S,V]=svd(dudx0AnalysisAllTraj,'econ');
            end
        case '2' % 'balanced' reduction
            if options.IncludeAdjointSens
                [~,S,V]=svd(cat(1,dydthAnalysisAllTraj,dudx0AnalysisAllTraj),'econ');
            else
                [~,S,V]=svd(cat(1,dydthAnalysisAllTraj),'econ');
            end
    end

    % Process SVD results
    SingularValues=diag(S)';
    if Observability
        rankSens=sum(SingularValues>(max(size(dydthAnalysis))*eps(S(1)))); % from Matlab rank documentation
    else
        rankSens=sum(SingularValues>(max(size(dudx0Analysis))*eps(S(1)))); % from Matlab rank documentation
    end

    NumberZeroSingularValues=nthA-rankSens;

    logSingularValues=log10(SingularValues);
    logSingularValues(isinf(logSingularValues))=log10(eps); % for hard zeros
    diffLogSingVal=abs(diff(logSingularValues));

    ThresholdGapSize=0; % adjust this criterion if needed
    IndexSingularValuesAfterGap=nthA-NumberZeroSingularValues+1;

    GapSize=max(diffLogSingVal);

    if (NumberZeroSingularValues>0) && (GapSize>ThresholdGapSize)
        Iden=2;
    else
        Iden=1;
    end

    nullspace=V(:,IndexSingularValuesAfterGap:end);

    if (NumberZeroSingularValues==1)
        CorParIndex=find(abs(nullspace')>2e-2);
    elseif (NumberZeroSingularValues>1)
        CorParIndex=find(any(abs(nullspace')>2e-2));
    else
        CorParIndex=[];
    end

    switch ChoiceAnalysis
        case '1'
            if NumberZeroSingularValues==0
                Message = sprintf('%s\n','No Non-Identifiable Parameters Found');
            else
                d1 = sprintf('%s\n\n','Non-Identifiable Parameters: '); d2=[];
                for k=1:numel(CorParIndex)
                    d2=cat(2,d2,sprintf('%s ',varNames{ParIndex(CorParIndex(k))}));
                end
                Message=cat(2,d1,d2,sprintf('\n\n'));
            end
        case '2'
            Message='';
        case '3'
            if (rankSens<nthA)
                RowTable=table;
                allSensors=1:nx;
                RowTable.MissingSensors={setdiff(allSensors,sensorIndex)}; % missing sensors
                RowTable.RankDef=nthA-rankSens; % computed rank deficiency
                RowTable.CorPar={CorParIndex}; % correlated parameters {strjoin(varNames(ParIndex(CorParTrue)))};
                RowTable.nPar=numel(RowTable.CorPar{1});
                Symmetries=cat(1,RowTable,Symmetries);
            end
    end
end

toc

% push results up into the Workspace
assignin('base','dydth',dydthAnalysisAllTraj);
assignin('base','V',V);
assignin('base','Sdiag',diag(S));
assignin('base','theta',theta);
assignin('base','T',timeForward);
    
switch true
    case StateReduction
        assignin('base','kRed',kRed);
    % Simulate the REDUCED model
    theta=thetaNom;

    yRed=[]; yFull=[];

    x0=modelx0(theta); % in case IC depends on theta
    z0=V(:,1:kRed)'*x0;

    NtSim=50; % number of points on which solution is generated
    timeForward=linspace(0,options.Tsim,NtSim);
    [~,X] = integrator(@(t,x) xdotTemplate(t,x,uInput,theta,fDyn,vAlgebra),timeForward,x0,opts);
    [~,Z] = integrator(@(t,z) zdotReduceTemplate(t,z,uInput,theta,fDyn,vAlgebra,V,kRed),timeForward,z0,opts);

    yRedStep=zeros(NtSim,ny); yFullStep=zeros(NtSim,ny);
    for j=1:NtSim
        yRedStep(j,:)=yReduceTemplate(timeForward(j),Z(j,:)',uInput,theta,h,vAlgebra,V,kRed,sensorIndex)';
        yFullStep(j,:)=yTemplate(timeForward,X(j,:)',uInput,theta,h,vAlgebra,sensorIndex);
    end
    yRed=cat(1,yRed,yRedStep);
    yFull=cat(1,yFull,yFullStep);

% cast symmetries variable, that contains all detected correlations for each sensor combination, into the workspace
    case ProbeNetwork
        if not(isempty(Symmetries))
            assignin('base','sym',sortrows(Symmetries,{'RankDef','nPar'}));
            d1=sprintf('Missing Sensor(s) \t Non-Identifiable Parameters \n\n');
            CorParNames=cellfun(@(x) sprintf('%s ',varNames{x}),Symmetries.CorPar,'UniformOutput',false);
            MissingSensorNames=cellfun(@(ind) sprintf('%s, ',varNames{np+ind}),Symmetries.MissingSensors,'UniformOutput',false);
            % add a colon and parentheses to MissingSensorNames
            MissingSensorNames=cellfun(@(S) S(1:end-2), MissingSensorNames, 'Uniform', 0);
            MissingSensorNames=strcat('(',MissingSensorNames,'):');
            d=transpose(horzcat(MissingSensorNames,CorParNames));
            d2=sprintf('%s \t %s \n',d{:});
        else
            d1=sprintf('No Symmetries found with %2d missing sensor(s).\n ',nyMissing);
            d2=[];
        end
        Message = cat(2,d1,d2,sprintf('%s\n',''));
end

% 8.5 identifiability graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
logSingularValues(isinf(logSingularValues))=log10(eps);

switch ChoiceAnalysis
    case '1'
        fig2=figure(2); fig2.Resize='off';
        subplot(2,1,1);
        plot(logSingularValues,'.','MarkerSize',28);
        xlabel('Index singular values (\sigma_i) in descending order');
        ylabel('Log_{10}(\sigma_i)');
        title('Singular Values on a Logarithmic Scale','FontWeight','bold');
        axis([0 nthA+1 -Inf Inf]);
        xticks(0:nthA+1);
        minS=min(SingularValues(:,end));
        textbp(['\sigma_{',num2str(nthA),'}=',num2str(minS,'%10.5e\n')],'Fontname','FixedWidth');

        subplot(2,1,2);
        if NumberZeroSingularValues>=1
            stem(V(:,IndexSingularValuesAfterGap:end));
        else
            stem(V(:,end));
        end
        axis([0 nthA+1 -1 1])
        xticks(0:nthA+1);
        labels=cell(nthA+2,1);

        if nthA<15
            for i=1:nthA
                labels{i+1}=varNames{ParIndex(i)};
            end
            xticklabels(labels);
        end

        xlabel('Parameters');
        ylabel('Components in Singular Vector(s)');
        if (NumberZeroSingularValues<2)
            title('Last Column of V','FontWeight','bold');
        else
            title(['Last ', num2str(NumberZeroSingularValues) ,' Columns of V'],'FontWeight','bold');
        end
    case '2'
        fig2=figure(2); fig2.Resize='off';
        plot(timeForward,yFull,timeForward,yRed,'ro')%,'MarkerSize',3);
        title(['Output Reduced Model with k_{Red} = ',num2str(kRed)]);
        assignin('base','time',timeForward);
        assignin('base','yRed',yRed);
        assignin('base','yFull',yFull);

%         dfdx = admDiffComplex(@(xVar) xdotTemplate(0,xVar,uInput,thetaNom,fDyn,vAlgebra),1,x0);
%         fig3=figure(3); fig3.Resize='off';
%         plot(eig(dfdx),'bx');
%         title('Eigenvalues Jacobi matrix at t_0')
%         xlabel('Real Part');
%         ylabel('Imaginary Part');
end