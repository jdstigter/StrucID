function [model,options] = StrucIDreadmodel(filename)
% filename refers to a text file that contains the model definition in a
% specific format (see example file)
% The result Model is a structure with fields

% run ADiMAT startup script for automatic differentiation. This is to make
% sure that the ADiMAT routines are available when running the sensitivity
% computations
if (~exist('admDiffComplex.m','file') && ~isdeployed)
    run('../adimat/ADiMat_startup');
end

%% 1. ENTER MODEL _ GO AND FETCH FILE TO ATTACH

if ~ischar(filename)
    error('Error. Argument to this function must be a string of characters.');
else
    delimiterIn = '!'; %delimiter used in code to detect different sections
    tree = importdata(filename,delimiterIn);
end

close all; % close all figures if these are still open
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[algebraicEqns,algebraicEqnsTxt,ODEs,ODEsTxt,ICStateEqns,ICStateEqnsTxt,icNames,stateNames,stateValues,...
    parNames,parValues,UnknownVarIndex,inputNames,inputEqns,inputEqnsTxt,outputNames,...
    outputEqns,outputEqnsTxt] = process_model_file(tree);

ny = numel(outputNames);

% Build functions of the algebraic equations, ODEs, Yout, Inputs, and x0.
% The function handles + definitions will be passed on as outputs/fields to
% the model structure variable with name 'model'

% algebraic relations
vAlgebra = str2func(['@(x,p,u)[',strjoin(algebraicEqns),']']);

% ODEs or dynamics f(t,x,u,p,v)
Xdot = str2func(['@(t,x,u,p,v)[',strjoin(ODEs),']']);

% Initial conditions
x0State = str2func(['@(p)[',strjoin(ICStateEqns),']']);

% Input signal
U = str2func(['@(t,p)[',strjoin(inputEqns),']']);

% Output signal
Y = str2func(['@(t,x,u,p,v)[',strjoin(outputEqns),']']);

np = numel(parNames);
nx = numel(stateNames);
FixedVarIndex = true(1,np+nx); % We fix all variables to their values

ModelEqnsTxt = struct('AlgebraicEqnsTxt',{algebraicEqnsTxt},'ODEsTxt',{ODEsTxt},...
    'ICStateEqnsTxt',{ICStateEqnsTxt},'YTxt',{outputEqnsTxt},'UTxt',{inputEqnsTxt});

model = struct('fDyn',Xdot,'vAlg',vAlgebra,'x0',x0State,'U',U,'Y',Y,...
    'inputNames',{inputNames},'outputNames',{outputNames},...
    'icNames',{icNames},'stateNames',{stateNames},'stateValues',stateValues,...
    'parNames',{parNames},'parValues',parValues,'EqnsTxt',ModelEqnsTxt,...
    'UnknownVarIndex',UnknownVarIndex,'FixedVarIndex',FixedVarIndex,...
    'OutputSelect',(1:ny)');

options = StrucIDset;





