% PottersWheel model definition file

function m = cascade2level()

m = pwGetEmptyModel();

%% Meta information

m.name        = 'Two-level cascade';
m.description = '';
m.authors     = {'MD'};
m.dates       = {'2015-03-11'};
m.modelFormat = 3;

%% Default sampling time points
m.t = 0:0.1:100;
m.tStart = 0;
m.tPlotStart = 0;


%% X - Dynamic variables
% m = pwAddX(m, *ID, *startValue, fitSetting, minValue, maxValue, unit, compartment, name, description, typeOfStartValue, designerProps, classname)
m = pwAddX(m, 'g1',  1, 'global', .1, 10, [], 'cytosol', []  , []  , [] , [] );
m = pwAddX(m, 'g1p', 0, 'fix'   , 0 ,  [], [], 'cytosol', []  , []  , [] , [] );

m = pwAddX(m, 'g2',  1, 'global', .1, 10, [], 'cytosol', []  , []  , [] , [] );
m = pwAddX(m, 'g2p', 0, 'fix'   , 0 ,  [], [], 'cytosol', []  , []  , [] , [] );

% negative feedback component: the phosphatase
m = pwAddX(m, 'X',    1, 'global', .1, 10, [], 'cytosol', []  , []  , [] , [] );
m = pwAddX(m, 'Xact', 0, 'fix'   , 0 ,  [], [], 'cytosol', []  , []  , [] , [] );


%% K - Dynamic parameters
% m = pwAddK(m, *ID, *value, fitSetting, minValue, maxValue, unit, name, description)
m = pwAddK(m, 'v1',  0.5, 'global', .1, 10);
m = pwAddK(m, 'v2',    5, 'global',  1, 10);
m = pwAddK(m, 'v3',    5, 'global',  1, 10);
m = pwAddK(m, 'v4', 0.03, 'global', .01, 1);
m = pwAddK(m, 'v5',  0.1, 'global', .01, 1);
m = pwAddK(m, 'v6',  0.1, 'global', .01, 1);



m = pwAddK(m, 'k1', .1, 'global',  .01, 1);
m = pwAddK(m, 'k2', .1, 'global',  .01, 1);
m = pwAddK(m, 'k3', .1, 'global',  .01, 1);
m = pwAddK(m, 'k4', .1, 'global',  .01, 1);
m = pwAddK(m, 'k5',  1, 'global',  .1, 10);
m = pwAddK(m, 'k6', 10, 'global',  1, 100);

% coupling of the intermediate feedback component to g1p dephosphorylation
m = pwAddK(m, 'alpha', 1, 'global',  0, 1);

% coupling of the intermediate feedback component to g2p dephosphorylation
m = pwAddK(m, 'beta' , 0, 'global',  0, 1);


%% R - Reactions
% m = pwAddR(m, *ID, *reactants, *products, *modifiers, *type, *options, *rateSignature, *parameters, description, name, fast, compartments, parameterTrunks, designerPropsR, stoichiometry, reversible)
%m = pwAddR(m, 'R1', {'g1'}, {'g1p'}, {'EGF'}, 'C' , [] , 'k1*m1*r1/(k2 + r1)', {'v1','k1'});
m = pwAddR(m, 'R1', {'g1'}, {'g1p'}, {'EGF' },          'C'   , [] , 'k1*m1*r1/(k2 + r1)'                         , {'v1','k1'});

% in order to fit the same model to WT and knockdown experimental data I'm
% scaling Xact concentration with driving parameter sirna. It's not
% entirely correct but I cannot see a better way of doing it now.
m = pwAddR(m, 'R2', {'g1p'}, {'g1'}, {'Xact', 'sirna'}, 'C'   , [] , 'k1*r1/(k2 + r1) * (1 + k3 * (m1*m2 - 1))'   , {'v2','k2', 'alpha'});
m = pwAddR(m, 'R3', {'g2'}, {'g2p'}, {'g1p' },          'C'   , [] , 'k1*m1*r1/(k2 + r1)'                         , {'v3','k3'});
m = pwAddR(m, 'R4', {'g2p'}, {'g2'}, {'Xact', 'sirna'}, 'C'   , [] , 'k1*r1/(k2 + r1) * (1 + k3 * (m1*m2 - 1))'   , {'v4','k4', 'beta'});
m = pwAddR(m, 'R5', {'X'}, {'Xact'}, {'g2p'},           'C'   , [] , 'k1*m1*r1 /(k2 + r1) '                       , {'v5','k5'});
m = pwAddR(m, 'R6', {'Xact'}, {'X'}, {     },           'C'   , [] , 'k1*r1/(k2 + r1)'                            , {'v6','k6'});


%% C - Compartments
% m = pwAddC(m, *ID, *size, outside, spatialDim, name, unit, constant, designerProps, classname, description)
m = pwAddC(m, 'cytosol', 1);


%% U - Driving inputs
% m = pwAddU(m, *ID, *uType, *uTimes, *uValues, compartment, name, description, u2Values, alternativeIDs, designerProps, classname, referenceXID, unit, uFormula)
m = pwAddU(m, 'EGF'    , 'steps', [-1 0]  , [0 0.1]  , 'cytosol');

% this parameter scales the level of X protein to mimick knockdown
% experiment
m = pwAddU(m, 'sirna'  , 'steps', [-1 0]  , [0 1]  , 'cytosol');


%% A - Algebraic equations
% m = pwAddA(m, *ID, *rhs, description, name, designerProps, classname, targetType, compartment, valueType)
%m = pwAddA(m, 'Xnorm', 'X * sirna', [], [], [], [], 'species', 'cytosol');


%% Y - Observables
% m = pwAddY(m, *ID, *rhs, errorModelRhs, noiseType, unit, name, description, alternativeIDs, designerProps, classname)
m = pwAddY(m, 'ag2p', 'scale_ag2p * g2p');


%% S - Observation parameters
% m = pwAddS(m, *ID, *value, fitSetting, minValue, maxValue, unit, name, description, usedInTimeTransformation)
m = pwAddS(m, 'scale_ag2p', 1, 'global', .1, 10);
%m = pwAddS(m, 'timeshift' ,  1e-09,  'local',       0,     5);



