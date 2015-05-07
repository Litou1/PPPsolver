clear;
close all;
addpath('./funLib/')

% FD -- frequency domain
% TD -- time domain
Domain='TD';

% SolverOptions=0 : L only
% SolverOptions=1 : L and C
% SolverOptions=2 : L R C G, enforce to this for TD
% SolverOptions=3 : R only
SolverOptions=2;

% set ESLflag to 1 to include via inductance
ESLflag=1;

% set viewMesh to 1 to visualize the meshing
viewMesh=0;

% set roundViaFlag to use rounded via model
roundViaFlag=1;

f1='tc';v1=0.034;% copper thickness mm
f2='sigma';v2=5.96e7;% conductivity S/m
f3='er';v3=4.3;% define the dielectric constant
f4='tan';v4=0.025;% loss tangent
f5='planeSizeX';v5=50;% board size--X/mm
f6='planeSizeY';v6=50;% board size--Y/mm
f7='h';v7=0.6;% board separation--h/mm
f8='viaSize';v8=0.25;% via size, here, 0.5mmx0.5mm
f9='short';v9=[25 12.5];% shorting via location
f10='source';v10=[25 37.5];% power via location

% putting all dummy variables in a struct, matprty
matprty=struct(f1,v1,f2,v2,f3,v3,f4,v4,f5,v5,f6,v6,f7,v7,f8,v8,f9,v9,f10,v10);

% mesh settings
% k -- number of nodes along a via side length
% subMeshLevel -- 0 for uniform mesh
% meshT
% meshT=1, submesh region = 3x3 /uniform mesh length
% meshT=2, submesh region = 5x5 /uniform mesh length
f1='k';v1=2;
f2='subMeshLevel';v2=4;
f3='meshT';v3=1;
meshSettings=struct(f1,v1,f2,v2,f3,v3);

%% time domain example
% triangle pulses for current excitation
% unit step for voltage excitation
sourceType='current'; 
tr= 1;%ns 
tEnd= 30;%ns
width= 6;%s pulse width, make sure this is a divisor of tEnd 
Nsteps= 100;

f1='ESR';v1=0.1e-3;%kohm, 0.1ohm
f2='ESL';v2=0.5e-3;%uH 0.5nH
f3='Cd';v3=1e6;%decaps, Cd=1uF=1e6pF
f4='Rs';v4=1e-3;% voltage source series resistance, 1ohm=1e-3kohm;
f5='sourceType';v5=sourceType;
f6='riseTime';v6=tr;
f7='tEnd';v7=tEnd;
f8='Nsteps';v8=Nsteps;
f9='width';v9=width;
TDsettings=struct(f1,v1,f2,v2,f3,v3,f4,v4,f5,v5,f6,v6,f7,v7,f8,v8,f9,v9);

[unkX,~]=PPPsolver(Domain,SolverOptions,TDsettings,ESLflag,roundViaFlag,...
    viewMesh,matprty,meshSettings);
[Vin_t,Is_t]=TD_unKx_To_portVoltage_Current(sourceType,unkX,...
    size(matprty.short,2)+size(matprty.source,2));
TDplot_util(Vin_t,Is_t,linspace(0,tEnd*1e-9,Nsteps))
animation % run animation of J
%% frequency domain example
% FD: return |Z| and arg(Z) vectors
% TD: return Vt and It vectors
% can add f -- frequency vector (Ghz) as the addtional frequency vector in the
% end for frequency response

f=linspace(0.01,5,100);
[Z,arg]=PPPsolver(Domain,SolverOptions,TDsettings,ESLflag,roundViaFlag,...
    viewMesh,matprty,meshSettings,f);
Jplot_util % plot current density