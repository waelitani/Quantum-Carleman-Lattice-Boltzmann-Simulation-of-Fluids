%% Clear your desk to focus
clear;
clc;
close all;

markerstyles = {'_','.','o','+','x','*','|'};

maxQC = 4;
minQC = 2;
Dt = 1;
tend = 10;
counter = 0;

%% Set parameters
nbDims = 1;
nbVars = 3;
tarray = 0:Dt:tend;

taurange = 3:1:6;
taurange = 10.^taurange;

for taucount = 1:1:length(taurange)
%% Define the function to solve
tau = taurange(taucount);
intau = -1/tau;

%D1Q3

funcToSolve{1} = @(f) 0.5*intau*(f{1}+f{3}-(f{1}-f{3})^2-(1/3));
funcToSolve{2} = @(f) intau*(f{2}+(f{1}-f{3})^2-(2/3));
funcToSolve{3} = @(f) 0.5*intau*(f{1}+f{3}-(f{1}-f{3})^2-(1/3));

funcDivergence = (nbVars-nbDims)*intau;
funcNorm = @(t) exp(funcDivergence*t);

%% Set the initial values
uu = 0.1;
initialValues = [1/6-uu/2;2/3;1/6+uu/2];

%% Obtain classical solution
classicalSol = classicalNumericalSolution(funcToSolve,initialValues,nbVars,tend,Dt);

%% Plot classical solution
subplot(2,nbVars,1);
ylabel('Simulated Quantum Solution');
subplot(2,nbVars,nbVars+1);
ylabel('First-Order in Time Classical Solution');

for v = 1:1:nbVars    
    subplot(2,nbVars,nbVars+v);
    p = plot(tarray,classicalSol(v,:),'k','DisplayName',['tau = ',num2str(tau)]);
    p.Marker = markerstyles{taucount};
    hold on;
    xlabel('Time (sec)');
end


%% Obtain and plot the quantum solution
% Set quantum parameters
% quantumOptions(normalized,kowalski,positionmomentum,lastlevelpush,taylorEvolution,taylorOrder)
options = quantumOptions(true,false,false,0);

for qc = minQC:1:maxQC
counter = counter + 1;
hold on;

[a,ad,q,p,nn] = operators(qc);
[initialState] = initialize(a,ad,q,p,nn,initialValues,tau,qc,nbVars,options);
[evolutionOp,opa,opad] = linearembedding(a,ad,q,p,nn,funcToSolve,Dt,tau,qc,options,nbVars,funcNorm);

sol = quantumSolution(evolutionOp,opa,opad,initialState,nbVars,tend,Dt,tau,funcNorm,options);

for v = 1:1:nbVars    
    subplot(2,nbVars,v);
    hold on;
    lbltxt = [];
    
    if options.TaylorExpandEvolutionOperator == true
        lbltxt = ['qc = ',num2str(qc),' , Taylor Order = ',num2str(options.TaylorExpansionOrder)];
    else
        lbltxt = ['qc = ',num2str(qc)];
    end
    
    p = plot(tarray,sol(v,:),'DisplayName',lbltxt);
    p.Marker = markerstyles{taucount};
    if (max(max(sol)) > 1)
        ylim([0 1]);
    end
    xlim([0 tend]);
    title(['f',num2str(v)])
    xlabel('Time (sec)');
end
end
end

subplot(2,nbVars,1);
legend('Position',[0.4675 0.4775 0.1 0.03],'Orientation','Horizontal');
subplot(2,nbVars,nbVars+1);
legend('Position',[0.4675 0.4425 0.1 0.03],'Orientation','Horizontal');
sgtitle(['Repeat Collision Process in a Single D1Q3 cell | u(t = 0) = ',num2str(uu),'-> | \Delta t = ',num2str(Dt)]);

saveas(gcf,[num2str(maxQC),'qc',num2str(Dt),'dt',num2str(tend),'Trun'],'eps');
