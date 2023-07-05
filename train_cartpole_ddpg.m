clear all; clc; close all

ObservationInfo = rlNumericSpec([4 1]);
ObservationInfo.Name = "CartPole States";
ObservationInfo.Description = 'x, dx, theta, dtheta';

actInfoCont = rlNumericSpec([1 1], 'LowerLimit', -1, 'UpperLimit', 1);

type myResetFunction.m

type myStepFunction2.m


% Acceleration due to gravity in m/s^2
envConstants.Gravity = 9.8;
% Mass of the cart
envConstants.MassCart = 1.0;
% Mass of the pole
envConstants.MassPole = 0.1;
% Half the length of the pole
envConstants.Length = 0.5;
% Max force the input can apply
envConstants.MaxForce = 1;
% Sample time
envConstants.Ts = 0.02;
% Angle at which to fail the episode
envConstants.ThetaThresholdRadians = 12 * pi/180;
% Distance at which to fail the episode
envConstants.XThreshold = 2.4;
% Reward each time step the cart-pole is balanced
envConstants.RewardForNotFalling = 1;
% Penalty when the cart-pole fails to balance
envConstants.PenaltyForFalling = -1;

StepHandle = @(Action,LoggedSignals) myStepFunction2(Action,LoggedSignals,envConstants);

ResetHandle = @() myResetFunction;


env = rlFunctionEnv(ObservationInfo,actInfoCont,StepHandle,ResetHandle);
obsInfo = getObservationInfo(env);
actInfo = getActionInfo(env);
Ts = 0.01;
Tf = 10;

% Define state path
statePath = [
    featureInputLayer( ...
        obsInfo.Dimension(1), ...
        Name="obsPathInputLayer")
    fullyConnectedLayer(40)
    reluLayer
    fullyConnectedLayer(30,Name="spOutLayer")
    ];

% Define action path
actionPath = [
    featureInputLayer( ...
        actInfo.Dimension(1), ...
        Name="actPathInputLayer")
    fullyConnectedLayer(30, ...
        Name="apOutLayer", ...
        BiasLearnRateFactor=0)
    ];

% Define common path
commonPath = [
    additionLayer(2,Name="add")
    reluLayer
    fullyConnectedLayer(1)
    ];

% Create layergraph, add layers and connect them
criticNetwork = layerGraph();
criticNetwork = addLayers(criticNetwork,statePath);
criticNetwork = addLayers(criticNetwork,actionPath);
criticNetwork = addLayers(criticNetwork,commonPath);
criticNetwork = connectLayers(criticNetwork,"spOutLayer","add/in1");
criticNetwork = connectLayers(criticNetwork,"apOutLayer","add/in2");

criticNetwork = dlnetwork(criticNetwork);
summary(criticNetwork)
plot(criticNetwork)

critic = rlQValueFunction(criticNetwork, ...
    obsInfo,actInfoCont, ...
    ObservationInputNames="obsPathInputLayer", ...
    ActionInputNames="actPathInputLayer");

actorNetwork = [
    featureInputLayer(obsInfo.Dimension(1))
    fullyConnectedLayer(40)
    reluLayer
    fullyConnectedLayer(30)
    reluLayer
    fullyConnectedLayer(1)
    tanhLayer
    scalingLayer(Scale=max(actInfoCont.UpperLimit))
    ];


actorNetwork = dlnetwork(actorNetwork);
summary(actorNetwork)

actor = rlContinuousDeterministicActor(actorNetwork,obsInfo,actInfoCont);

criticOpts = rlOptimizerOptions(LearnRate=1e-03,GradientThreshold=1);
actorOpts = rlOptimizerOptions(LearnRate=1e-04,GradientThreshold=1);

agentOpts = rlDDPGAgentOptions(...
    SampleTime=Ts,...
    CriticOptimizerOptions=criticOpts,...
    ActorOptimizerOptions=actorOpts,...
    ExperienceBufferLength=1e6,...
    DiscountFactor=0.99,...
    MiniBatchSize=128);

agentOpts.NoiseOptions.Variance = 0.6;
agentOpts.NoiseOptions.VarianceDecayRate = 1e-5;
agent = rlDDPGAgent(actor,critic,agentOpts);
maxepisodes = 5000;
maxsteps = ceil(Tf/Ts);
trainOpts = rlTrainingOptions(...
    MaxEpisodes=maxepisodes,...
    MaxStepsPerEpisode=maxsteps,...
    ScoreAveragingWindowLength=5,...
    Verbose=false,...
    Plots="training-progress",...
    StopTrainingCriteria="AverageReward",...
    StopTrainingValue=400,...
    SaveAgentCriteria="EpisodeReward",...
    SaveAgentValue=250);

    trainingStats = train(agent,env,trainOpts);

simOptions = rlSimulationOptions(MaxSteps=500);
experience = sim(env,agent,simOptions);

obs = squeeze(experience.Observation.CartPoleStates.Data);
figure(1)
plot(obs(3, :)), grid on, xlabel('time [s]'), ylabel('\theta [deg]')

obs = squeeze(experience.Reward);
figure(2)
plot(obs(1, :)), grid on, xlabel('time'), ylabel('Reward')

obs = squeeze(experience.Action.act1);
figure(3)
plot(obs(1, :)), grid on, xlabel('action'), ylabel('value')


