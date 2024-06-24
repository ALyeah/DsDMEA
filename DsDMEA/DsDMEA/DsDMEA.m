classdef DsDMEA< ALGORITHM
% <multi> <real> <multimodal>


    methods
        function main(Algorithm,Problem)
            %% Generate the sampling points and random population
            Population = Problem.Initialization();
            [CR,F] = Algorithm.ParameterSet(0.9,0.5);
            ReferencePointobj = Population;
            ReferencePointdec = Population;

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = randi(Problem.N,1,2*Problem.N);
                Offspring  = OperatorDE(Population,Population(MatingPool(1:end/2)),Population(MatingPool(end/2+1:end)),{CR,F,0,0});
                ReferencePointobj = obj([ReferencePointobj,Offspring],5.*Problem.N);
                ReferencePointdec = dec([ReferencePointdec,Offspring],5.*Problem.N);
                Population = EnvironmentalSelection([Population,Offspring],ReferencePointobj,ReferencePointdec,Problem.N);
            end
        end
    end
end