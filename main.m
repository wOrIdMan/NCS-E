%% cec2020
clear;
close;
addpath benchmark/CEC2020/
global initial_flag
Date=date;
Date(Date=='-') = '_';
MaxRun = 25;
epsim=1e-4;
% algorithms = {'NCS_FR', 'eNCS','vNCS','MOEAD_NCS'};
algorithms = {'NCS_E'};
% for i = 1:length(algorithms)
%     addpath(sprintf('%s/',algorithms{i}));
% end
for algorithm = algorithms
    algo_name = string(algorithm);
    epoch = 10;
    % f = fopen(sprintf('res/%s_%s2020.txt',Date,algo_name),'a');
    for funcNum = [1:57]
        [par] = Cal_par(funcNum);
        D = par.n;
        lb = par.xmin;
        ub = par.xmax;
        if D<=10
            MaxFes = 1e5;
            r = 0.95;
        elseif D<=30
            MaxFes = 2*1e5;
            r = 0.98;
        elseif D<=50
            MaxFes = 4*1e5;
            r = 0.99;
        elseif D<=150
            MaxFes = 8*1e5;
            r = 0.994;
        else 
            MaxFes = 1e6;
            r = 0.994;
        end
        if funcNum==3
            r=0.991;
        end
        funcNum
        initial_flag = 0;

        
        func = @(x,funcNum) cec20_func(x,funcNum);
        bestx = nan(MaxRun,D,10);
        FeasibleObj1st = nan(1,MaxRun);
        FeasibleFes1st = nan(1,MaxRun);
        input.r = r;
        input.epoch = epoch;
        for run = 1:MaxRun
            [bestx(run,:,:), FeasibleObj1st(run), FeasibleFes1st(run)] = feval(algo_name,MaxFes,func,funcNum,lb,ub,10,D,input);
        end
        % save(sprintf('res/%s_x_%s2020_RC%d.mat',Date,algo_name,funcNum),"bestx");
        [besty,bestg,besth] = cec20_func(bestx(:,:,end),funcNum);
        Vio = calVio(MaxRun,bestg,besth,epsim);
        fprintf('%e ',besty);fprintf('\n');fprintf('%.6e ',Vio);fprintf('\n');fprintf('%f ',FeasibleObj1st);fprintf('\n');fprintf('%d ',FeasibleFes1st);fprintf('\n');
%         fprintf(f,'%e ',besty);fprintf(f,'\n');fprintf(f,'%.6e ',Vio);fprintf(f,'\n');fprintf(f,'%f ',FeasibleObj1st);fprintf(f,'\n');fprintf(f,'%d ',FeasibleFes1st);fprintf(f,'\n');
    end
%     fclose(f);

        
end