#!/usr/bin/env python3 

import os, sys 
import numpy as np 
from Sampling.sampler import SamplingVirtial
from sample import Sample
from copy import deepcopy 
from uuid import uuid4
import pandas as pd 
import dill 

class Diver(SamplingVirtial):
    def __init__(self) -> None:
        super().__init__()
        self.load_schema_file()
        self.method = "Diver"

    def set_logger(self, logger) -> None:
        super().set_logger(logger)
        self.logger.warning("Sampling method initializaing ...")

    def set_factory(self, factory) -> None:
        self.factory = factory
        self.logger.warning("WorkerFactory is ready for Dynesty sampler")

    def load_schema_file(self):
        self.schema = self.path['DiverSchema']

    def set_config(self, config_info) -> None:
        self.config = config_info
        self.init_generator()

    def init_generator(self) -> None:
        self.load_variable()
        self._D = len(self.vars)
        self.set_parameters()
        self.param_assign()

    def set_parameters(self):
        self._mode = self.config['Sampling']["run"].get('mode', 'current')
        if self._mode == "jDE":
            # valid_bndry_options = ['not enforced', 'brick wall', 'random re-initialization', 'reflection']
            self._bndry = self.config['Sampling']["run"].get('bndry', 'brick wall')
            self._removeDuplicates = self.config['Sampling']["run"].get('removeDuplicates', False)
            self._doBayesian = self.config['Sampling']["run"].get('doBayesian', False)

            self._maxNodePop = self.config['Sampling']["run"].get('maxNodePop', 100)
            self._Ztolerance = self.config['Sampling']["run"].get('Ztolerance', 1e-3)

            self._init_population_strategy = self.config['Sampling']["run"].get('init_population_strategy', 0)
            self._max_initialisation_attempts = self.config['Sampling']["run"].get('max_initialisation_attempts', 5)
            self._max_acceptable_value = self.config['Sampling']["run"].get('max_acceptable_value', float('inf'))

            # 设置种群大小和最大代数等其他参数
            self._strategy = "self-adaptive rand/1/bin (jDE)"
            self._Fsize = 1 
            self._F  = [0.0]
            self._Cr = None
            self._lambda = 0.0 
            self._expon = False
            self._maxciv = self.config['Sampling']["run"].get('maxciv', 1)
            self._maxgen = self.config['Sampling']["run"].get('maxgen', 500)
            self._convthresh = self.config['Sampling']["run"].get('convthresh', 1e-3)
            self._convsteps = self.config['Sampling']["run"].get("convsteps", 10)

            self._NP = self.config['Sampling']["run"].get('NP', 100)

            self.seed = self.config['Sampling']["run"].get('seed', None)
            if self.seed is not None:
                np.random.seed(self.seed)

        if self._mode == "lambdajDE":
            """
            参数分配函数 - 适用于 lambdajDE 模式
            """
            self._strategy = 'self-adaptive rand-to-best/1/bin (lambdajDE)'
            self._Fsize = 1
            self._F = [0.0]
            self._Cr = None
            self._lambda = self.config["Sampling"]["run"].get("lambda", 0.0)
            self._expon = False
            self._maxciv = self.config['Sampling']["run"].get('maxciv', 1)
            self._maxgen = self.config['Sampling']["run"].get('maxgen', 500)
            self._convthresh = self.config['Sampling']["run"].get('convthresh', 1e-3)
            self._convsteps = self.config['Sampling']["run"].get("convsteps", 10)

        elif self._mode == "current":
            """
            参数分配函数 - 适用于 current 模式
            """
            # 基础参数初始化
            self._maxciv    = self.config['Sampling']["run"].get('maxciv', 1)
            self._maxgen    = self.config['Sampling']["run"].get('maxgen', 500)
            self._NP        = self.config["Sampling"]["run"].get("NP", 100)
            self._F         = self.config["Sampling"]["run"].get("F",  [0.7])
            self._Cr        = self.config["Sampling"]["run"].get("Cr", 0.9)
            self._convthresh = self.config['Sampling']["run"].get('convthresh', 1e-3)
            self._convsteps = self.config['Sampling']["run"].get("convsteps", 10)



        self.population = np.random.uniform(0, 1, (self._NP, self.config['Sampling']['run'].get('dim', 2)))
        self.fitness = np.full(self._NP, float('inf'))

        # 随机种子




    















        if self._mode == "lambdajDE":
            self._lambda = self.config['Sampling']["run"].get('lambda', 0.0) # default rand/1/bin 

        self._expon = self.config['Sampling']["run"].get('expon', False)

        # valid_bndry_options = ['not enforced', 'brick wall', 'random re-initialization', 'reflection']
        self._bndry = self.config['Sampling']["run"].get('bndry', 'brick wall')
        self._convthresh = self.config['Sampling']["run"].get('convthresh', 1e-6)
        self._removeDuplicates = self.config['Sampling']["run"].get('removeDuplicates', False)
        self._doBayesian = self.config['Sampling']["run"].get('doBayesian', False)

        self._maxNodePop = self.config['Sampling']["run"].get('maxNodePop', 100)
        self._Ztolerance = self.config['Sampling']["run"].get('Ztolerance', 1e-3)

        self._init_population_strategy = self.config['Sampling']["run"].get('init_population_strategy', 0)
        self._max_initialisation_attempts = self.config['Sampling']["run"].get('max_initialisation_attempts', 5)
        self._max_acceptable_value = self.config['Sampling']["run"].get('max_acceptable_value', float('inf'))

        self._Cr = self.config['Sampling']['run'].get("Cr", 0.7)
        self._F  = self.config['Sampling']['run'].get("F", 0.8)

        # 设置种群大小和最大代数等其他参数
        self._NP = self.config['Sampling']["run"].get('NP', 100)
        self._maxgen = self.config['Sampling']["run"].get('maxgen', 500)
        self._maxciv = self.config['Sampling']["run"].get('maxciv', 1)
        self.population = np.random.uniform(0, 1, (self._NP, self.config['Sampling']['run'].get('dim', 2)))
        self.fitness = np.full(self._NP, float('inf'))

        # 随机种子
        self.seed = self.config['Sampling']["run"].get('seed', None)
        if self.seed is not None:
            np.random.seed(self.seed)

        self.logger.warning("Diver Sampling\nCopyright Elinore Roebber and Pat Scott")
        self.logger.warning("If you use Diver Sampling method in your work, please cite \n\t1.) arXiv:1705.07959. 2.) arXiv:2101.04525. ")
        self.check_parameter()

    def check_parameter(self, F=None, Cr=None, current=False, expon=False, lambda_val=None):
        """
        在 jDE 或 lambdajDE 模式下设置差分进化参数，并处理 F, Cr, lambda, current, expon 等参数的警告。
        """
        # 仅当 mode 设置为 'jDE' 或 'lambdajDE' 时执行
        if self._mode == 'jDE':
            if self._F is not None:
                self.logger.warning('WARNING: value set for F not used during jDE run')
            if self._Cr is not None:
                self.logger.warning('WARNING: value set for Cr not used during jDE run')
            if self._expon:
                self.logger.warning('WARNING: jDE uses binary crossover. Value set for expon will be ignored.')

            # jDE 模式下 F 的大小设为 1，并初始化为 0
            self._Fsize = 1
            self._F = np.zeros(self._Fsize)

        if self._mode == 'lambdajDE':
            self._lambda = 0.0  
            if lambda_val is not None and self._verbose >= 1:
                print('WARNING: value set for lambda not used during lambdajDE run')
        else:
            if lambda_val is not None:
                if self._verbose >= 1:
                    if lambda_val < 0.0:
                        print('WARNING: lambda < 0. DE may not converge properly.')
                    if lambda_val > 1.0:
                        print('WARNING: lambda > 1. DE may not converge properly.')
                self._lambda = lambda_val
            else:
                self._lambda = 0.0  # 默认值为 0，即 rand/1/bin

            # 输出突变策略的信息
            if self._mode == 'lambdajDE':
                self._DEstrategy = 'self-adaptive rand-to-best/1/bin (lambdajDE)'
            elif self._lambda == 0.0:
                self._DEstrategy = 'self-adaptive rand/1/bin (jDE)'
            elif self._lambda == 1.0:
                self._DEstrategy = 'self-adaptive best/1/bin (jDE with additional fixed lambda)'
            else:
                self._DEstrategy = 'self-adaptive rand-to-best/1/bin (jDE with additional fixed lambda)'

            """
            在非 jDE 模式下设置差分进化 (DE) 参数，包括 F, Cr, lambda, current, expon 等，并处理相关警告。
            """
        else:
            # 设置 F 参数
            if F is not None:
                if np.any(F <= 0.0):
                    if self._verbose >= 1:
                        self.logger.warning('WARNING: some elements of F are 0 or negative. DE may not converge properly.')
                elif np.any(F >= 1.0):
                    if self._verbose >= 1:
                        print('WARNING: some elements of F are 1 or greater. DE may not converge properly.')
                self._Fsize = len(F)
                self._F = F
            else:
                # 默认 F 值为 0.7
                self._Fsize = 1
                self._F = np.array([0.7])

            # 设置 Cr 参数
            if Cr is not None:
                if Cr < 0.0:
                    if self._verbose >= 1:
                        print('WARNING: Cr < 0. Using Cr = 0.')
                    self._Cr = 0.0
                elif Cr > 1.0:
                    if self._verbose >= 1:
                        print('WARNING: Cr > 1. Using Cr = 1.')
                    self._Cr = 1.0
                else:
                    self._Cr = Cr
            else:
                self._Cr = 0.9  # 默认 Cr 值为 0.9

            # 设置 lambda 参数
            if lambda_val is not None:
                if self._verbose >= 1:
                    if lambda_val < 0.0:
                        print('WARNING: lambda < 0. DE may not converge properly.')
                    if lambda_val > 1.0:
                        print('WARNING: lambda > 1. DE may not converge properly.')
                self._lambda = lambda_val
            else:
                self._lambda = 0.0  # 默认 lambda 值为 0

            # 设置突变策略 current 和 expon
            self._current = current
            self._expon = expon

            # 确定突变策略
            if self._lambda == 0.0:
                if self._current:
                    DEstrategy = 'current/'
                else:
                    DEstrategy = 'rand/'
            elif self._lambda == 1.0:
                DEstrategy = 'best/'
            else:
                if self._current:
                    DEstrategy = 'current-to-best/'
                else:
                    DEstrategy = 'rand-to-best/'

            # 添加 Fsize 到策略名称
            DEstrategy += str(self._Fsize)

            # 确定交叉策略（exponential 或 binary）
            if self._expon:
                DEstrategy += '/exp'
            else:
                DEstrategy += '/bin'

            self._DEstrategy = DEstrategy


        if self._mode == "lambdajDE":
            self._removeDuplicates = True
        else: 
            self._removeDuplicates = False

    def param_assign(self):
        bestFitParams = kwargs.get('bestFitParams', None)
        bestFitDerived = kwargs.get('bestFitDerived', None)
        discrete = kwargs.get('discrete', [])
        partitionDiscrete = kwargs.get('partitionDiscrete', False)
        maxciv = kwargs.get('maxciv', 1)
        maxgen = kwargs.get('maxgen', 300)
        NP = kwargs.get('NP', max(10 * len(lowerbounds), 7))  # Conservative default
        F = kwargs.get('F', [0.7])
        Cr = kwargs.get('Cr', 0.9)
        lambda_ = kwargs.get('lambda', 0.0)
        current = kwargs.get('current', False)
        expon = kwargs.get('expon', False)
        bndry = kwargs.get('bndry', 1)
        jDE = kwargs.get('jDE', True)
        lambdajDE = kwargs.get('lambdajDE', True)
        convthresh = kwargs.get('convthresh', 1e-3)
        convsteps = kwargs.get('convsteps', 10)
        removeDuplicates = kwargs.get('removeDuplicates', False)
        doBayesian = kwargs.get('doBayesian', False)
        maxNodePop = kwargs.get('maxNodePop', None)
        Ztolerance = kwargs.get('Ztolerance', 0.01)
        savecount = kwargs.get('savecount', 1)
        disableIO = kwargs.get('disableIO', False)
        outputRaw = kwargs.get('outputRaw', True)
        outputSam = kwargs.get('outputSam', True)
        init_population_strategy = kwargs.get('init_population_strategy', 0)
        discard_unfit_points = kwargs.get('discard_unfit_points', False)
        max_initialisation_attempts = kwargs.get('max_initialisation_attempts', 10000)
        max_acceptable_value = kwargs.get('max_acceptable_value', 1e6)
        seed = kwargs.get('seed', -1)
    

    # Initialize run_params dictionary
    run_params.update({
        'D': D,
        'DE': {},
        'maxNodePop': maxNodePop if doBayesian else None,
        'tol': Ztolerance if doBayesian else None,
        'numciv': maxciv,
        'savefreq': savecount,
        'calcZ': doBayesian,
        'init_population_strategy': init_population_strategy,
        'discard_unfit_points': discard_unfit_points,
        'max_initialisation_attempts': max_initialisation_attempts,
        'max_acceptable_value': max_acceptable_value,
        'disableIO': disableIO,
        'outputRaw': outputRaw and not disableIO,
        'context': kwargs.get('context', None),
        'convthresh': set_if_positive_real('convthresh', convthresh, 1e-3),
        'convsteps': set_if_positive_int('convsteps', convsteps, 10),
    })

    # Set DE strategy
    if self._mode == "jDE":
        if lambdajDE:
            run_params['DE']['strategy'] = 'self-adaptive rand-to-best/1/bin (lambdajDE)'
            run_params['DE']['lambdajDE'] = True
        else:
            run_params['DE']['strategy'] = 'self-adaptive rand/1/bin (jDE)'
            run_params['DE']['lambdajDE'] = False
        run_params['DE']['jDE'] = True
        run_params['DE']['current'] = False
        run_params['DE']['expon'] = False
        run_params['DE']['Fsize'] = 1
        run_params['DE']['F'] = [0.0]  # Dummy value for F
    else:
        run_params['DE']['Fsize'] = len(F)
        run_params['DE']['F'] = F
        run_params['DE']['Cr'] = Cr
        run_params['DE']['lambda'] = lambda_
        run_params['DE']['current'] = current
        run_params['DE']['expon'] = expon
        run_params['DE']['strategy'] = f"rand/1/bin" if lambda_ == 0 else f"rand-to-best/{len(F)}/bin"

    # Check and adjust NP (population size)
    run_params['DE']['NP'] = set_if_positive_int('NP', NP, max(10 * D, 7))
    if run_params['DE']['NP'] < (2 * len(F) + 3):
        run_params['DE']['NP'] = 2 * len(F) + 3

    # Handle partitionDiscrete and discrete variables
    run_params['discrete'] = discrete
    run_params['partitionDiscrete'] = partitionDiscrete and len(discrete) > 0
    if run_params['partitionDiscrete']:
        # Logic to handle partitioning (omitted for simplicity)
        pass

    # Print parameters if verbose
    if verbose >= 1:
        print(f"Selected DE strategy: {run_params['DE']['strategy']}")
        print(f"NP = {run_params['DE']['NP']}, F = {run_params['DE']['F']}, Cr = {run_params['DE'].get('Cr')}, lambda = {run_params['DE'].get('lambda')}")
        if run_params['partitionDiscrete']:
            print(f"Discrete dimensions: {run_params['discrete']}")
        print(f"Max generations: {maxgen}, Max civilizations: {maxciv}")
        if doBayesian:
            print(f"Bayesian evidence calculation enabled. MaxNodePop = {maxNodePop}, Z tolerance = {Ztolerance}")

    return run_params


    async def optimize(self, func, path="de_save.pkl"):
        """执行差分进化算法的核心逻辑"""
        for civ in range(self._max_civ):
            print(f"Starting civilization {civ}")
            for gen in range(self._max_gen):
                self._gen = gen
                await self.evolve_population(func)
                print(f"Generation {gen} complete. Best fitness: {self._best_fit}")

                # 定期保存进度
                if gen % self._save_interval == 0:
                    self.save_state(path)

                # 检查是否满足退出条件
                if self.check_convergence():
                    print("Convergence reached.")
                    break
            self._civ += 1

    async def evolve_population(self, func):
        """对种群进行突变、交叉、选择，更新种群"""
        tasks = []
        for individual in self._population:
            tasks.append(self.mutate_and_select(individual, func))
        results = await asyncio.gather(*tasks)
        self._population = np.array(results)

    async def mutate_and_select(self, individual, func):
        """突变和选择过程"""
        trial_vector = self.mutate(individual)
        fitness = func(trial_vector)
        if fitness < self._best_fit:
            self._best_fit = fitness
            self._best_params = trial_vector
        return trial_vector

    def check_convergence(self):
        """检查是否满足收敛条件"""
        return self._best_fit < 1e-6

    def save_state(self, path):
        """使用dill库保存当前的算法状态"""
        with open(path, 'wb') as f:
            dill.dump({
                'population': self._population,
                'best_fit': self._best_fit,
                'best_params': self._best_params,
                'gen': self._gen,
                'civ': self._civ
            }, f)
        print(f"State saved to {path}")

    def load_state(self, path):
        """使用dill库恢复算法状态"""
        with open(path, 'rb') as f:
            state = dill.load(f)
            self._population = state['population']
            self._best_fit = state['best_fit']
            self._best_params = state['best_params']
            self._gen = state['gen']
            self._civ = state['civ']
        print(f"State loaded from {path}")


def initialize_run_params(self, seed=None):
    """
    初始化运行参数并分配种群向量和自适应 DE 内存。
    """
    # 初始化随机数生成器
    self._seed = seed if seed is not None else np.random.randint(1, 1e6)
    np.random.seed(self._seed)

    # 分配种群向量的存储空间
    self._X = np.zeros((self._NP, self._D))  # 主种群向量
    self._Xnew = np.zeros((self._NP, self._D))  # 新种群向量
    self._Xsub = np.zeros((self._NP, self._D))  # 子种群向量（用于离散参数）
    self._FjDE = None
    self._CrjDE = None
    self._lambdajDE = None

    # 如果使用自适应 DE (jDE)，分配额外的内存
    if self._mode == 'jDE':
        self._FjDE = np.zeros(self._NP)  # 自适应 F
        self._CrjDE = np.zeros(self._NP)  # 自适应 Cr
        if self._mode == 'lambdajDE':
            self._lambdajDE = np.zeros(self._NP)  # 自适应 lambda

    # 设置收敛标准（默认是 meanimprovement，可以根据需要更改）
    self._convergence_criterion = 'meanimprovement'

    # 检查贝叶斯推理是否启用，如果启用则需要 prior 函数
    if self._doBayesian and self._prior is None:
        raise ValueError("Error: evidence calculation requested without specifying a prior.")

    # 如果禁用了 IO，则不能进行证据计算
    if self._doBayesian and self._disableIO:
        raise ValueError("Error: evidence calculation is not possible with IO disabled.")




    # Core of Diver algorithm 
    def __iter__(self):
        """make the diver as an iterator object"""
        return self

    def __next__(self):
        """Reture an individual and its fitness value for each iteration """
        # Stop iteration if all generation is done 
        if self.current_gen >= self.max_gen:
            raise StopIteration

        # Generate a new generation, if all individuals are finished
        if self.current_index >= self.population_size:
            self.evolve_population()  # 生成新一代
            self.current_index = 0
            self.current_gen += 1

        # Get the current individual and its fitness value 
        individual = self.population[self.current_index]
        fitness_value = self.fitness[self.current_index]

        # update the optimized value 
        if fitness_value < self.best_score:
            self.best_score = fitness_value
            self.best_fit = individual

        # Add new index, preparting the next individual 
        self.current_index += 1

        return individual, fitness_value, self.best_fit, self.best_score

    def evolve_population(self):
        """种群进化（突变、交叉、选择）"""
        new_population = np.copy(self.population)
        new_fitness = []

        for i in range(self.population_size):
            mutant = self.mutate(new_population)
            trial = self.crossover(mutant, new_population[i])
            fitness_value = self.func(trial)
            new_fitness.append(fitness_value)

            # 选择逻辑
            if fitness_value < self.fitness[i]:
                new_population[i] = trial

        self.population = new_population
        self.fitness = np.array(new_fitness)

    def mutate(self, population):
        """突变操作"""
        idxs = np.random.choice(len(population), 3, replace=False)
        a, b, c = population[idxs]
        F = self.kwargs.get('F', 0.8)
        return a + F * (b - c)

    def crossover(self, mutant, target):
        """交叉操作"""
        dim = len(target)
        trial = np.copy(target)
        Cr = self.kwargs.get('Cr', 0.7)
        for i in range(dim):
            if np.random.rand() < Cr:
                trial[i] = mutant[i]
        return trial
