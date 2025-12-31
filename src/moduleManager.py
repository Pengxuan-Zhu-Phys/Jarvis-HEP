#!/usr/bin/env python3 
from __future__ import annotations
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading
from modulePool import ModulePool
from copy import deepcopy
class ModuleManager:
    _instance = None
    _lock = threading.Lock()

    def __new__(cls, *args, **kwargs):
        with cls._lock:
            if cls._instance is None:
                cls._instance = super(ModuleManager, cls).__new__(cls)
                cls._instance.module_pools = {}
        return cls._instance

    # ModuleManager的其他方法

    def set_config(self, config):
        self.config = config

    def set_logger(self, logger):
        self.logger = logger
    
    def set_max_worker(self, nmax):
        self._max_worker = nmax

    def set_likelihood(self):
        from Module.likelihood import LogLikelihood
        self.loglikelihood = LogLikelihood(self.config['Sampling']['LogLikelihood'])

    def set_funcs(self, funcs):
        from inner_func import update_funcs
        self._funcs = update_funcs(funcs)

    @property
    def database(self):
        return self._database

    @property
    def funcs(self):
        return self._funcs

    @property
    def max_worker(self):
        return self._max_worker

    def execute_workflow(self, sample_info):
        """Execute the pre-configured workflow, updating the observables dictionary.

        Args:
            observables (dict): The dictionary of observables.
            uuid (str): The unique identifier for the sample, used for logging or other purposes.

        Returns:
            float: The calculated likelihood value.
        """
        # Execute according to the workflow's layer sequence
        self.logger.info(f"Start execute workflow for sample -> {sample_info['uuid']}")
        observables = deepcopy(sample_info['observables'])
        
        for layer in sorted(self.workflow.keys()):
            module_names = self.workflow[layer]
            with ThreadPoolExecutor(max_workers=len(module_names)) as executor:
                # Submit all modules in the current layer for parallel execution, and update observables
                future_to_module = {executor.submit(self.execute_module, module_name, observables.copy(), sample_info): module_name for module_name in module_names}
                for future in as_completed(future_to_module):
                    try:
                        # Get the returned updated observables
                        updated_observables = future.result()
                        # Merge the update results here, paying attention to concurrent update conflicts
                        observables.update(updated_observables)
                    except Exception as exc:
                        module_name = future_to_module[future]
                        self.logger.error(f'Module {module_name} generated an exception: {exc}')
                        
            passed = self.nuisance_check(observables, sample_info)
            if not passed: 
                return 1. 
        # After all modules have executed, calculate likelihood based on the final observables
        if self.config['Sampling'].get('LogLikelihood', False):
            observables = self.calculate_likelihood(observables, sample_info)
            sample_info['observables'] = observables
            self.database.add_data(observables)
            # print("Manager Line 82 ->", observables['LogL'])
            return observables['LogL']
        else: 
            sample_info['observables'] = observables 
            # print(observables)
            self.database.add_data(observables)
            return 1.

    def nuisance_check(self, observables: dict, sample_info: dict) -> bool:
        """Check nuisance passconditions and decide whether to early-stop workflow.

        Contract:
        - Uses dict-input callable interface: entry['expr'](observables) -> bool
        - Uses deps checker: entry['checker'](observables.keys()) -> (ok, missing)
        - Only evaluates a condition when all deps are available.
        - If any evaluated condition returns False, triggers early stop.

        Returns:
            (stop: bool, reason: str)
        """
        
        pcdt = True
        if not self._with_nuisance: 
            return True
        
        if self.nuisance_loglikelihoods: 
            for ll, loglike in self.nuisance_loglikelihoods.items(): 
                if loglike['checker'](observables.keys())[0]:
                    value = loglike['expr'](observables)
                    observables[ll] = value 
                    sample_info['nuisance']['active']['LogLikelihoods'][ll] = value 
                    sample_info['logger'].info("Evaluating Nuisance Loglikelihood -> \n\tname \t-> {}\n\texpr \t-> {}\n\tOutput \t-> {}".format(ll, loglike['expression'], value))
                
        if self.nuisance_passconditions: 
            for pc, passcond in self.nuisance_passconditions.items(): 
                if passcond['checker'](observables.keys())[0]: 
                    value = passcond['expr'](observables)
                    observables[pc] = value 
                    sample_info['nuisance']['active']['PassConditions'][pc] = value
                    sample_info['logger'].info("Evaluating Nuisance PassConditions -> \n\tname \t-> {}\n\texpr \t-> {}\n\tOutput \t-> {}".format(pc, passcond['expression'], value))
                    pcdt = pcdt and value 
            sample_info['nuisance']['active']['pass'] = pcdt

        if pcdt: 
            observables['uuid'] = sample_info['uuid']
            
        observables['NAttempt'] = sample_info['nuisance']['NAttempt']
        sample_info['observables'] = observables
        
        return pcdt



    def execute_module(self, module_name, observables, sample_info):
        """Execute a single module's computation task, updating the observables dictionary.

        Args:
            module_name (str): The name of the module.
            observables (dict): The dictionary of observables.
            uuid (str): The unique identifier for the sample.

        Returns:
            dict: The updated observables dictionary.
        """
        module_pool = self.module_pools.get(module_name)
        if not module_pool:
            self.logger.error(f"Module pool for '{module_name}' not found.")
            return {}

        updated_observables = module_pool.execute(observables, sample_info)
        return updated_observables

    def calculate_likelihood(self, observables, sample_info):
        """Calculate the likelihood based on the updated observables.

        Args:
            observables (dict): The final observables dictionary.

        Returns:
            float: The calculated likelihood value.
        """
        # Calculate likelihood based on observables, this is pseudocode for the calculation logic

        logl = self.loglikelihood.calculate(observables, sample_info)
        observables.update(logl)
        return observables

    def add_module_pool(self, module):
        if module.name not in self.module_pools:
            self.logger.warning(f"Manager adding ModulePool {module.name}. ")
            self.module_pools[module.name] = ModulePool(module, max_workers=self.max_worker)
            self.module_pools[module.name].set_funcs(self.funcs)
            self.module_pools[module.name].set_logger()
            # self.module_pools[module.name].load_installed_instances()
