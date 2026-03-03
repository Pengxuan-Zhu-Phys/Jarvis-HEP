#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The wrapper around multiprocessing pool that can be helpful
with dynesty since it avoids some overhead that one would get
with standard pool
"""

import multiprocessing as mp
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait

__all__ = ['Pool', 'JarvisFactoryAsyncPool']


class FunctionCache:
    """
    Singleton class to cache the functions and optional arguments between calls
    """


def initializer(loglike, prior_transform, logl_args, logl_kwargs, ptform_args,
                ptform_kwargs):
    """
    Initialized function used to initialize the
    singleton object inside each worker of the pool
    """
    FunctionCache.loglike = loglike
    FunctionCache.prior_transform = prior_transform
    FunctionCache.logl_args = logl_args
    FunctionCache.logl_kwargs = logl_kwargs
    FunctionCache.ptform_args = ptform_args
    FunctionCache.ptform_kwargs = ptform_kwargs


def loglike_cache(x):
    """
    Likelihood function call
    """
    return FunctionCache.loglike(x, *FunctionCache.logl_args,
                                 **FunctionCache.logl_kwargs)


def prior_transform_cache(x):
    """
    Prior transform call
    """
    return FunctionCache.prior_transform(x, *FunctionCache.ptform_args,
                                         **FunctionCache.ptform_kwargs)


class Pool:
    """
    The multiprocessing pool wrapper class
    It is intended to be used as a context manager for dynesty sampler only.

    Parameters
    ----------
    njobs: int
        The number of multiprocessing jobs/processes
    loglike: function
        ln(likelihood) function
    prior_transform: function
        Function transforming from a unit cube to the parameter
        space of interest according to the prior
    logl_args: tuple(optional)
        The optional arguments to be added to the likelihood
        function call
    logl_kwargs: tuple(optional)
        The optional keywords to be added to the likelihood
        function call
    ptform_args: tuple(optional)
        The optional arguments to be added to the prior transform
        function call
    ptform_kwargs: tuple(optional)
        The optional keywords to be added to the prior transform
        function call

    Attributes
    ----------
    loglike: function
        ln(likelihood) function
    prior_transform: function
        Function transforming from a unit cube to the parameter
        space of interest according to the prior

    Examples
    --------
    To use the dynest pool you have to use it with the context manager::

        with dynesty.pool.Pool(16, loglike, prior_transform) as pool:
            dns = DynamicNestedSampler(pool.loglike, pool.prior_transform, ndim,
                                     pool=pool)

    Also note that you have to provide the .loglike/.prior_transform attributes
    from the pool object to the Nested samper rather than your original
    functions!
    """

    def __init__(self,
                 njobs,
                 loglike,
                 prior_transform,
                 logl_args=None,
                 logl_kwargs=None,
                 ptform_args=None,
                 ptform_kwargs=None):
        self.logl_args = logl_args
        self.logl_kwargs = logl_kwargs
        self.ptform_args = ptform_args
        self.ptform_kwargs = ptform_kwargs
        self.njobs = njobs
        self.loglike_0 = loglike
        self.prior_transform_0 = prior_transform
        self.loglike = loglike_cache
        self.prior_transform = prior_transform_cache
        self.pool = None

    def __enter__(self):
        """
        Activate the pool
        """
        initargs = (self.loglike_0, self.prior_transform_0, self.logl_args
                    or (), self.logl_kwargs or {}, self.ptform_args
                    or (), self.ptform_kwargs or {})
        self.pool = mp.Pool(self.njobs, initializer, initargs)
        initializer(*initargs)
        # running this in the master process seems to help with
        # restoration of the sampler ( #403)
        return self

    def map(self, F, x):
        """ Apply the function F to the list x

        Parameters
        ==========

        F: function
        x: iterable
        """
        return self.pool.map(F, x)

    def __exit__(self, exc_type, exc_val, exc_tb):
        try:
            self.pool.terminate()
        except:  # noqa
            pass
        try:
            del (FunctionCache.loglike, FunctionCache.prior_transform,
                 FunctionCache.logl_args, FunctionCache.logl_kwargs,
                 FunctionCache.ptform_args, FunctionCache.ptform_kwargs)
        except:  # noqa
            pass

    @property
    def size(self):
        """
        Return the number of processes in the pool
        """
        return self.njobs

    def close(self):
        self.pool.close()

    def join(self):
        self.pool.join()


class JarvisFactoryAsyncPool:
    """Thread-backed map pool for Jarvis-HEP dynesty integration."""

    def __init__(self, njobs):
        self.njobs = max(1, int(njobs))
        self._executor = ThreadPoolExecutor(max_workers=self.njobs)

    @property
    def size(self):
        return self.njobs

    def submit(self, func, item):
        if self._executor is None:
            raise RuntimeError("JarvisFactoryAsyncPool has been shut down")
        return self._executor.submit(func, item)

    def wait_first_completed(self, futures):
        if self._executor is None:
            raise RuntimeError("JarvisFactoryAsyncPool has been shut down")
        return wait(set(futures), return_when=FIRST_COMPLETED)

    def map(self, func, iterable):
        if self._executor is None:
            raise RuntimeError("JarvisFactoryAsyncPool has been shut down")

        iterator = iter(iterable)
        results = []
        pending = {}
        next_idx = 0
        exhausted = False

        while pending or not exhausted:
            while not exhausted and len(pending) < self.njobs:
                try:
                    item = next(iterator)
                except StopIteration:
                    exhausted = True
                    break
                future = self.submit(func, item)
                pending[future] = next_idx
                results.append(None)
                next_idx += 1

            if not pending:
                break

            done, _ = self.wait_first_completed(pending.keys())
            for future in done:
                i = pending.pop(future)
                results[i] = future.result()

        return results

    def shutdown(self, wait_for_tasks=True, cancel_futures=True):
        if self._executor is not None:
            self._executor.shutdown(wait=wait_for_tasks, cancel_futures=cancel_futures)
            self._executor = None

    def close(self):
        self.shutdown(wait_for_tasks=True, cancel_futures=False)

    def join(self):
        return None
