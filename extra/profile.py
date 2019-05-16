import time
from functools import wraps


def timing(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        elapsed_time = time.time() - start_time

        if 'log_time' in kwargs:
            name = kwargs.get('log_name', func.__name__.upper())
            kwargs['log_time'][name] = elapsed_time
        else:
            print('Elapsed time of function {}: {} seconds'.format(func.__name__, elapsed_time))

        return result

    return wrapper
