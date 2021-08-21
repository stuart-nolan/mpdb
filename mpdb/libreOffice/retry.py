#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mpdb: Material Property Data Base as python module
retry.py: retry decorator based on ref below.

Ref: https://wiki.python.org/moin/PythonDecoratorLibrary#Retry
Revision Date: 2018/04/28
"""
from time import sleep
from functools import wraps

def retry(tries=2, delay=3, backoff=2, logger=None, exceptions=(Exception),
          msgf=print):
    """Retry calling the decorated function using an exponential backoff.

    http://www.saltycrane.com/blog/2009/11/trying-out-retry-decorator-python/
    original from: http://wiki.python.org/moin/PythonDecoratorLibrary#Retry

    :param tries: number of times to try (not retry) before giving up
    :type tries: int
    :param delay: initial delay between retries in seconds
    :type delay: int
    :param backoff: backoff multiplier e.g. value of 2 will double the delay
        each retry
    :type backoff: int
    :param logger: logger to use. If None, print
    :type logger: logging.Logger instance
    :param exceptions: tuple of exceptions to check
    :type exceptions: tuple of exceptions
    :param msgf: python function that takes 1 str arg msg (an error message) 
    :type msgf: python function that takes 1 str arg
    """
    def deco_retry(f):

        @wraps(f)
        def f_retry(*args, **kwargs):
            mtries, mdelay = tries, delay
            while mtries > 1:
                try:
                    return f(*args, **kwargs)
                except exceptions as e:
                    msg = "%s, Retrying in %d seconds..." % (str(e), mdelay)
                    if logger:
                        logger.warning(msg)
                    else:
                        msgf(msg)
                    sleep(mdelay)
                    mtries -= 1
                    mdelay *= backoff
            return f(*args, **kwargs)

        return f_retry  # true decorator

    return deco_retry

@retry(tries=4)
def test_fail(text):
    raise Exception("Fail")

if __name__ == '__main__':
    test_fail("it works!")
