#!/usr/bin/python
# -*- coding: utf-8 -*-
import time
import random
import pytest
from kaic.tools.lru import lru_cache
from __builtin__ import classmethod


class TestDecorators:
    def test_lru(self):
        class LRUTest(object):
            """class"""

            def __init__(self):
                self.num = 0

            @lru_cache(maxsize=10, timeout=3)
            def test_method(self, num):
                """test_method_doc"""
                self.num += num
                return self.num

        @lru_cache(maxsize=10, timeout=3)
        def test_func(num):
            """test_func_doc"""
            return num

        @lru_cache(maxsize=10, timeout=3)
        def test_func_time(num):
            """test_func_time_doc"""
            return time.time()

        @lru_cache(maxsize=10, timeout=None)
        def test_func_args(*args, **kwargs):
            return random.randint(1, 10000000)

        # Init vars:
        c1 = LRUTest()
        c2 = LRUTest()
        m1 = c1.test_method
        m2 = c2.test_method
        f1 = test_func

        # Test basic caching functionality:
        assert m1(1) == m1(1)
        assert c1.num == 1  # c1.num now equals 1 - once cached, once real
        assert f1(1) == f1(1)

        # Test caching is different between instances - once cached, once not cached:
        assert m1(2) != m2(2)
        assert m1(2) != m2(2)

        # Validate the cache_clear funcionality only on one instance:
        prev1 = m1(1)
        prev2 = m2(1)
        prev3 = f1(1)
        m1.cache_clear()
        assert m1(1) != prev1
        assert m2(1) == prev2
        assert f1(1) == prev3

        # Validate the docstring and the name are set correctly:
        assert m1.__doc__ == "test_method_doc"
        assert f1.__doc__ == "test_func_doc"
        assert m1.__name__ == "test_method"
        assert f1.__name__ == "test_func"

        # Test the limit of the cache, cache size is 10, fill 15 vars,
        # the first 5 will be overwritten for each and the other 5 are untouched. Test that:
        c1.num = 0
        c2.num = 10
        m1.cache_clear()
        m2.cache_clear()
        f1.cache_clear()
        temp_list = map(lambda i: (test_func_time(i), m1(i), m2(i)), range(15))

        for i in range(5, 10):
            assert temp_list[i] == (test_func_time(i), m1(i), m2(i))
        for i in range(0, 5):
            assert temp_list[i] != (test_func_time(i), m1(i), m2(i))
        # With the last run the next 5 vars were overwritten, now it should have only 0..4 and 10..14:
        for i in range(5, 10):
            assert temp_list[i] != (test_func_time(i), m1(i), m2(i))

        # Test different vars don't collide:
        assert test_func_args(1) != test_func_args('1')
        assert test_func_args(1.0) != test_func_args('1.0')
        assert test_func_args(1.0) != test_func_args(1)
        assert test_func_args(None) != test_func_args('None')
        assert test_func_args(test_func) == test_func_args(test_func)
        assert test_func_args(LRUTest) == test_func_args(LRUTest)
        assert test_func_args(object) == test_func_args(object)
        assert test_func_args(1, num=1) != test_func_args(1, num='1')
        # Test the sorting of kwargs:
        assert test_func_args(1, aaa=1, bbb=2) == test_func_args(1, bbb=2, aaa=1)
        assert test_func_args(1, aaa='1', bbb=2) != test_func_args(1, bbb=2, aaa=1)

        # Sanity validation of values
        c1.num = 0
        c2.num = 10
        m1.cache_clear()
        m2.cache_clear()
        f1.cache_clear()
        assert (f1(0), m1(0) == m2(0)), (0, 0, 10)
        assert (f1(0), m1(0) == m2(0)), (0, 0, 10)
        assert (f1(1), m1(1) == m2(1)), (1, 1, 11)
        assert (f1(2), m1(2) == m2(2)), (2, 3, 13)
        assert (f1(2), m1(2) == m2(2)), (2, 3, 13)
        assert (f1(3), m1(3) == m2(3)), (3, 6, 16)
        assert (f1(3), m1(3) == m2(3)), (3, 6, 16)
        assert (f1(4), m1(4) == m2(4)), (4, 10, 20)
        assert (f1(4), m1(4) == m2(4)), (4, 10, 20)

        # Test timeout - sleep, it should refresh cache, and then check it was cleared:
        prev_time = test_func_time(0)
        assert test_func_time(0) == prev_time
        assert m1(4) == 10
        assert m2(4) == 20
        time.sleep(3.5)
        assert test_func_time(0) != prev_time
        assert m1(4) != 10
        assert m2(4) != 20
