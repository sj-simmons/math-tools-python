#!/usr/bin/env python3

'''Provides basic number theory functions.

Functions for computing gcds, factoring, etc.
'''

import time, math

__author__ = "Scott Simmons (2014)"
__status__ = "Production"

def gcd(n, m):
    return m if n % m == 0 else gcd(m, n % m)

def xgcd(a, b):
    ''' a,b positive integers. Returns list [gcd(a,b),m,n] satisfying:
        gcd(a,b) = m * a + n * b '''
    s = [1, 0]; t = [0, 1]
    while True:
        quot = -(a // b)
        a = a % b
        s[0] += quot * s[1]; t[0] += quot * t[1]
        if a == 0:
            return [b, s[1], t[1]]
        quot = -(b // a)
        b = b % a
        s[1] += quot * s[0]; t[1] += quot * t[0]
        if b == 0:
            return [a, s[0], t[0]]

def sieve(n = 1000000):
    ''' Sieve of Eratosthenes. Returns list of primes <= n.  '''
    primes = []
    sieve = [True] * (n + 1)
    for p in range (2, n+1):
        if sieve[p]:
            primes.append(p)
            for i in range(p * p, n + 1, p):
                sieve[i] = False
    return primes

def isprime(n):
    ''' isprime(n) - Test whether n is prime using a variety of pseudoprime tests.'''
    if n in [2,3,5,7,11,13,17,19,23,29]: return True
    return isprimeE(n, 2) and isprimeE(n, 3) and isprimeE(n, 5)

def isprimeF(n, b):
    ''' isprimeF(n) - Test whether n is prime or a Fermat pseudoprime to base b.'''
    return pow(b, n-1, n) == 1

def isprimeE(n, b):
    ''' isprimeE(n) - Test whether n is prime or an Euler pseudoprime to base b.'''
    if not isprimeF(n, b): return False
    r = n-1
    while r % 2 == 0: r //= 2
    c = pow(b, r, n)
    if c == 1: return True
    while True:
        if c == 1: return False
        if c == n-1: return True
        c = pow(c, 2, n)

def factor(n):
    ''' factor(n) - Find a prime factor of n using a variety of methods.'''
    if isprime(n): return n
    for fact in [2,3,5,7,11,13,17,19,23,29]:
        if n % fact == 0: return fact
    return factorPR(n)  # Needs work - no guarantee that a prime factor will be returned

def factors(n):
    ''' factors(n) - Return a sorted list of the prime factors of n.'''
    if isprime(n):
        return [n]
    fact = factor(n)
    if fact == 1: return "Unable to factor "+str(n)
    facts = factors(n // fact) + factors(fact)
    facts.sort()
    return facts

def factorPR(n):
    ''' factorPR(n) - Find a factor of n using the Pollard Rho method.
    Note: This method will occasionally fail.'''
    for slow in [2,3,4,6]:
        numsteps = 2 * math.floor(math.sqrt(math.sqrt(n))); fast=slow; i=1
        while i < numsteps:
            slow = (slow*slow + 1) % n
            i = i + 1
            fast = (fast*fast + 1) % n
            fast = (fast*fast + 1) % n
            g = gcd(fast-slow,n)
            if (g != 1):
                if (g == n):
                    break
                else:
                    return g
    return 1

def eulerphi(n):
    ''' eulerphi(n) - Computer Euler's Phi function of n - the number of integers
    strictly less than n which are coprime to n.  Otherwise defined as the order
    of the group of integers mod n.'''
    thefactors = factors(n)
    thefactors.sort()
    phi = 1
    oldfact = 1
    for fact in thefactors:
        if fact==oldfact:
            phi = phi * fact
        else:
            phi = phi * (fact-1)
            oldfact = fact
    return phi

class Timer:
    def __init__(self, startTime = time.time(), endTime = time.time()):
        self.__startTime = startTime
        self.__endTime = startTime       #            Timer
                                         #  =============================
    def getStartTime(self):              #      -startTime : float
        return self.__startTime          #        -endTime : float
                                         #  -----------------------------
    def getEndTime(self):                #     Timer(startTime : float,
        return self.__endTime            #             endTime : float)
                                         #     getStartTime() : float
    def start(self):                     #       getEndTime() : float
        self.__startTime = time.time()   #            start() : None
                                         #             stop() : float
    def stop(self):                      #   getElapsedTime() : float
        self.__endTime = time.time()

    def getElapsedTime(self):
        return "Elapsed time: "+str(format(self.__endTime - self.__startTime,".4f"))+\
                " seconds. ("+str(format((self.__endTime - self.__startTime)/60,".2f"))+" minutes.)"

