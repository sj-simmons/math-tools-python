
from sjLib import *

n=10000000
count = 0;

for i in range(2,n):
    if isprime(i):
        count += 1
print("number of primes less than",n,"is",count)
