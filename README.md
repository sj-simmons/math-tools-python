# math-tools-python

Assorted math-related utilities in Python.

---

#### Notes
  * The file [polynomial.py](polynomial.py) is a Python Polynomial implementation. To see the documentation:
    ```
    >>> import polynomial
    >>> help(polynomial)
    ```

  * [bernoulli.py](bernoulli.py) computes Bernoulli numbers via their generating series, so via polynomial
    multiplication.   This is not the most efficient way to compute those numbers.  In a Data Structures and Algorithms
    course taught in Shanghai, we mainly used this code to compare the speed of our c++ vs Python implementations of
    polynomial mulitplication.

    The Python version can be run from the command line with `python3 bernoulli.py 100`, which computes B_100.

