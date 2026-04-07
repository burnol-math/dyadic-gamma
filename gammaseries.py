"""
Computation of gamma for a given target decimal precision and level

Author: Jean-François Burnol
Created: April 2, 2026 v0.1.0
This version: April 7, 2026 v0.2.0

- add silent and trunc booolean optional parameters,
- return a string in place of printing it.

TODO: should I return an object of type mpf rather?

© Jean-François Burnol, 2026

This script is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License.
Full license text: https://creativecommons.org/licenses/by-sa/4.0/

ANY RE-USE OR PLAGIARIZING BY AN ARTIFICIAL INTELLIGENCE WITHOUT
PROPER ATTRIBUTION IS STRICTLY FORBIDDEN AND WILL GET PUNISHED

"""

__version__ = "0.2.0"
__date__    = "2026-04-07"
__author__  = "Jean-François Burnol"

from math import log, ceil
from mpmath import mp, mpf, nstr

def gamma_ell(N:int,
              ell:int = 8,
              silent:bool = False,
              trunc:bool = False) -> str:
    """Computes gamma with N decimal digits, using level ell

    :param N int: requested number of decimal places
    :param ell int: (optional, defaults to 8) the level ell
    :param silent bool: (optional, defaults to False) whether to print
                        some extra info
    :param trunc bool: (optional, defaults to False) whether to round or
                       truncate the return value
    :rtype: str
    :return: String with N significant figures of the decimal expansion
             of Euler's constant, either rounded or truncated as per trunc
             optional parameter.

    """
    P_target = ceil((N + 1) * log(10, 2))
    P_eps = P_target + min(ell, 10)
    P_max = P_eps + 12

    mp.prec = P_max
    current_prec = P_max
    epsilon = mpf(2) ** -P_eps

    ell_m_one = ell - 1
    two_to_the_ell_m_one = 2 ** ell_m_one
    two_to_the_ell = 2 * two_to_the_ell_m_one

    m = 0
    two_to_the_m_p_one = 2

    S = sum(mpf(1)/n for n in range(two_to_the_ell_m_one - 1, 0, -1))
    S -= (ell-1) * mp.log(mpf(2))

    invpowers = [ 1/mpf(n) for n in range(two_to_the_ell - 1,
                                          two_to_the_ell_m_one - 1,
                                          -1)
                 ]

    cms = [ 0 ]
    ems = [ 0 ]


    while True:
        m += 1
        two_to_the_m_p_one *= 2

        mp.prec = current_prec

        em = mpf(0)
        binomj = 1  # binomial coefficients of type int
        for j in range(1, m):  # e0 = 0 de toute façon
            binomj = binomj * (m - j + 2) // j
            em += mpf(binomj) * ems[-j]

        em += two_to_the_m_p_one
        em /= two_to_the_m_p_one - 2
        ems.append(em)

        # keep inverse powers from smallest to largest
        i = 0
        for n in range(two_to_the_ell - 1,
                       two_to_the_ell_m_one - 1,
                       -1):
            invpowers[i] /= n
            i += 1

        # would it be more efficient to use invpowers[::-1]
        # and keep in ordering induced from natural numbers?
        cm = em / ( m + 1 ) * sum(x for x in invpowers)

        if cm < epsilon:
            break

        cms.append(cm)

        mp.prec = P_max

        S += (-1)**(m-1) * cm

        # perhaps 30 binary digits is a bit overkill.
        current_prec = max(30, current_prec - ell_m_one)
        mp.prec = current_prec

    mp.prec = P_max

    if not silent:
        print(f"ell is {ell}")
        print(f"Last used: m = {m-1}, |cm| = {nstr(cms[-1], 4, strip_zeros=False)}")
        print(f" not used: m = {m}, |cm| = {nstr(cm, 4, strip_zeros=False)}")

        # ATTENTION ! nstr supprime trailing zeros par défaut !!
        S_string = nstr(S, N + 3, strip_zeros=False)
        # print last digits + 3 more
        print(f"Last digits are {S_string[-12:-3]}({S_string[-3:]})")

        print(f"gamma {'truncated' if trunc else 'rounded'} "
              f"to {N} significant figures is (as a string):")

    if trunc:
        # we hope this gives correct truncation...
        return nstr(S, N + 6, strip_zeros=False)[:-6]
    else:
        return nstr(S, N, strip_zeros=False)


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 2:
        print(gamma_ell(int(sys.argv[1]), int(sys.argv[2])))
    elif len(sys.argv) > 1:
        print(gamma_ell(int(sys.argv[1])))
    else:
        print("\ngamma_ell(N[,ell]) gives gamma to N decimal digits.  Default: ell=8")
