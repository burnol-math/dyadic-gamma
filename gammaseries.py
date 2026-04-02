"""
Computation of gamma for a given target decimal precision and level

Author: Jean-François Burnol
Created: April 2, 2026

This script is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License.
Full license text: https://creativecommons.org/licenses/by-sa/4.0/
"""

__version__ = "0.1.0"
__date__    = "2026-04-02"
__author__  = "Jean-François Burnol"

__all__ = ["gamma_ell", "__version__", "__date__", "__author__"]

from math import log
from mpmath import mp, mpf, nstr

def gamma_ell(N: int, ell: int = 8) -> None:
    """Computes gamma with N decimal digits, using level ell
    """
    P = int(N * log(10) / log(2)) + 1
    P_threshold = P + 10
    P_max = P + 20

    mp.prec = P_max

    ell_m_one = ell - 1
    two_to_the_ell_m_one = 2 ** ell_m_one
    two_to_the_ell = 2 * two_to_the_ell_m_one

    m = 0
    two_to_the_m_p_one = 2

    S = sum(mpf(1)/n for n in range(1, two_to_the_ell_m_one))
    S -= (ell-1) * mp.log(mpf(2))

    invpowers = [ 1/mpf(n) for n in range(two_to_the_ell_m_one, two_to_the_ell) ]

    cms = [ 0 ]
    ems = [ 0 ]

    current_prec = P_max
    threshold = mpf(2) ** (-P_threshold)

    while True:
        m += 1
        two_to_the_m_p_one *= 2

        mp.prec = current_prec

        i = 0
        for n in range(two_to_the_ell_m_one, two_to_the_ell):
            invpowers[i] /= n
            i += 1

        em = 0
        bj = 1
        for j in range(1, m ):  # e0 = 0 de toute façon
            bj = bj * (m - j + 2) // j
            em += mpf(bj) * ems[-j]

        em += two_to_the_m_p_one
        em /= two_to_the_m_p_one - 2
        ems.append(em)

        cm = em / ( m + 1 ) * sum(x for x in invpowers)

        if cm < threshold:
            break

        cms.append(cm)

        mp.prec = P_max

        S += (-1)**(m-1) * cm

        current_prec = max(30, current_prec - ell_m_one)
        mp.prec = current_prec

    mp.prec = P_max

    # Print results to stdout
    print(f"ell is {ell}")
    print(f"Last used: m = {m-1}, |cm| = {nstr(cms[-1], 4, strip_zeros=False)}")
    print(f" not used: m = {m}, |cm| = {nstr(cm, 4, strip_zeros=False)}")

    # ATTENTION ! nstr supprime trailing zeros par défaut !!
    S_string = nstr(S, N + 3, strip_zeros=False)
    print(f"gamma = {S_string[:N+2]}({S_string[-3:]})")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 2:
        gamma_ell(int(sys.argv[1]), int(sys.argv[2]))
    elif len(sys.argv) > 1:
        gamma_ell(int(sys.argv[1]))
    else:
        print("\ngamma_ell(N[,ell]) gives gamma to N decimal digits.  Default: ell=8")
