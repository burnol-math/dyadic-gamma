"""
Computation of gamma with the \\ell = 4 series, for a given target decimal precision

Author: Jean-François Burnol
Created: April 2, 2026

This script is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License.
Full license text: https://creativecommons.org/licenses/by-sa/4.0/
"""

__version__ = "0.1.0"
__date__    = "2026-04-02"
__author__  = "Jean-François Burnol"

__all__ = ["foo", "__version__", "__date__", "__author__"]

from mpmath import mp, mpf, log, nstr

def gamma_ell4(N: int) -> None:
    """Computes gamma with N decimal digits
    """
    log2_10 = mp.log(10) / mp.log(2)
    P = int(N * log2_10) + 1
    Pthreshold = P + 10
    Pmax = P + 20

    mp.prec = Pmax

    S = sum(mpf(1)/n for n in range(1, 8)) - 3 * log(mpf(2))

    cms = [ 0 ]
    ems = [ 0 ]

    current_prec = Pmax
    threshold = mpf(2) ** (-Pthreshold)

    invpowers = [ 1/mpf(n) for n in range(8,16) ]
    m = 0
    twotothemplusone = 2

    while True:
        m += 1
        twotothemplusone *= 2

        mp.prec = current_prec

        i = 0
        for n in range(8, 16):
            invpowers[i] /= n
            i += 1

        em = 0
        bj = 1
        for j in range(1, m ):  # e0 = 0 de toute façon
            bj = bj * (m - j + 2) // j
            em += mpf(bj) * ems[-j]

        em += twotothemplusone
        em /= twotothemplusone - 2
        ems.append(em)

        cm = em / ( m + 1 ) * sum(x for x in invpowers)

        if cm < threshold:
            break

        cms.append(cm)

        mp.prec = Pmax

        S += (-1)**(m-1) * cm

        current_prec = max(20, current_prec - 3)
        mp.prec = current_prec

    mp.prec = Pmax

    # Print results to stdout
    print(f"Last used: m = {m-1}, |cm| = {nstr(cms[-1], 4)}")
    print(f" not used: m = {m}, |cm| = {nstr(cm, 4)}")

    Sstring = nstr(S, N + 3)
    # print(Sstring)
    print(f"gamma = {Sstring[:N+2]}({Sstring[-3:]})")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        gamma_ell4(int(sys.argv[1]))
    else:
        print("\ngamma_ell4(N) gives gamma to N decimal digits")
