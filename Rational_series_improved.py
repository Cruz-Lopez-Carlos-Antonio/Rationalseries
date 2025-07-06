import itertools
from mpmath import mp, mpf, fac, binomial, lerchphi, polygamma
from sympy import Symbol, Poly, diff, lambdify

# Precision
mp.dps = 50
t = Symbol('t')  # Defining the symbolic variable

def DP(pol_coeff, omega):#Version of the code using a symbolic approach
    p = Poly(pol_coeff[::-1], t)  # Coefficients of the Polynomial
    derivatives_list = []
    for k in range(omega + 1):
        deriv = diff(p.as_expr(), t, k)
        f = lambdify(t, deriv, modules='mpmath')  # Symbolic derivative
        derivatives_list.append(lambda x, f=f: mpf(f(x)))  
    return derivatives_list

def DPsiZ(n, z, x):
    z, x = mpf(z), mpf(x)
    if abs(z) < 1:
        return (-1)**n * fac(n) * lerchphi(z, n + 1, x)
    elif z == 1:
        return -polygamma(n, x)
    elif z == -1:
        term1 = polygamma(n, (x + 1)/2)
        term2 = polygamma(n, x / 2)
        return (1 / (2**(n + 1))) * (term1 - term2)
    else:
        raise ValueError("Not defined case")

def listPsiZDerivatives(z, omega):
    return [lambda x, n=n: DPsiZ(n, z, x) for n in range(omega + 1)]

def ComputeF(k, a_i, z, derivatives_P, derivatives_Psi_z):
    F = mpf(0)
    a_i = mpf(a_i)
    for l in range(k + 1):
        C = binomial(int(k), int(l))
        P_l = derivatives_P[l](-a_i)
        Psi_k_l = derivatives_Psi_z[k - l](a_i)
        F += C * ((-1) ** l) * P_l * Psi_k_l
    return F

def SubsetDiophantine(parts, total):
    for tup in itertools.product(range(total + 1), repeat=parts):
        if sum(tup) == total:
            yield list(tup)

def DiophantineS(num_ai, maxim_mi):
    return [list(SubsetDiophantine(num_ai, k)) for k in range(maxim_mi + 1)]

def OmegaF(i, j, m_i, a_i, D):
    c = mpf(0)
    Aux_L = m_i[:i] + m_i[i+1:]
    Aux2 = a_i[:i] + a_i[i+1:]
    for u in D[j]:
        a = mpf(1)
        for k in range(len(Aux_L)):
            b_f = binomial(int(Aux_L[k] + u[k]), int(Aux_L[k]))
            denom = (mpf(a_i[i]) - mpf(Aux2[k]))**u[k]
            a *= b_f / denom if denom != 0 else mpf('inf')
        c += a
    return c

def DenominatorProduct(i, a, m):
    a_i = mpf(a[i])
    product = mpf(1)
    for aj, mj in zip(a[:i] + a[i+1:], m[:i] + m[i+1:]):
        product *= (mpf(aj) - a_i)**(mj + 1)
    return product

def Generalization(Pol, z, ai_list, mi_list):
    max_order = max(mi_list)
    PDerivs = DP(Pol, max_order)
    PreSol = DiophantineS(len(ai_list) - 1, max_order)
    PsiDerivs = listPsiZDerivatives(z, max_order)
    suma = mpf(0)
    for i in range(len(ai_list)):
        suma1 = mpf(0)
        denom = DenominatorProduct(i, ai_list, mi_list)
        for k in range(mi_list[i] + 1):
            term = ((-1) ** k) / fac(k)
            term *= ComputeF(k, ai_list[i], z, PDerivs, PsiDerivs)
            term *= OmegaF(i, mi_list[i] - k, mi_list, ai_list, PreSol)
            suma1 += term
        suma += suma1 / denom
    return suma

Polynomial = [0]*13 + [1] #From least to Greatest"
ai=[1,3,5]
mi = [2,4,6]
zeta=1
Calculation =  Generalization(Polynomial,zeta,ai,mi)
print("Result", Calculation)
