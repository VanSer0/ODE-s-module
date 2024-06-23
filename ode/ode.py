"""Provee de métodos numéricos para resolver ODE's.

Este módulo le permite al usuario resolver numéricamente ecuaciones diferenciales ordinarias.

Examples:
    >>> import ode
    >>> import numpy as np
    >>> def f(x, t):
    ...     return x + t
    >>> ode.RK2(f, 10, 0.0, 0.0, 10.0)
    array([0.00000000e+00, 5.00000000e-01, 3.41666667e+00, 1.23750000e+01,
           3.64375000e+01, 9.82604167e+01, 2.54484375e+02, 6.46710937e+02,
           1.62894401e+03, 4.08619336e+03])



El módulo contiene las siguientes funciones:

- `Euler(f, N, xi, ti, tf)` - Retorna un arreglo con los resultados numéricos de la ODE dada por `f` utilizando el método de `Euler`.
- `RK2(f, N, xi, ti, tf)` - Retorna un arreglo con los resultados numéricos de la ODE dada por `f` utilizando el método de `Runge-Kutta 2`.
- `RK4(f, N, xi, ti, tf)` - Retorna un arreglo con los resultados numéricos de la ODE dada por `f` utilizando el método de `Runge-Kutta 4`.
"""


import numpy as np

def Euler(f, N, xi, ti, tf):

    """Calcula los resultados numéricos de una ODE usando el método de Euler

    Examples:
        >>> import numpy as np
        >>> def f(x, t):
        ...     return x + t
        >>> Euler(f, 10 , 0.0, 0.0, 10.0)
        array([  0.        ,   0.        ,   1.11111111,   4.44444444,
                12.22222222,  28.88888889,  63.33333333, 133.33333333,
               274.44444444, 557.77777778])

    Args:
        f (function): La función a ser resuelta.
        N (int): Número de pasos en la variable independiente.
        xi (float): Valor inicial de la variable dependiente.
        ti (float): Valor inicial de la variable independiente.
        tf (float): Valor final de la variable independiente.

    Returns:
        ndarray: Un arreglo de numpy con los resultados numéricos de la variable dependiente para cada vaalor de la variable independiente.
    """

    h = (tf-ti)/N
    t = np.linspace(ti, tf, N)
    x = np.zeros(len(t))
    x[0] = xi

    for i in range(len(t)-1):

        x[i+1] = x[i] + h * f(x[i], t[i])

    return x


def RK2(f, N, xi, ti, tf):

    """Calcula los resultados numéricos de una ODE usando el método de Euler

    Examples:
        >>> import numpy as np
        >>> def f(x, t):
        ...     return x + t
        >>> RK2(f, 10, 0.0, 0.0, 10.0)
        array([0.00000000e+00, 5.00000000e-01, 3.41666667e+00, 1.23750000e+01,
               3.64375000e+01, 9.82604167e+01, 2.54484375e+02, 6.46710937e+02,
               1.62894401e+03, 4.08619336e+03])

    Args:
        f (function): La función a ser resuelta.
        N (int): Número de pasos en la variable independiente.
        xi (float): Valor inicial de la variable dependiente.
        ti (float): Valor inicial de la variable independiente.
        tf (float): Valor final de la variable independiente.

    Returns:
        ndarray: Un arreglo de numpy con los resultados numéricos de la variable dependiente para cada vaalor de la variable independiente.
    """

    h = (tf-ti)/N
    t = np.linspace(ti, tf, N)
    x = np.zeros(len(t))
    x[0] = xi

    for i in range(len(t)-1):

        k1 = h*f(x[i],t[i])
        k2 = h*f(x[i] + k1/2, t[i] + h/2)

        x[i+1] = x[i] + k2

    return x


def RK4(f, N, xi, ti, tf):

    """Calcula los resultados numéricos de una ODE usando el método de Euler

    Examples:
        >>> import numpy as np
        >>> def f(x, t):
        ...     return x + t
        >>> RK4(f, 10, 0.0, 0.0, 10.0)
        array([0.00000000e+00, 7.08333333e-01, 4.52488426e+00, 1.67595245e+01,
               5.17931566e+01, 1.48574058e+02, 4.12587149e+02, 1.12952075e+03,
               3.07311407e+03, 8.33891079e+03])

    Args:
        f (function): La función a ser resuelta.
        N (int): Número de pasos en la variable independiente.
        xi (float): Valor inicial de la variable dependiente.
        ti (float): Valor inicial de la variable independiente.
        tf (float): Valor final de la variable independiente.

    Returns:
        ndarray: Un arreglo de numpy con los resultados numéricos de la variable dependiente para cada vaalor de la variable independiente.
    """

    h = (tf-ti)/N
    t = np.linspace(ti, tf, N)
    x = np.zeros(len(t))
    x[0] = xi

    for i in range(len(t)-1):

        k1 = h*f(x[i], t[i])
        k2 = h*f(x[i] + k1/2, t[i] + h/2)
        k3 = h*f(x[i] + k2/2, t[i] + h/2)
        k4 = h*f(x[i] + k3, t[i]+h)

        x[i+1] = x[i] + (1/6)*(k1+2*k2+2*k3+k4)

    return x
