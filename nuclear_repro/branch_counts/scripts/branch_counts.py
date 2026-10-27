#!/usr/bin/env python3
"""Recuento de alelos derivados compartidos por rama (método de Meyer 2016 y Massilani 2026).

Define tres estadísticos sobre un individuo de prueba X, dos genomas de referencia A y B
del mismo linaje, y un humano moderno M. El estado ancestral se toma de los primates.

  f_AB  = de los sitios donde A y B llevan el derivado y M el ancestral, fracción con X derivado
  f_A   = de los sitios donde solo A lo lleva (B y M ancestrales),       fracción con X derivado
  f_B   = de los sitios donde solo B lo lleva (A y M ancestrales),       fracción con X derivado

Massilani et al. 2026 publica, con X = Denisova 17, A = Denisova 5 y B = Vindija 33.19:
f_AB = 91 %, f_A = 45 %, f_B = 20 %.

Entrada: matrices de genotipos pseudohaploides codificados 0 = ancestral, 1 = derivado,
-1 = sin dato. Todas del mismo largo y en el mismo orden de posiciones.
"""
from __future__ import annotations

import numpy as np

__all__ = ["branch_counts", "jackknife_fraction"]


def branch_counts(x, a, b, m, blocks=None):
    """Devuelve los tres estadísticos con su recuento y, si se dan bloques, su error típico.

    x, a, b, m: arrays de enteros con 0 (ancestral), 1 (derivado) o -1 (sin dato).
    blocks:     array de enteros con el índice de bloque de cada posición, para el jackknife.

    Devuelve dict con clave 'AB', 'A' y 'B'; cada valor es un dict con n (sitios del
    conjunto de referencia que X cubre), k (de esos, cuántos lleva derivados), f = k/n
    y ee (error típico por jackknife, o None).
    """
    x, a, b, m = (np.asarray(v) for v in (x, a, b, m))
    if not (len(x) == len(a) == len(b) == len(m)):
        raise ValueError("todos los arrays deben tener el mismo largo")
    for name, v in (("x", x), ("a", a), ("b", b), ("m", m)):
        if not np.isin(v, (-1, 0, 1)).all():
            raise ValueError(f"{name} contiene valores distintos de -1, 0 y 1")

    # Los sitios de referencia exigen dato en A, B y M. X puede faltar: eso es lo que
    # distingue el denominador (sitios de referencia que X cubre) del conjunto total.
    ref = (a >= 0) & (b >= 0) & (m == 0)
    sets = {
        "AB": ref & (a == 1) & (b == 1),
        "A": ref & (a == 1) & (b == 0),
        "B": ref & (a == 0) & (b == 1),
    }
    out = {}
    for key, sel in sets.items():
        cov = sel & (x >= 0)          # sitios del conjunto que X cubre
        der = cov & (x == 1)
        n, k = int(cov.sum()), int(der.sum())
        f = k / n if n else np.nan
        ee = None
        if blocks is not None and n:
            ee = jackknife_fraction(der, cov, np.asarray(blocks))
        out[key] = {"n": n, "k": k, "f": f, "ee": ee}
    return out


def jackknife_fraction(der, cov, blocks):
    """Error típico de k/n por jackknife de bloques ponderado (Busing et al. 1999)."""
    nb = int(blocks.max()) + 1
    kb = np.bincount(blocks[der], minlength=nb).astype(float)
    nb_ = np.bincount(blocks[cov], minlength=nb).astype(float)
    K, N = kb.sum(), nb_.sum()
    keep = nb_ > 0
    if keep.sum() < 2:
        return np.nan
    kb, nb_ = kb[keep], nb_[keep]
    g = len(nb_)
    theta = K / N
    # estimador dejando fuera cada bloque
    denom = N - nb_
    with np.errstate(invalid="ignore", divide="ignore"):
        theta_j = np.where(denom > 0, (K - kb) / denom, np.nan)
    h = N / nb_                                   # peso de Busing
    tau = h * theta - (h - 1) * theta_j           # pseudovalores
    ok = np.isfinite(tau)
    if ok.sum() < 2:
        return np.nan
    tau, h = tau[ok], h[ok]
    tau_bar = tau.mean()
    var = np.sum((tau - tau_bar) ** 2 / (h - 1)) / g
    return float(np.sqrt(var)) if var > 0 else np.nan
