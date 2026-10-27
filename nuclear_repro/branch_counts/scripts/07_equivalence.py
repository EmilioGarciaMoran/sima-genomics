import numpy as np
from scipy.stats import norm
# Sitios informativos proyectados (branch_counts/README.md, apartado 3 bis)
N = {"Femur XIII":31, "AT-5431":161, "AT-5482":274, "AT-5431+AT-5482":436}
A = 0.05   # alfa de la prueba de equivalencia (TOST): IC del 90%
Z = norm.ppf(1-A)   # 1.645
print("Con D=0 el error tipico es raiz(1/n). Los sitios estan muy separados en el genoma,")
print("asi que apenas hay ligamiento: el binomial vale aqui (no el factor 2,7 de datos densos).\n")
print(f"{'especimen':20s} {'n':>6s} {'EE':>7s} {'delta minimo':>13s}  {'P(equivalencia | D=0)':>22s}")
for k,n in N.items():
    se = np.sqrt(1/n)
    dmin = Z*se                      # delta minimo alcanzable (solo si D observada = 0)
    # probabilidad de establecer equivalencia con delta=0.10 si la verdad es D=0
    for d in (0.10,):
        lim = d - Z*se               # |D_obs| debe ser menor que esto
        p = 0 if lim<=0 else (2*norm.cdf(lim/se)-1)
    print(f"{k:20s} {n:6d} {se:7.4f} {dmin:13.3f}  {p:21.0%}")

print("\n--- Cuantos sitios harian falta para un delta dado (con D observada = 0) ---")
for d in (0.15,0.10,0.08,0.05,0.03):
    print(f"  delta = {d:.2f}  ->  n > {int(np.ceil((Z/d)**2)):,} sitios informativos")

print("\n--- Si la verdad es la prediccion de Meyer (D = 0.68) ---")
n=436; se=np.sqrt((1-0.68**2)/n)
print(f"  AT-5431+AT-5482: D=0.68, EE={se:.3f}, Z={0.68/se:.1f}")
print(f"  probabilidad de rechazar D=0 (Z>3): {1-norm.cdf(3*np.sqrt(1/n), loc=0.68, scale=se):.1%}")
