import numpy as np
from scipy.stats import norm
# Coherencia interna de los porcentajes de Meyer: los dos especimenes deben dar el mismo cociente
# si miden la misma cantidad subyacente. Cobertura por rama: mediana simulada, escenario 470 ka.
sp={"AT-5431":(0.43,0.09,397,384),"AT-5482":(0.39,0.07,677,652)}
out={}
for k,(fN,fD,cN,cD) in sp.items():
    kN,kD=fN*cN,fD*cD
    # EE del log-cociente por delta: var(log f) = (1-f)/(f*n)
    v=(1-fN)/(fN*cN)+(1-fD)/(fD*cD)
    lr=np.log(fN/fD); se=np.sqrt(v)
    out[k]=(fN/fD,lr,se)
    print(f"{k}: f_N/f_D = {fN/fD:.2f}  IC95 = [{np.exp(lr-1.96*se):.2f}, {np.exp(lr+1.96*se):.2f}]  (n_N={cN}, n_D={cD})")
(a,la,sa),(b,lb,sb)=out["AT-5431"],out["AT-5482"]
z=(la-lb)/np.sqrt(sa**2+sb**2); p=2*(1-norm.cdf(abs(z)))
print(f"\nDiferencia entre los dos cocientes: Z = {z:.2f}, p = {p:.2f}")
print("-> coherentes entre si" if p>0.05 else "-> discrepantes")
# cociente combinado
lc=(la/sa**2+lb/sb**2)/(1/sa**2+1/sb**2); sc=np.sqrt(1/(1/sa**2+1/sb**2))
print(f"Cociente combinado = {np.exp(lc):.2f}  IC95 = [{np.exp(lc-1.96*sc):.2f}, {np.exp(lc+1.96*sc):.2f}]")
print(f"\nSimulacion: 1,15 (440 ka) / 1,76 (470 ka) / 5,03 (550 ka)")
for esc,v in (("440",1.15),("470",1.76),("550",5.03)):
    zz=(lc-np.log(v))/sc
    print(f"  escenario {esc} ka: Z = {zz:.1f} -> {'compatible' if abs(zz)<3 else 'INCOMPATIBLE'}")
