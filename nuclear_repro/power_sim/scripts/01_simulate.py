#!/usr/bin/env python3
"""Paso 1: simula genomas autosómicos completos (longitudes hg19) bajo la demografía del README.

uso: python3 scripts/01_simulate.py <T_ND_ka> <semilla> <dir_salida>
salida: <dir_salida>/chrN.npz con 'pos' (posiciones) y 'hap' (sitios × 10 haplotipos, uint8 0/1)
        y <dir_salida>/meta.json

Orden de haplotipos: mbuti(0,1) den3(2,3) vindija(4,5) chag(6,7) sima(8,9)

El chimpancé no se simula: con el modelo de sitios infinitos, el alelo 0 es siempre el ancestral,
así que el grupo externo es un "chimpancé ideal" que polariza sin errores. Esto evita simular la
rama de 6,5 Ma, que multiplicaba el tiempo de cálculo por ~20.
"""
import json
import sys
from pathlib import Path

import msprime
import numpy as np

GEN = 29.0
MU = 1.25e-8
REC = 1e-8
# Longitudes de los autosomas en hg19 (GRCh37)
HG19 = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,
        146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540,
        102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566]


def demography(t_nd_ka):
    g = lambda ka: ka * 1000.0 / GEN  # noqa: E731  (ka -> generaciones)
    d = msprime.Demography()
    d.add_population(name="AFR", initial_size=15_000)
    d.add_population(name="DEN", initial_size=3_000)
    d.add_population(name="NEA", initial_size=3_000)
    d.add_population(name="ARC", initial_size=5_000)
    d.add_population(name="HOM", initial_size=10_000)
    d.add_population_split(time=g(t_nd_ka), derived=["DEN", "NEA"], ancestral="ARC")
    d.add_population_split(time=g(600), derived=["AFR", "ARC"], ancestral="HOM")
    return d, g


def main(t_nd_ka, seed, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    dem, g = demography(t_nd_ka)
    samples = [
        msprime.SampleSet(1, population="AFR", time=0, ploidy=2),
        msprime.SampleSet(1, population="DEN", time=g(65), ploidy=2),
        msprime.SampleSet(1, population="NEA", time=g(50), ploidy=2),
        msprime.SampleSet(1, population="NEA", time=g(80), ploidy=2),
        msprime.SampleSet(1, population="NEA", time=g(430), ploidy=2),
    ]
    rng = np.random.default_rng(seed)
    nsites = 0
    for i, L in enumerate(HG19, start=1):
        s1, s2 = (int(x) for x in rng.integers(1, 2**31 - 1, size=2))
        ts = msprime.sim_ancestry(samples=samples, demography=dem, sequence_length=L,
                                  recombination_rate=REC, random_seed=s1, discrete_genome=False)
        ts = msprime.sim_mutations(ts, rate=MU, random_seed=s2, discrete_genome=False,
                                   model=msprime.BinaryMutationModel())
        hap = ts.genotype_matrix().astype(np.uint8)          # sitios × 12
        pos = ts.tables.sites.position.astype(np.float64)
        np.savez_compressed(out / f"chr{i}.npz", pos=pos, hap=hap)
        nsites += len(pos)
        print(f"chr{i}: {len(pos):,} sitios", flush=True)
    json.dump({"T_ND_ka": t_nd_ka, "seed": seed, "genome_length": sum(HG19),
               "segregating_sites": nsites, "msprime": msprime.__version__,
               "mu": MU, "rec": REC, "gen_time": GEN},
              open(out / "meta.json", "w"), indent=1)
    print(f"Total: {nsites:,} sitios segregantes en {sum(HG19):,} pb")


if __name__ == "__main__":
    main(float(sys.argv[1]), int(sys.argv[2]), sys.argv[3])
