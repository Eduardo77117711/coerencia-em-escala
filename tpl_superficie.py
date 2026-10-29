"""tpl_superficie.py v1.0 - TPL passo 3: erro de bit e de fase juntos, com medicao ruidosa da sindrome.

Pai: tpl_correcao.py (pin f46209a8...c30f), que provou o limiar so para erro de fase,
com sindrome perfeita, e deixou MENCIONADO: "com medicao ruidosa o limiar real e menor".

Modelo: codigo de superficie rotacionado (stim), d rodadas de sindrome, ruido de circuito:
despolarizante p apos cada porta de dois qubits (bits e fases), despolarizante p nos dados
a cada rodada, erro p em cada medicao e em cada reset. Decodificador: emparelhamento de
peso minimo (pymatching) construido do proprio modelo de erro do circuito.
Memoria Z (protege contra erro de bit logico) e memoria X (contra erro de fase logico).

PROVA    - p=0: zero erros logicos (deterministico).
           Abaixo do limiar: p_L cai com d (3 -> 5 -> 7) com separacao estatistica.
           Acima do limiar: aumentar d nao ajuda.
           Mutantes do decodificador (nenhuma correcao, correcao aleatoria) detectados.
MOSTRA   - limiar estimado pelo cruzamento das curvas d=3 e d=7 (interpolacao log-log).
MENCIONA - extrapolacao para trilhoes de qubits logicos a partir do fator de supressao
           medido (Lambda = p_L(d) / p_L(d+2)); limites do modelo.
Sem booleanos de decisao: o resultado e um estado textual. Semente fixa: 1997.
"""
import math
import sys

import numpy as np
import pymatching
import stim

SEMENTE = 1997
DS = (3, 5, 7)
PS = (0.002, 0.004, 0.006, 0.008, 0.010, 0.012, 0.015)
SHOTS = 40_000
SHOTS_BAIXO = 400_000  # mais tiros onde os erros sao raros, para a separacao estatistica ter poder
PS_BAIXO = (0.002, 0.004)
ALVOS = (2.5e12, 5e12, 2e13)


def circuito(base, d, p):
    return stim.Circuit.generated(
        "surface_code:rotated_memory_" + base, distance=d, rounds=d,
        after_clifford_depolarization=p, before_round_data_depolarization=p,
        before_measure_flip_probability=p, after_reset_flip_probability=p)


def medir(base, d, p, shots=SHOTS, decod="mwpm", semente=SEMENTE):
    c = circuito(base, d, p)
    amostrador = c.compile_detector_sampler(seed=semente)
    det, obs = amostrador.sample(shots, separate_observables=True)
    if decod == "mwpm":
        m = pymatching.Matching.from_detector_error_model(c.detector_error_model(decompose_errors=True))
        pred = m.decode_batch(det)
    elif decod == "nenhum":
        pred = np.zeros_like(obs)
    else:  # aleatorio
        pred = np.random.default_rng(semente).integers(0, 2, size=obs.shape, dtype=obs.dtype)
    erros = int(np.any(pred != obs, axis=1).sum())
    return erros, shots


def intervalo(k, n, z=3.0):
    """Intervalo de Wilson a 3 sigma para k erros em n tentativas."""
    ph = k / n
    den = 1 + z * z / n
    centro = (ph + z * z / (2 * n)) / den
    meia = z * math.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / den
    return max(0.0, centro - meia), min(1.0, centro + meia)


def taxa(tabela, base, d, p):
    k, n = tabela[(base, d, p)]
    return k / n


def cruzamento(tabela, base):
    """p onde p_L(d=3) e p_L(d=7) se cruzam, interpolado em log-log."""
    for p0, p1 in zip(PS, PS[1:]):
        a0 = taxa(tabela, base, 7, p0) - taxa(tabela, base, 3, p0)
        a1 = taxa(tabela, base, 7, p1) - taxa(tabela, base, 3, p1)
        if a0 < 0 <= a1:
            def lg(dd, pp):
                k, n = tabela[(base, dd, pp)]
                return math.log(max(k, 0.5) / n)
            f0 = lg(7, p0) - lg(3, p0)
            f1 = lg(7, p1) - lg(3, p1)
            t = f0 / (f0 - f1)
            return math.exp(math.log(p0) + t * (math.log(p1) - math.log(p0)))
    return None


def main():
    desvios = []
    # PROVA 1: sem ruido, nenhum erro logico
    for base in ("x", "z"):
        k, n = medir(base, 5, 0.0, shots=2000)
        if k != 0:
            desvios.append("p=0 gerou %d erros na memoria %s" % (k, base))
    # Tabela completa
    tabela = {}
    for base in ("x", "z"):
        for d in DS:
            for p in PS:
                tabela[(base, d, p)] = medir(base, d, p, shots=SHOTS_BAIXO if p in PS_BAIXO else SHOTS)
    # PROVA 2 e 3: supressao abaixo do limiar, sem ganho acima
    for base in ("x", "z"):
        for p in PS_BAIXO:
            iv = [intervalo(*tabela[(base, d, p)]) for d in DS]
            if not (iv[1][1] < iv[0][0] and iv[2][1] < iv[1][0]):
                desvios.append("supressao nao separada estatisticamente em %s p=%.3f" % (base, p))
        k3, k7 = tabela[(base, 3, 0.015)][0], tabela[(base, 7, 0.015)][0]
        if k7 < k3:
            desvios.append("acima do limiar d=7 ainda melhorou em %s" % base)
    # PROVA 4: mutantes do decodificador, a p=0.004 e d=5
    mut = {}
    k_ok, n_ok = tabela[("z", 5, 0.004)]
    for nome in ("nenhum", "aleatorio"):
        k, n = medir("z", 5, 0.004, decod=nome)
        mut[nome] = (k / n, intervalo(k, n)[0] > intervalo(k_ok, n_ok)[1])
    if desvios:
        raise SystemExit("DIVERGENTE: " + "; ".join(desvios))

    print("PROVA: codigo de superficie rotacionado, ruido de circuito (bit+fase), medicao e reset ruidosos, %d rodadas = d" % DS[-1])
    print("PROVA: p=0 -> 0 erros logicos nas memorias X e Z (2000 tiros cada)")
    print("PROVA: p_L por distancia, %d tiros por ponto (%d em p=0.002 e 0.004), semente %d (intervalo a 3 sigma)" % (SHOTS, SHOTS_BAIXO, SEMENTE))
    for base in ("z", "x"):
        print("  memoria %s (%s logico)" % (base.upper(), "bit" if base == "z" else "fase"))
        for p in PS:
            celulas = []
            for d in DS:
                k, n = tabela[(base, d, p)]
                lo, hi = intervalo(k, n)
                celulas.append("d=%d: %.2e [%.1e,%.1e]" % (d, k / n, lo, hi))
            print("    p=%.3f  " % p + "  ".join(celulas))
    print("PROVA: abaixo do limiar (p=0.002 e 0.004) cada aumento de d derruba p_L, intervalos separados")
    print("PROVA: acima do limiar (p=0.015) aumentar d nao ajuda")
    for nome, (r, ok) in mut.items():
        print("PROVA: mutante do decodificador '%s' -> %s (p_L=%.2e contra %.2e do decodificador real)"
              % (nome, "MORTO" if ok else "SOBREVIVEU", r, k_ok / n_ok))
    print()
    lims = {b: cruzamento(tabela, b) for b in ("z", "x")}
    for b, v in lims.items():
        print("MOSTRA: limiar estimado (memoria %s, cruzamento d=3 x d=7): %s" % (b.upper(), "%.4f" % v if v else "fora da faixa medida"))
    print()
    base, p = "z", 0.002
    r3, r5, r7 = (taxa(tabela, base, d, p) for d in DS)
    lam = math.sqrt((r3 / r5) * (r5 / r7))
    pl7 = r7
    print("MENCIONA (extrapolacao, nao simulada): em p=%.3f o fator de supressao medido e Lambda=%.2f por +2 de distancia" % (p, lam))
    for n_log in ALVOS:
        d, pl = 7, pl7
        while n_log * pl > math.log(2) and d < 201:
            d += 2
            pl /= lam
        print("  N=%.1e qubits logicos, 1 rodada a 50%%: d~%d, ~2d^2=%d fisicos por logico, total ~%.1e" % (n_log, d, 2 * d * d, n_log * 2 * d * d))
    print("  limites: ruido de circuito padrao do stim, erros independentes, sem vazamento nem erros correlacionados; uma rodada logica.")
    if not all(ok for _, ok in mut.values()):
        raise SystemExit("ESTADO: BLOQUEADO (mutante sobreviveu)")
    print()
    print("ESTADO: PROVADO o limiar com bit+fase e sindrome ruidosa; MOSTRADO o valor do limiar; MENCIONADA a escala")


if __name__ == "__main__":
    main() if sys.argv[1:] == [] else SystemExit("uso: python3 tpl_superficie.py")
