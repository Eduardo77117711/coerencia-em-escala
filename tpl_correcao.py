"""tpl_correcao.py v1.0 - TPL: o limiar da correcao de erros, provado e nao so mencionado.

Passo 2 depois de tpl_sincronia.py (v1.1, pin 373c9d10...96f0), que deixou o
limiar como MENCIONADO. Aqui o limiar e medido para o ruido daquele modelo:
ruido de fase independente, Z com probabilidade p = (1 - c) / 2 por qubit.

Codigo: repeticao de fase (base X), distancia d impar, decodificador por maioria.
PROVA    - erro logico p_L(d, p) por enumeracao exata dos 2^d padroes de erro
           (d = 1..15), conferido contra a soma binomial fechada.
           Monte Carlo com semente fixa conferido dentro de 5 desvios-padrao.
           Limiar observado: p < 1/2 -> p_L cai com d; p > 1/2 -> p_L sobe com d.
MOSTRA   - coerencia logica de N qubits logicos: C = (1 - 2 p_L)^N, o mesmo c^n
           de tpl_sincronia, agora com c_L = 1 - 2 p_L; distancia minima para
           2,5e12 / 5e12 / 2e13 qubits logicos a 50%.
MENCIONA - limite: este codigo corrige so erro de fase. Qubits reais tem erro
           de fase e de bit; para os dois e preciso codigo de superficie.
           Supoe sindrome medida sem erro e erros independentes.
Sem booleanos de decisao: o resultado e um estado textual.
"""
import math
import random
import sys

DS = tuple(range(1, 16, 2))
PS = (0.01, 0.05, 0.1, 0.2, 0.3, 0.45, 0.55, 0.7)
TOL = 1e-12
ALVOS = (2.5e12, 5e12, 2e13)
SEMENTE = 1997


def decodifica_maioria(padrao, d):
    """Falha logica quando mais da metade dos qubits sofreu Z."""
    return sum(padrao) > d // 2


def pl_enumeracao(d, p, decod=decodifica_maioria):
    total = 0.0
    for m in range(2 ** d):
        padrao = [(m >> i) & 1 for i in range(d)]
        k = sum(padrao)
        if decod(padrao, d):
            total += p ** k * (1 - p) ** (d - k)
    return total


def pl_fechada(d, p):
    return sum(math.comb(d, k) * p ** k * (1 - p) ** (d - k) for k in range(d // 2 + 1, d + 1))


def pl_monte_carlo(d, p, amostras, rng):
    falhas = sum(1 for _ in range(amostras) if sum(rng.random() < p for _ in range(d)) > d // 2)
    return falhas / amostras


def prova(decod=decodifica_maioria):
    """Devolve o primeiro desvio, ou None."""
    for d in DS:
        for p in PS:
            if abs(pl_enumeracao(d, p, decod) - pl_fechada(d, p)) > TOL:
                return "enumeracao != formula em d=%d p=%.2f" % (d, p)
    for p in PS:
        serie = [pl_fechada(d, p) for d in DS]
        cai = all(b < a for a, b in zip(serie, serie[1:]))
        sobe = all(b > a for a, b in zip(serie, serie[1:]))
        if (p < 0.5 and not cai) or (p > 0.5 and not sobe):
            return "limiar nao observado em p=%.2f" % p
    return None


def prova_monte_carlo():
    rng = random.Random(SEMENTE)
    n = 200_000
    piores = 0.0
    for d, p in ((3, 0.1), (5, 0.2), (7, 0.3), (9, 0.45)):
        exato = pl_fechada(d, p)
        sigma = math.sqrt(exato * (1 - exato) / n)
        z = abs(pl_monte_carlo(d, p, n, rng) - exato) / sigma
        piores = max(piores, z)
        if z > 5:
            return None, "Monte Carlo fora de 5 sigma em d=%d p=%.2f" % (d, p)
    return piores, None


MUTANTES = {
    "maioria deslocada (>= d//2)": lambda pad, d: sum(pad) >= d // 2,
    "ignora o ultimo qubit": lambda pad, d: sum(pad[:-1]) > d // 2,
    "falha com um erro so": lambda pad, d: sum(pad) >= 1,
}


def d_minimo(p, n_logicos, alvo=0.5):
    """Menor d impar com (1 - 2 p_L)^N >= alvo. Calculo em log para N gigante."""
    for d in range(1, 400, 2):
        pl = pl_fechada(d, p)
        if pl < 0.5 and n_logicos * math.log1p(-2 * pl) >= math.log(alvo):
            return d, pl
    return None, None


def main():
    desvio = prova()
    if desvio is not None:
        raise SystemExit("DIVERGENTE: " + desvio)
    z, erro_mc = prova_monte_carlo()
    if erro_mc is not None:
        raise SystemExit("DIVERGENTE: " + erro_mc)
    res = {nome: prova(m) for nome, m in MUTANTES.items()}
    print("PROVA: enumeracao exata de 2^d padroes, d=1..15, %d valores de p; igual a soma binomial (tol %g)" % (len(PS), TOL))
    print("PROVA: Monte Carlo (semente %d, 200000 amostras x 4) dentro de %.2f desvios-padrao" % (SEMENTE, z))
    print("PROVA: limiar observado em 1/2: p<1/2 p_L cai com d, p>1/2 p_L sobe com d")
    for nome, d in res.items():
        print("PROVA: mutante do decodificador '%s' -> %s" % (nome, "MORTO em " + d if d else "SOBREVIVEU"))
    print()
    print("PROVA: p_L por distancia (valores exatos)")
    for p in (0.01, 0.1, 0.3, 0.45, 0.55):
        print("  p=%.2f  " % p + "  ".join("d=%d:%.2e" % (d, pl_fechada(d, p)) for d in (1, 3, 5, 9, 15)))
    print()
    print("MOSTRA: distancia minima para N qubits logicos manterem 50% de coerencia")
    for p in (1e-3, 1e-2, 0.1):
        for n in ALVOS:
            d, pl = d_minimo(p, n)
            print("  p=%.0e  N=%.1e  ->  d=%d (p_L=%.2e)  qubits fisicos = N*d = %.2e" % (p, n, d, pl, n * d))
    print()
    print("MENCIONA: este codigo corrige so erro de fase (o ruido de tpl_sincronia).")
    print("  Qubits reais sofrem fase e bit; para os dois: codigo de superficie (~2d^2 por qubit logico).")
    print("  Supoe medicao de sindrome perfeita e erros independentes; com medicao ruidosa o limiar real e menor.")
    if any(d is None for d in res.values()):
        raise SystemExit("ESTADO: BLOQUEADO (mutante sobreviveu)")
    print()
    print("ESTADO: PROVADO o limiar para ruido de fase; MOSTRADA a escala; MENCIONADO o erro de bit")


if __name__ == "__main__":
    main() if sys.argv[1:] == [] else SystemExit("uso: python3 tpl_correcao.py")
