"""tpl_sincronia.py v1.1 - TPL: quanto da sincronia sobrevive quando a escala cresce.
v1.1 (pai f46e2c58...309e): mutacao real do canal; pin da dependencia;
"qubits fisicos" em vez de "atomos"; limiar so como escala heuristica.

Base: decoerencia.py (Drive). Para 1 qubit, o ambiente que registra o caminho
com sobreposicao c = cos(theta) deixa visibilidade V = c.
Pergunta de escala: n qubits em sincronia total (GHZ), cada um com seu ambiente.

PROVA    - matriz densidade 2^n x 2^n, canal de Kraus aplicado qubit a qubit;
           coerencia medida C(n) = 2 rho[0..0, 1..1] (com sinal) comparada com c^n.
           Para n=1 deve coincidir com a visibilidade do decoerencia.py.
MOSTRA   - formula fechada C(n) = c^n; para manter C >= alvo, cada qubit precisa
           de erro eps <= -ln(alvo)/n (aprox. 1a ordem), conferido de forma exata.
MENCIONA - trilhoes de qubits e a escala de 10^80; e a brecha conhecida:
           correcao de erros, onde a coerencia logica cresce com a escala
           quando o erro fisico fica abaixo de um limiar (escala heuristica
           p_L ~ A (p/p_th)^((d+1)/2), codigo de superficie, d impar).
Sem booleanos de decisao: o resultado e um estado textual.
"""
import math
import sys

import numpy as np

import hashlib
import os

AQUI = os.path.dirname(os.path.abspath(__file__))
PIN_DECOERENCIA = "ef378ac634e434f0088d31554cfbd06a923e28dcebf9fca592ea7795fc66e9d4"
try:
    with open(os.path.join(AQUI, "decoerencia.py"), "rb") as f:
        _dep = hashlib.sha256(f.read()).hexdigest()
except FileNotFoundError:
    raise SystemExit("BLOQUEADO: decoerencia.py ausente (item n=1 NAO EXECUTADO)") from None
if _dep != PIN_DECOERENCIA:
    raise SystemExit("BLOQUEADO: decoerencia.py fora do pin (dependencia nao verificada)")
sys.path.insert(0, AQUI)
import decoerencia  # noqa: E402

TOL = 1e-10
N_MAX = 8
ALVO = 0.5
ESCALAS = (1e3, 1e6, 1e9, 2.5e12, 5e12, 2e13, 1e80)


def ghz_rho(n):
    v = np.zeros(2 ** n)
    v[0] = v[-1] = 1 / math.sqrt(2)
    return np.outer(v, v).astype(np.complex128)


def canal_ambiente(rho, n, i, c):
    """Kraus: K0 = sqrt((1+c)/2) I, K1 = sqrt((1-c)/2) Z no qubit i."""
    op = _op(n, i, Z)
    return (1 + c) / 2 * rho + (1 - c) / 2 * (op @ rho @ op)


def coerencia(rho):
    """Elemento com sinal e fase: 2*rho[0..0,1..1]. O modulo sozinho nao ve um Z global."""
    return complex(2 * rho[0, -1])


def medir(n, c, canal=None):
    canal = canal or canal_ambiente
    rho = ghz_rho(n)
    for i in range(n):
        rho = canal(rho, n, i, c)
    if abs(np.trace(rho).real - 1) > TOL:
        raise SystemExit("DIVERGENTE: traco nao conservado")
    return coerencia(rho)


def eps_max(n, alvo):
    """Erro maximo por qubit (eps = 1 - c) para coerencia >= alvo. Exato."""
    return -math.expm1(math.log(alvo) / n)


def prova(canal=None):
    """Devolve o primeiro desvio encontrado, ou None se tudo confere."""
    cs = (1.0, 0.99, 0.9, math.cos(0.7), 0.5, 0.0)
    for c in cs:
        v1 = decoerencia.visibilidade(math.acos(c))
        if abs(abs(medir(1, c, canal)) - v1) > 1e-6:
            return "n=1 nao coincide com decoerencia.py (c=%.3f)" % c
        for n in range(1, N_MAX + 1):
            if abs(medir(n, c, canal) - c ** n) > TOL:
                return "n=%d c=%.3f" % (n, c)
    return None


def _op(n, i, m):
    op = np.array([[1.0]], dtype=np.complex128)
    for j in range(n):
        op = np.kron(op, m if j == i else np.eye(2))
    return op


Z = np.array([[1, 0], [0, -1]], dtype=np.complex128)
X = np.array([[0, 1], [1, 0]], dtype=np.complex128)
MUTANTES = {
    "peso trocado": lambda r, n, i, c: (1 - c) / 2 * r + (1 + c) / 2 * (_op(n, i, Z) @ r @ _op(n, i, Z)),
    "ultimo qubit esquecido": lambda r, n, i, c: r if i == n - 1 else canal_ambiente(r, n, i, c),
    "ruido X no lugar de Z": lambda r, n, i, c: (1 + c) / 2 * r + (1 - c) / 2 * (_op(n, i, X) @ r @ _op(n, i, X)),
    "c ao quadrado": lambda r, n, i, c: canal_ambiente(r, n, i, c * c),
}


def mutacao():
    """Mutacao real: o canal e trocado e a prova inteira roda de novo."""
    return {nome: prova(m) for nome, m in MUTANTES.items()}


def limiar(p, p_th=1e-2, d=25, a=0.1):
    return a * (p / p_th) ** ((d + 1) / 2)


def main():
    desvio = prova()
    if desvio is not None:
        raise SystemExit("DIVERGENTE: " + desvio)
    res = mutacao()
    m = sum(1 for d in res.values() if d is not None)
    print("PROVA: matriz densidade n=1..%d, 6 forcas de ambiente; coerencia medida = c^n (tol %g)" % (N_MAX, TOL))
    print("PROVA: n=1 coincide com decoerencia.py (pin %s...)" % PIN_DECOERENCIA[:8])
    for nome, d in res.items():
        print("PROVA: mutante do canal '%s' -> %s" % (nome, "MORTO em " + d if d else "SOBREVIVEU"))
    print()
    print("MOSTRA: para manter %.0f%% da sincronia, erro maximo por qubit:" % (ALVO * 100))
    for n in ESCALAS:
        e = eps_max(n, ALVO)
        print("  n = %.1e  ->  eps <= %.3e   ln C calculado = %.12f (alvo %.12f)" % (n, e, n * math.log1p(-e), math.log(ALVO)))
    print("  (ponto flutuante de 64 bits; erro relativo ~1e-16 por operacao)")
    print()
    print("MENCIONA (heuristica, nao simulada): codigo de superficie, d=25, p_th=1e-2, A=0.1")
    for p in (5e-3, 1e-3, 1e-4):
        pl = limiar(p)
        print("  erro fisico %.0e -> erro logico ~%.1e -> escala de qubits logicos ~%.1e (se erros independentes)" % (p, pl, math.log(2) / pl))
    print("  cada qubit logico d=25 usa ~2*d^2 = %d qubits fisicos (estimativa de recursos)" % (2 * 25 ** 2))
    print("  erro fisico atual: MENCIONADO sem fonte neste codigo")
    if m != len(MUTANTES):
        raise SystemExit("ESTADO: BLOQUEADO")
    print()
    print("ESTADO: PROVADO ate n=%d; MOSTRADO c^n para todo n; MENCIONADO o limiar" % N_MAX)


if __name__ == "__main__":
    main()
