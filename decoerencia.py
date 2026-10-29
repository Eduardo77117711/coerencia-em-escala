"""decoerencia.py - interferometro com marcador de caminho (so biblioteca padrao).
Foton em dois caminhos; o ambiente registra o caminho com fidelidade s = sin(theta).
Estado: (|A>|e_A> + |B>|e_B>)/sqrt(2), com <e_A|e_B> = cos(theta).
Na saida, com fase phi entre os caminhos: P(phi) = (1 + cos(theta) cos(phi)) / 2."""
import math, sys


def amplitudes(theta, phi):
    # componentes da saida 0 para cada estado do ambiente (base |0>, |1> do ambiente)
    eA = (1.0, 0.0)
    eB = (math.cos(theta), math.sin(theta))
    r = 1 / math.sqrt(2)
    # saida 0 recebe (caminho A + e^{i phi} caminho B)/sqrt(2)
    return [complex(r * r * eA[k]) + r * r * complex(math.cos(phi), math.sin(phi)) * eB[k] for k in range(2)]


def prob_saida(theta, phi):
    return sum(abs(a) ** 2 for a in amplitudes(theta, phi))


def visibilidade(theta, passos=721):
    ps = [prob_saida(theta, 2 * math.pi * i / (passos - 1)) for i in range(passos)]
    return (max(ps) - min(ps)) / (max(ps) + min(ps))


if sys.argv[1:] == ["autoteste"]:
    def ok(c, m):
        if not c:
            raise SystemExit("FALHA " + m)
    e = 1e-9
    for theta in (0.0, 0.3, 0.7, 1.0, math.pi / 4, 1.3, math.pi / 2):
        ok(abs(visibilidade(theta) - abs(math.cos(theta))) < 1e-6, "visibilidade = |<eA|eB>| em %.2f" % theta)
        for phi in (0.0, 0.5, 2.0, math.pi):
            ok(abs(prob_saida(theta, phi) - (1 + math.cos(theta) * math.cos(phi)) / 2) < e, "formula")
    ok(abs(prob_saida(0.0, 0.0) - 1) < e and abs(prob_saida(0.0, math.pi)) < e, "onda pura: franja total")
    ok(abs(prob_saida(math.pi / 2, 0.0) - 0.5) < e and abs(prob_saida(math.pi / 2, math.pi) - 0.5) < e, "caminho conhecido: sem franja")
    print("autoteste: OK")
