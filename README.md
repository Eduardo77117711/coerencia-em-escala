# coerencia-em-escala

Quatro programas curtos que respondem, com prova executável, a uma pergunta só: **quanto da coerência quântica sobrevive quando a escala cresce — e o que a correção de erros muda nisso.**

Cada resultado é marcado pelo que o código realmente faz:
**PROVA** (calculado e conferido), **MOSTRA** (fórmula fechada conferida), **MENCIONA** (dito, não simulado).

| Arquivo | Pergunta | Resultado |
|---|---|---|
| `decoerencia.py` | Um fóton, dois caminhos, um ambiente que "anota" o caminho. Quanto da franja de interferência fica? | Visibilidade = \|⟨e_A\|e_B⟩\| = cos θ. Saber o caminho apaga a franja. |
| `tpl_sincronia.py` | E com *n* qubits em sincronia total (GHZ), cada um com seu ambiente? | Coerência = cᴺ, provada com matriz densidade até n = 8; erro máximo por qubit para manter 50% até n = 10⁸⁰. |
| `tpl_correcao.py` | A correção de erros vence a escala? | Limiar provado em p = ½ para ruído de fase (enumeração exata dos 2ᵈ padrões, d até 15, mais Monte Carlo); distância mínima para trilhões de qubits lógicos. |
| `tpl_superficie.py` | E com erro de bit **e** de fase, com a própria medição errando? | Código de superfície (stim + pymatching), d = 3, 5, 7: abaixo do limiar p_L cai com d; limiar medido ≈ 0,6–0,7%; mutantes do decodificador derrubados. |

Os programas se encadeiam: `tpl_sincronia.py` só roda se `decoerencia.py` tiver exatamente o hash fixado (pin `ef378ac6…`). Cada um traz **mutantes** — versões erradas de propósito — e só declara estado se todos forem derrubados.

## Rodar

```bash
pip install numpy stim pymatching
python3 -I decoerencia.py autoteste
python3 -I tpl_sincronia.py
python3 -I tpl_correcao.py
python3 -I tpl_superficie.py
```

Os quatro rodam no GitHub Actions a cada envio.

## Limites (ditos pelo próprio código)

- `tpl_superficie.py` usa o ruído de circuito padrão do stim: erros independentes, sem vazamento nem erros correlacionados; a escala para trilhões é extrapolação.
- `tpl_correcao.py` supõe medição de síndrome perfeita; `tpl_superficie.py` remove essa suposição.
- A escala do código de superfície em `tpl_sincronia.py` é heurística, não simulada.

---
Parte do ecossistema [TPL](https://github.com/Eduardo77117711/-tpl-).

Criado por Carlos Eduardo Acevedo Rodrigues. Copyright 2026.
Licenciado sob a [GNU AGPL-3.0](LICENSE): livre para usar, estudar, modificar e redistribuir. Toda versão derivada, inclusive a oferecida como serviço pela internet, continua livre sob a mesma licença e deve manter esta atribuição e o arquivo [NOTICE](NOTICE).
