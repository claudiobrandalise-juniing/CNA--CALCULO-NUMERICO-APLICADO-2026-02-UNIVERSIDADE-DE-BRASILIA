"""
executar_tudo.py
================
Roda os testes da biblioteca e TODOS os programas do repositório, regenerando
as pastas outputs/. Serve para conferir que tudo funciona depois de qualquer
mudança (reprodutibilidade).

    python executar_tudo.py

Leva cerca de 1 minuto. Cada programa roda num processo separado; no final
aparece um resumo com o tempo de cada um e se terminou sem erro.
"""

import subprocess
import sys
import time
from pathlib import Path

RAIZ = Path(__file__).resolve().parent

# ordem: primeiro os testes da biblioteca, depois as atividades na ordem do curso
PROGRAMAS = (
    ["testes/testar_tudo.py"]
    + sorted(str(p.relative_to(RAIZ)) for p in RAIZ.glob("Aulas_Praticas/AP*/main.py"))
    + sorted(str(p.relative_to(RAIZ)) for p in RAIZ.glob("PPC_*/main.py"))
    + ["APC_03/verificar_apc3.py"]
    + sorted(str(p.relative_to(RAIZ)) for p in RAIZ.glob("Exemplos_Processador/E*.py"))
    + ["Exemplos_Processador/modelo_em_branco.py"]
)


def main():
    resumo = []
    for programa in PROGRAMAS:
        inicio = time.perf_counter()
        processo = subprocess.run([sys.executable, programa], cwd=RAIZ, capture_output=True, text=True)
        duracao = time.perf_counter() - inicio
        ok = processo.returncode == 0
        resumo.append((programa, ok, duracao))
        print(f"[{'OK' if ok else 'ERRO'}] {programa}  ({duracao:.1f} s)")
        if not ok:
            print(processo.stderr[-2000:])
    falhas = [p for p, ok, _ in resumo if not ok]
    print(f"\n{len(resumo) - len(falhas)}/{len(resumo)} programas terminaram sem erro "
          f"(tempo total {sum(d for _, _, d in resumo):.0f} s).")
    sys.exit(1 if falhas else 0)


if __name__ == "__main__":
    main()
