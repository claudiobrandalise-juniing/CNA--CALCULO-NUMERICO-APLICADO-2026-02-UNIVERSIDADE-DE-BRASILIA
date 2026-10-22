n = 5

# Matriz A e Vetor b corrigidos conforme o enunciado (pág. 3)
A = [
[ 6.0, 0.0, -1.0, 0.0, 0.0],
[-3.0, 3.0, 0.0, 0.0, 0.0],
[ 0.0, -1.0, 9.0, 0.0, 0.0],
[ 0.0, -1.0, -8.0, 11.0, -2.0],
[-3.0, -1.0, 0.0, 0.0, 4.0]
]

# Atenção para b[1] = 10 (Q01*c01 - Q12*... conforme balanços)
b = [50.0, 0.0, 160.0, 0.0, 0.0]

# Réplicas para checagem final do resíduo
A_orig = [row[:] for row in A]
b_orig = b[:]

mult_div_count = 0
add_sub_count = 0

# 1. Eliminação Progressiva
for k in range(n - 1):
    pivo = A[k][k]
    if abs(pivo) < 1e-12:
        raise ValueError(f"Pivô nulo/muito pequeno na posição ({k},{k}).")

    for i in range(k + 1, n):
        fator = A[i][k] / pivo
        mult_div_count += 1

        for j in range(k + 1, n):
            A[i][j] -= fator * A[k][j]
            mult_div_count += 1
            add_sub_count += 1

        b[i] -= fator * b[k]
        mult_div_count += 1
        add_sub_count += 1

        A[i][k] = 0.0

print(f"\nMatriz A após a etapa k={k+1}:")
for row in A:
    print([round(val, 4) for val in row])

# 2. Substituição Retroativa
c = [0.0] * n
c[n - 1] = b[n - 1] / A[n - 1][n - 1]
mult_div_count += 1

for i in range(n - 2, -1, -1):
    soma = b[i]
    for j in range(i + 1, n):
        soma -= A[i][j] * c[j]
        mult_div_count += 1
        add_sub_count += 1
    c[i] = soma / A[i][i]
    mult_div_count += 1

# 3. Resíduo
residuo = []
for i in range(n):
    val_calc = sum(A_orig[i][j] * c[j] for j in range(n))
    residuo.append(val_calc - b_orig[i])

norma_residuo = max(abs(r) for r in residuo)

# 4. Saída
print("\n" + "="*40)
print("RESULTADOS FINAIS:")
for idx, val in enumerate(c, start=1):
    print(f"c{idx} = {val:.4f} mg/m³")

print(f"\nResíduo Máximo (||r||_max): {norma_residuo:.2e}")
print("\nCONTAGEM DE FLOPs:")
print(f"Multiplicações / Divisões: {mult_div_count}")
print(f"Somas / Subtrações: {add_sub_count}")
print(f"Total de FLOPs: {mult_div_count + add_sub_count}")