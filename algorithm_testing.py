
def find_inverse(element : int):
    for i in range(0, 256):
        raw = element * i
        result = raw & 255
        if (result == 1):
            return i
    return 0

def gf_2_mul_mod(a : int, b : int, q : int):
    result : int = 0
    for i in range(0, 8):
        if (a >> i) & 1 == 1:
            result = result ^ b
        b = b << 1
        if (b >> 8) & 1 == 1:
            b = b ^ q
    return result

def gf_2_square_mod(a : int, q : int):
    result : int = 0
    return result

def gf_2_inverse_it(a : int, q :int):
    temp = gf_2_mul_mod(a, a, q)
    temp = gf_2_mul_mod(temp, a, q)
    temp = gf_2_mul_mod(temp, temp, q)

    arm1 = gf_2_mul_mod(temp, a, q)
    ar = gf_2_mul_mod(arm1, a, q)
    ar_inverse = ar

    return arm1#gf_2_mul_mod(ar_inverse, arm1, q)

def gf_2_create_exp_table(generator : int, irreducible_polynomial : int):
    table = []
    t : int = 1
    for i in range(0, 255):
        table.append(t)
        t = gf_2_mul_mod(generator, t, irreducible_polynomial)
    return table

aes_irreducible_polynomial = 0b100011011
ccsds_irreducible_polynomial = 0b110000111
#print(gf_2_create_exp_table(0b00000011, aes_irreducible_polynomial))
#print(gf_2_create_exp_table(0b00000011, 0b100011101))

test_inverse = gf_2_inverse_it(57, aes_irreducible_polynomial)
test_result = gf_2_mul_mod(57, test_inverse, aes_irreducible_polynomial)
print(test_inverse)
print(test_result)

ref_square = gf_2_mul_mod(57, 57, aes_irreducible_polynomial)
test_square = 0

def forney_algorithm(msg, syndromes, roots, locator, evaluator):
    error_mag = []
    error_pos = []
    for i in range(0, len(roots)):
        error_pos = len(msg - roots[i])
        roots[i]
        locator[i]
        error_mag = 0
        for j in range(0, len(locator)):
            return
    return

def gf_2_polynomial_mod(dividend : int, divisor : int):
    result : int = dividend
    for i in range(8, 0):
        if (dividend >> i) & 1 == 1:
            dividend = dividend ^ (divisor << i)
        return 0

def calc_reduction_matrix(q : int):
    matrix = [[0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]]
    for j in range(0, 8):
        matrix[j][0] = (q >> j) & 1

    for j in range(0, 8):
        for i in range(1, 7):
            if j - 1 >= 0:
                matrix[j][i] = matrix[j - 1][i - 1] ^ matrix[7][i - 1]
            else:
                matrix[j][i] = matrix[7][i - 1]
    return matrix

def collapse_matrix(matrix):
    collapsed = []
    for j in range(0, len(matrix)):
        row : int = 0
        for i in range(0, len(matrix[0])):
            row |= matrix[j][i] << i
        collapsed.append(row)
    return collapsed

multiplication_matrix = [
    [ 0, -1, -1, -1, -1, -1, -1, -1],
    [ 1,  0, -1, -1, -1, -1, -1, -1],
    [ 2,  1,  0, -1, -1, -1, -1, -1],
    [ 3,  2,  1,  0, -1, -1, -1, -1],
    [ 4,  3,  2,  1,  0, -1, -1, -1],
    [ 5,  4,  3,  2,  1,  0, -1, -1],
    [ 6,  5,  4,  3,  2,  1,  0, -1],
    [ 7,  6,  5,  4,  3,  2,  1,  0],
    [-1,  7,  6,  5,  4,  3,  2,  1],
    [-1, -1,  7,  6,  5,  4,  3,  2],
    [-1, -1, -1,  7,  6,  5,  4,  3],
    [-1, -1, -1, -1,  7,  6,  5,  4],
    [-1, -1, -1, -1, -1,  7,  6,  5],
    [-1, -1, -1, -1, -1, -1,  7,  6],
    [-1, -1, -1, -1, -1, -1, -1,  7]]
def gf_2_mul_mastrovito(a : int, b : int, reduction_matrix):
    d : int = 0
    #for i in range(0, len(multiplication_matrix)):
    #    for j in range(0, 8):
    #        if multiplication_matrix[i][j] != -1:
    #            val = ((a >> multiplication_matrix[i][j]) & 1) & ((b >> j) & 1)
    #            d ^= val << i
    for i in range(0, 8):
        if (a >> i) & 1 == 1:
            d = d ^ (b << i)
    
    c : int = d & 0xFF
    for j in range(0, 8):
        for i in range(0, 7):
            bit = ((d >> (i+8)) & 1)
            val = bit & reduction_matrix[j][i]
            c ^= val << j
    return c

def gf_2_create_test_table(generator : int, reduction_matrix):
    table = []
    t : int = 1
    for i in range(0, 255):
        table.append(t)
        t = gf_2_mul_mastrovito(generator, t, reduction_matrix)
    return table

aes_matrix = calc_reduction_matrix(aes_irreducible_polynomial)
ccsds_matrix = calc_reduction_matrix(ccsds_irreducible_polynomial)

def print_matrix(matrix):
    for row in matrix:
        print(row)
    print('\n')

print_matrix(aes_matrix)
print_matrix(ccsds_matrix)
print_matrix(calc_reduction_matrix(0b100101011))
print_matrix(calc_reduction_matrix(0b101011111))
print_matrix(calc_reduction_matrix(0b111110101))

#print(gf_2_create_test_table(0b00000011, aes_matrix))