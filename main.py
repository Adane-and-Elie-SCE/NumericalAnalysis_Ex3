# Adane Adego 315721969
# Elie Bracha 204795900

# https://github.com/Adane-and-Elie-SCE/NumericalAnalysis_Ex3

def determinant_calc(mat):
    if len(mat) == 2:
        ans = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]
        return ans

    minor = [[0 for i in range(len(mat) - 1)] for j in range(len(mat) - 1)]
    determinant = 0

    for k in range(len(mat)):
        i, j = 0, 0
        while i < len(mat):
            if i != k:
                minor[j] = mat[i][1:]
                j += 1
            i += 1
        determinant += ((-1) ** k) * mat[k][0] * determinant_calc(minor)
    return determinant


def elementary_delete(mat, tag):
    i, j = tag
    ans = [[int(i == j) for j in range(len(mat))] for i in range(len(mat))]
    ans[i][j] = -1 * (mat[i][j] / mat[j][j])
    return ans


def elementary_switch_rows(mat, tag):
    ans = [[int(i == j) for j in range(len(mat))] for i in range(len(mat))]
    ans[tag[0]], ans[tag[1]] = ans[tag[1]], ans[tag[0]]
    return ans


def elementary_mul_row(mat, index, scalar):
    e = [[int(i == j) for j in range(len(mat))] for i in range(len(mat))]
    mat[index][index] = scalar
    return mat


def matrix_pivoting(mat, index):
    size = len(mat)
    max_value = mat[index][index]
    max_index = index

    for i in range(index + 1, size):
        if mat[i][index] > max_value:
            max_value = mat[i][index]
            max_index = i

    return matrix_mul(elementary_switch_rows(mat, [index, max_index]), mat)


def fix_approximation(mat):
    epsilon = 2 ** (-26)
    for i in range(len(mat)):
        for j in range(len(mat)):
            if abs(mat[i][j]) < epsilon:
                mat[i][j] = 0
    return mat


def matrix_mul(a, b):
    size = len(a)
    ans = [[0 for i in range(size)] for j in range(size)]
    for i in range(size):
        for j in range(size):
            for k in range(size):
                ans[i][j] += a[i][k] * b[k][j]
    return fix_approximation(ans)


def matrix_vector_mul(mat, b):
    size = len(mat)
    ans = [0 for i in range(size)]

    for i in range(size):
        for j in range(size):
            ans[i] += mat[i][j] * b[j]
    return ans


def print_matrix(mat):
    for row in mat:
        print(row)
    print()


def print_vector(v):
    size = len(v)
    for i in range(size):
        print('[' + str(v[i]) + ']')


def inverse_by_gauss(mat):
    if not determinant_calc(mat):
        print('no inverse')

    size = len(mat)
    temp = mat
    ans = [[int(i == j) for j in range(size)] for i in range(size)]

    # below diagonal
    for j in range(0, size - 1):

        for i in range(j + 1, size):
            e = elementary_delete(temp, [i, j])
            temp = matrix_mul(e, temp)
            ans = matrix_mul(e, ans)

    # above diagonal
    for j in range(size - 1, -1, -1):

        for i in range(j - 1, -1, -1):
            e = elementary_delete(temp, [i, j])
            temp = matrix_mul(e, temp)
            ans = matrix_mul(e, ans)

    # final step
    for i in range(size):
        for j in range(size):
            ans[i][j] /= temp[i][i]

    return ans


def lu_decomposition(mat):
    size = len(mat)
    u = mat
    l = [[int(i == j) for j in range(size)] for i in range(size)]

    # below diagonal
    for j in range(0, size):
        for i in range(j + 1, size):
            e = elementary_delete(u, [i, j])
            u = matrix_mul(e, u)
            e[i][j] = -1 * e[i][j]
            l = matrix_mul(l, e)

    print("Triangle Matrix L:")
    print_matrix(l)
    print("Triangle Matrix U:")
    print_matrix(u)


# --------------------------------------------------------------


# D matrix
def diagonal_matrix(mat):
    diagonal_mat = [[0 for i in range(len(mat))] for j in range(len(mat))]

    for i in range(len(mat)):
        for j in range(len(mat)):
            diagonal_mat[i][j] = mat[i][j]

    for i in range(len(mat)):
        for j in range(len(mat)):
            if i != j:
                diagonal_mat[i][j] = 0
    return diagonal_mat


# L matrix
def low_matrix(mat):
    low_mat = [[0 for i in range(len(mat))] for j in range(len(mat))]

    for i in range(len(mat)):
        for j in range(len(mat)):
            low_mat[i][j] = mat[i][j]

    for i in range(len(mat)):
        for j in range(len(mat)):
            if i <= j:
                low_mat[i][j] = 0
    return low_mat


# U matrix
def upper_matrix(mat):
    upper_mat = [[0 for i in range(len(mat))] for j in range(len(mat))]

    for i in range(len(mat)):
        for j in range(len(mat)):
            upper_mat[i][j] = mat[i][j]

    for i in range(len(mat)):
        for j in range(len(mat)):
            if i >= j:
                upper_mat[i][j] = 0
    return upper_mat


# add between two matrix
def matrix_add(mat1, mat2):
    mat3 = [[0 for i in range(len(mat1))] for j in range(len(mat2))]
    for i in range(len(mat3)):
        for j in range(len(mat3)):
            mat3[i][j] = mat1[i][j] + mat2[i][j]
    return mat3


def mul_scalar_matrix(scalar, mat):
    mat1 = [[0 for i in range(len(mat))] for j in range(len(mat))]
    for i in range(len(mat)):
        for j in range(len(mat)):
            mat1[i][j] = mat[i][j] * scalar
    return mat1


# G with Jacobian
def G_Jacobian(mat):
    L = low_matrix(mat)
    U = upper_matrix(mat)
    D = diagonal_matrix(mat)
    LU = matrix_add(L, U)
    LU = pivoted_matrix(LU)
    print(determinant_calc(LU))
    G = matrix_mul(mul_scalar_matrix(-1, inverse_by_gauss(D)), inverse_by_gauss(LU))
    return G


def pivoted_matrix(mat):
    mat1 = [[0 for i in range(len(mat))] for j in range(len(mat))]
    for i in range(len(mat)):
        for j in range(len(mat)):
            mat1[i][j] = mat[i][j]

    for i in range(len(mat)):
        mat1 = matrix_pivoting(mat1, i)
    return mat1


def H_Jacobian(mat):
    H = inverse_by_gauss(mat)
    return H


# If the matrix have a dominant diagonal (tnai maspik)
def dominant_diagonal(mat):
    for k in range(len(mat)):
        matrix_pivoting(mat, k)

    for i in range(len(mat)):
        sum = 0
        pivot = abs(mat[i][i])
        for j in range(len(mat)):
            if j != i:
                sum += abs(mat[i][j])
        if sum >= pivot:
            return False
    return True


def print_list(z, y, x):
    print("[ Zr+1 , Yr+1 , Xr+1 , Number of interaction ]")
    for i in range(len(x)):
        print("[ {0} , {1} , {2} , {3} ]".format(z[i], y[i], x[i], i))


def Jacobian_method(mat, b):
    epsilon = 0.00001
    X, Y, Z = [0], [0], [0]
    i = 0
    while True:
        X1 = round(abs((b[0] - mat[0][1] * Y[i] - mat[0][2] * Z[i]) / mat[0][0]), 6)
        Y1 = round(abs((b[1] - mat[1][0] * X[i] - mat[1][2] * Z[i]) / mat[1][1]), 6)
        Z1 = round(abs((b[2] - mat[2][0] * X[i] - mat[2][1] * Y[i]) / mat[2][2]), 6)
        X.append(X1)
        Y.append(Y1)
        Z.append(Z1)
        if epsilon > abs(X1 - X[i]):
            print("----------------Jacobian Method----------------")
            print_list(Z, Y, X)
            break
        i += 1


def Gauss_Seidel_method(mat, b):
    epsilon = 0.00001
    X, Y, Z = [0], [0], [0]
    i = 0
    while True:
        X1 = round(abs((b[0] - mat[0][1] * Y[i] - mat[0][2] * Z[i]) / mat[0][0]), 6)
        Y1 = round(abs((b[1] - mat[1][0] * X1 - mat[1][2] * Z[i]) / mat[1][1]), 6)
        Z1 = round(abs((b[2] - mat[2][0] * X1 - mat[2][1] * Y1) / mat[2][2]), 6)
        X.append(X1)
        Y.append(Y1)
        Z.append(Z1)
        if epsilon > abs(X1 - X[i]):
            print("--------------Gauss_Seidel Method--------------")
            print_list(Z, Y, X)
            break
        i += 1


def Wrong_Jacobian_method(mat, b):
    epsilon = 0.00001
    X, Y, Z = [0], [0], [0]
    i = 0
    while True:
        X1 = round(abs((b[0] - mat[0][1] * Y[i] - mat[0][2] * Z[i]) / mat[0][0]), 6)
        Y1 = round(abs((b[1] - mat[1][0] * X[i] - mat[1][2] * Z[i]) / mat[1][1]), 6)
        Z1 = round(abs((b[2] - mat[2][0] * X[i] - mat[2][1] * Y[i]) / mat[2][2]), 6)
        X.append(X1)
        Y.append(Y1)
        Z.append(Z1)
        if i > 50:
            print("----------------Jacobian Method----------------")
            print_list(Z, Y, X)
            break
        i += 1


def Wrong_Gauss_Seidel_method(mat, b):
    epsilon = 0.00001
    X, Y, Z = [0], [0], [0]
    i = 0
    while True:
        X1 = round(abs((b[0] - mat[0][1] * Y[i] - mat[0][2] * Z[i]) / mat[0][0]), 6)
        Y1 = round(abs((b[1] - mat[1][0] * X1 - mat[1][2] * Z[i]) / mat[1][1]), 6)
        Z1 = round(abs((b[2] - mat[2][0] * X1 - mat[2][1] * Y1) / mat[2][2]), 6)
        X.append(X1)
        Y.append(Y1)
        Z.append(Z1)
        if i > 50:
            print("--------------Gauss_Seidel Method--------------")
            print_list(Z, Y, X)
            break
        i += 1


def Driver(mat, b):
    if dominant_diagonal(mat):
        print("This matrix has Diagonal Dominant")
        Jacobian_method(mat, b)
        print("\n\n\n\n")
        Gauss_Seidel_method(mat, b)
    else:
        print("This matrix don't has Diagonal Dominant but \nwe will do both method (Jacobian and Gauss-Seidel)")
        Jacobian_method(mat, b)
        print("\n\n\n\n")
        Gauss_Seidel_method(mat, b)


A = [[4, 2, 0], [2, 10, 4], [0, 4, 5]]

B =[[2, 10, 4], [4,2,0], [0, 4, 5]]

b = [2, 6, 5]

Driver(B, b)
