#Adane Adego 315721969
#Elie Bracha 204795900

def matrix_pivoting(mat, index):
    size = len(mat)
    max_value = mat[index][index]
    max_index = index

    for i in range(index + 1, size):
        if mat[i][index] > max_value:
            max_value = mat[i][index]
            max_index = i

    return matrix_mul(elementary_switch_rows(mat, [index, max_index]), mat)


def elementary_switch_rows(mat, tag):
    ans = [[int(i == j) for j in range(len(mat))] for i in range(len(mat))]
    ans[tag[0]], ans[tag[1]] = ans[tag[1]], ans[tag[0]]
    return ans


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


def print_matrix(mat):
    for row in mat:
        print(row)
    print()


# --------------------------------------------------------------

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
    print("Before pivoting")
    print(mat)
    for k in range(len(mat)):
        matrix_pivoting(mat, k)
    mat2 = mat
    print("After pivoting")
    print(mat2)

    if dominant_diagonal(mat2):
        print("This matrix has Diagonal Dominant")
        Jacobian_method(mat2, b)
        print("\n\n\n\n")
        Gauss_Seidel_method(mat2, b)
    else:
        print("This matrix don't has Diagonal Dominant but \nwe will do both method (Jacobian and Gauss-Seidel)")
        Wrong_Jacobian_method(mat2, b)
        print("\n\n\n\n")
        Wrong_Gauss_Seidel_method(mat2, b)


A = [[4, 2, 0], [2, 10, 4], [0, 4, 5]]
b = [2, 6, 5]

Driver(A, b)
