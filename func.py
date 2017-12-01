from sympy import solve, Symbol, eye, pprint


def diagonalize(mat):
    print('A =')
    pprint(mat)
    print()

    eigvals = get_eigvals(mat)
    for i in range(0, len(eigvals)):
        print('λ' + str(i + 1) + ' =')
        pprint(eigvals[i])
        print()

    eigvecs = get_eigvecs(mat, eigvals)
    for i in range(0, len(eigvecs)):
        print('K' + str(i + 1) + ' =')
        pprint(eigvecs[i])
        print()

    p_mat = get_p_mat(mat, eigvecs)
    print('P =')
    pprint(p_mat)
    print()

    print('P^-1 =')
    pprint(p_mat ** -1)
    print()

    diag_mat = (p_mat ** -1) * mat * p_mat
    print('D = P^-1 * A * P =')
    pprint(diag_mat)
    print()


def get_eigvals(mat):
    r = mat.rows
    lam = Symbol('lamda', real=True)
    det_mat = mat - lam * eye(r)
    ret = solve(det_mat.det())
    return ret


def get_eigvecs(mat, eigvals):
    r = mat.rows
    ret = []
    for eigval in eigvals:
        target_mat = mat - eigval * eye(r)
        nss = get_nullspace(target_mat)
        for i in range(0, len(nss)):
            ret.append(nss[i])
    return ret


def get_nullspace(mat):
    c = mat.cols
    rdd, pvs = mat.rref(simplify=False)

    variables = []
    for i in range(0, c):
        flg = True
        for j in range(0, len(pvs)):
            if pvs[j] == i:
                flg = False
                break
        if flg:
            variables.append(i)

    tmp = []
    for v in variables:
        can_result = [0] * c
        can_result[v] = 1
        for i in range(0, len(pvs)):
            for ty in pvs[(i + 1):] + (v,):
                can_result[pvs[i]] -= rdd[i, ty]
        tmp.append(can_result)

    ret = []
    for i in range(0, len(tmp)):
        ret.append(mat._new(mat.cols, 1, tmp[i]))

    return ret


def get_p_mat(mat, eigvecs):
    ret = eigvecs[0]
    for i in range(1, len(eigvecs)):
        ret = ret.row_join(eigvecs[i])
    return ret
