"""Network topology matrices with numpy. No GUI code here."""
from dataclasses import dataclass
import numpy as np


def branch_names(n):
    letters = "abcdefghijklmnopqrstuvwxyz"
    return [letters[i] if i < 26 else f"b{i + 1}" for i in range(n)]


def _clean(text):
    for ch in "[],;":
        text = text.replace(ch, " ")
    return text


def parse_matrix(text):
    rows = []
    for line in text.strip().splitlines():
        line = _clean(line)
        if line.strip():
            rows.append([float(x) for x in line.split()])
    if not rows:
        raise ValueError("Matrix is empty")
    if len({len(r) for r in rows}) != 1:
        raise ValueError("Rows have different lengths")
    return np.array(rows)


def parse_vector(text, n, name):
    text = _clean(text).strip()
    if not text:
        return None
    v = np.array([float(x) for x in text.split()])
    if v.size != n:
        raise ValueError(f"{name} needs {n} values, got {v.size}")
    return v


def spanning_tree(src, dst, n_nodes, order):
    """Greedy union-find. Returns branch indices that form a forest."""
    parent = list(range(n_nodes))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    tree = []
    for k in order:
        ra, rb = find(src[k]), find(dst[k])
        if ra != rb:
            parent[ra] = rb
            tree.append(k)
    return tree


@dataclass
class Result:
    N: int
    b: int
    names: list
    src: np.ndarray
    dst: np.ndarray
    tree: list
    links: list
    perm: list
    A: np.ndarray   # reduced incidence, tree columns first
    B: np.ndarray   # tie-set, tree columns first
    C: np.ndarray   # cut-set, tree columns first
    struct_err: float
    J: np.ndarray = None   # branch currents, original order
    V: np.ndarray = None   # branch voltages, original order
    kcl_err: float = None
    kvl_err: float = None
    cutset_agrees: bool = None


def analyze(A, tree_names="", Z=None, Vs=None, Is=None):
    A = np.asarray(A, float)
    if not np.isin(A, [-1, 0, 1]).all():
        raise ValueError("A may only contain -1, 0, 1")
    n, b = A.shape
    pos, neg = A == 1, A == -1
    if (pos.sum(0) > 1).any() or (neg.sum(0) > 1).any():
        raise ValueError("Each column may have at most one +1 and one -1")
    if (~(pos.any(0) | neg.any(0))).any():
        raise ValueError("A column is all zeros")

    full = bool((pos.any(0) & neg.any(0)).all())
    N = n if full else n + 1
    Ar = A[:-1] if full else A
    src = np.where(pos.any(0), pos.argmax(0), N - 1)  # +1 = current leaves
    dst = np.where(neg.any(0), neg.argmax(0), N - 1)  # -1 = current enters
    names = branch_names(b)

    toks = tree_names.replace(",", " ").split()
    if toks:
        bad = [t for t in toks if t not in names]
        if bad:
            raise ValueError(f"Unknown branch names: {bad}")
        given = [names.index(t) for t in toks]
        tree = spanning_tree(src, dst, N, given)
        if sorted(tree) != sorted(given) or len(tree) != N - 1:
            raise ValueError("Given branches do not form a spanning tree")
    else:
        tree = spanning_tree(src, dst, N, range(b))
        if len(tree) != N - 1:
            raise ValueError("Graph is not connected")

    links = [k for k in range(b) if k not in tree]
    perm = tree + links
    Ap = Ar[:, perm]
    AT, AL = Ap[:, :N - 1], Ap[:, N - 1:]

    CL = np.linalg.solve(AT, AL)                 # C_L = A_T^-1 A_L
    C = np.hstack([np.eye(N - 1), CL])           # C = A_T^-1 A
    B = np.hstack([-CL.T, np.eye(len(links))])   # B_T = -C_L^T
    struct_err = float(np.abs(Ap @ B.T).max()) if len(links) else 0.0

    res = Result(N, b, names, src, dst, tree, links, perm, Ap, B, C, struct_err)

    if Z is not None:
        Vs = np.zeros(b) if Vs is None else Vs
        Is = np.zeros(b) if Is is None else Is
        Zp, Vp, Ip = Z[perm], Vs[perm], Is[perm]
        Zm = np.diag(Zp)
        # loop method: B Z B^T IL = B Vs - B Z Is
        if len(links):
            IL = np.linalg.solve(B @ Zm @ B.T, B @ Vp - B @ Zm @ Ip)
        else:
            IL = np.zeros(0)
        J = B.T @ IL
        V = Zm @ (J + Ip) - Vp
        res.kcl_err = float(np.abs(Ap @ J).max())
        res.kvl_err = float(np.abs(B @ V).max()) if len(links) else 0.0
        # cut-set method: C Y C^T VT = C Is - C Y Vs (needs all Z nonzero)
        if (Zp != 0).all():
            Y = np.diag(1 / Zp)
            VT = np.linalg.solve(C @ Y @ C.T, C @ Ip - C @ Y @ Vp)
            V2 = C.T @ VT
            J2 = Y @ (V2 + Vp) - Ip
            res.cutset_agrees = bool(np.allclose(V, V2) and np.allclose(J, J2))
        Jo, Vo = np.empty(b), np.empty(b)
        Jo[perm], Vo[perm] = J, V
        res.J, res.V = Jo, Vo
    return res
