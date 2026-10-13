"""Qt front end. Run: python main.py"""

import sys

import matplotlib
import numpy as np

matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFontDatabase
from PyQt6.QtWidgets import (
    QApplication,
    QComboBox,
    QLabel,
    QLineEdit,
    QMainWindow,
    QPlainTextEdit,
    QPushButton,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

from cadaid.cad_core import analyze, parse_matrix, parse_vector
from cadaid.cad_plot import draw_graph

EXAMPLES = {
    "Lecture 1 graph (full A)": dict(
        A="-1 1 0 -1 0 0\n0 -1 1 0 1 0\n1 0 -1 0 0 1\n0 0 0 1 -1 -1",
        tree="d e f",
        Z="",
        Vs="",
        Is="",
    ),
    "Lecture 2 practice (reduced A)": dict(
        A="-1 1 1 1 0\n0 0 0 -1 0\n1 -1 0 0 -1", tree="", Z="", Vs="", Is=""
    ),
    "Lecture 3 circuit (solve)": dict(
        A="-1 1 0 0 0\n0 -1 1 1 0\n0 0 0 -1 1",
        tree="",
        Z="0 5 5 10 5",
        Vs="10 0 0 0 0",
        Is="0 0 0 0 0",
    ),
}


def fmt(M, cols):
    M = np.round(np.asarray(M, float), 4) + 0.0
    cells = [[f"{x:g}" for x in r] for r in M]
    w = max([len(c) for r in cells for c in r] + [len(n) for n in cols]) + 2
    lines = ["".join(n.rjust(w) for n in cols)]
    lines += ["".join(c.rjust(w) for c in r) for r in cells]
    return "\n".join(lines)


def report(res):
    nm = res.names
    cols = [nm[k] for k in res.perm]
    out = [
        f"Nodes N = {res.N}, branches B = {res.b}",
        "Tree:  " + " ".join(nm[k] for k in res.tree),
        "Links: "
        + " ".join(nm[k] for k in res.links)
        + f"   (B - N + 1 = {res.b - res.N + 1})",
        "",
        "A (reduced, tree first):",
        fmt(res.A, cols),
        "",
        "B (tie-set):",
        fmt(res.B, cols),
        "",
        "C (cut-set):",
        fmt(res.C, cols),
        "",
        f"Check max|A B^T| = {res.struct_err:.2e}",
    ]
    if res.J is not None:
        out += [
            "",
            "Branch results (loop method):",
            "".join(s.rjust(10) for s in ["branch", "J", "V"]),
        ]
        for k in range(res.b):
            out.append(f"{nm[k]:>10}{res.J[k]:>10.4f}{res.V[k]:>10.4f}")
        out.append(f"Check max|A J| (KCL) = {res.kcl_err:.2e}")
        out.append(f"Check max|B V| (KVL) = {res.kvl_err:.2e}")
        if res.cutset_agrees is None:
            out.append("Cut-set method skipped (a Z is zero)")
        else:
            out.append(f"Cut-set method agrees: {res.cutset_agrees}")
    return "\n".join(out)


class Window(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("CAD network topology tool")
        mono = QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont)

        self.example = QComboBox()
        self.example.addItems(["Load example..."] + list(EXAMPLES))
        self.example.currentIndexChanged.connect(self.load_example)

        self.A = QPlainTextEdit()
        self.A.setFont(mono)
        self.tree = QLineEdit()
        self.tree.setPlaceholderText("auto, or e.g. d e f")
        self.Z = QLineEdit()
        self.Vs = QLineEdit()
        self.Is = QLineEdit()
        for e in (self.Z, self.Vs, self.Is):
            e.setPlaceholderText("optional, one value per branch")
        btn = QPushButton("Compute and plot")
        btn.clicked.connect(self.compute)

        left = QWidget()
        lay = QVBoxLayout(left)
        lay.addWidget(self.example)
        for text, w in [
            ("Incidence matrix A (full or reduced)", self.A),
            ("Tree branches", self.tree),
            ("Z (impedance)", self.Z),
            ("Vs (voltage sources)", self.Vs),
            ("Is (current sources)", self.Is),
        ]:
            lay.addWidget(QLabel(text))
            lay.addWidget(w)
        lay.addWidget(btn)

        self.fig = Figure(figsize=(5, 4))
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasQTAgg(self.fig)
        self.out = QPlainTextEdit()
        self.out.setReadOnly(True)
        self.out.setFont(mono)

        right = QSplitter(Qt.Orientation.Vertical)
        right.addWidget(self.canvas)
        right.addWidget(self.out)
        right.setSizes([450, 300])

        main = QSplitter(Qt.Orientation.Horizontal)
        main.addWidget(left)
        main.addWidget(right)
        main.setSizes([340, 760])
        self.setCentralWidget(main)
        self.resize(1150, 780)

    def load_example(self, i):
        if i == 0:
            return
        ex = EXAMPLES[self.example.currentText()]
        self.A.setPlainText(ex["A"])
        self.tree.setText(ex["tree"])
        self.Z.setText(ex["Z"])
        self.Vs.setText(ex["Vs"])
        self.Is.setText(ex["Is"])
        self.compute()

    def compute(self):
        try:
            A = parse_matrix(self.A.toPlainText())
            b = A.shape[1]
            Z = parse_vector(self.Z.text(), b, "Z")
            Vs = parse_vector(self.Vs.text(), b, "Vs")
            Is = parse_vector(self.Is.text(), b, "Is")
            res = analyze(A, self.tree.text(), Z, Vs, Is)
        except Exception as e:
            self.out.setPlainText(f"Error: {e}")
            return
        draw_graph(self.ax, res)
        self.canvas.draw_idle()
        self.out.setPlainText(report(res))


def main() -> None:
    app = QApplication(sys.argv)
    w = Window()
    w.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
