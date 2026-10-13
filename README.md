# CADaid

CAD network topology tool: builds the reduced incidence matrix `A`, tie-set matrix `B`,
cut-set matrix `C`, and solves branch currents `J` and voltages `V`.

## Requirements

- Python **3.12+**
- One of: [`uv`](https://docs.astral.sh/uv/) (recommended) or `pip`

## Installation

### Option A — `uv` (recommended)

Install `uv`:
-check website [uv website](https://docs.astral.sh/uv/getting-started/installation/)

```bash
git clone <repo-url>
cd cadaid
uv sync                
uv run cadaid          # run
```

### Option B — `pip`

```sh
git clone <repo-url>
cd cadaid
py -3.12 -m venv .venv

.venv\Scripts\activate # Windows
source .venv/bin/activate # Mac / Linux

python -m pip install -U pip
pip install .
cadaid # run
```

## Usage notes

- Input A: full (N rows) or reduced (N-1 rows, last node is reference).
- Convention: +1 current leaves node, -1 enters node.
- Tree: leave empty for auto, or list branch letters (a b c ...).
- Z, Vs, Is: optional, one value per branch (letters order). Solves branch J and V.
- Branch voltage: V = Z (J + Is) - Vs
