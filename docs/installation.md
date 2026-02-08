## Installation

Current version: 0.1.1


### Install official release from PyPI

```bash

pip install sage_broccoliqn
## Quick test; if printing version then succeeded
python -c "import sage_broccoliqn; print(sage_broccoliqn.__version__)"

```


### Development install

```bash
git clone https://github.com/denvercal1234GitHub/sage_broccoliqn
cd sage_broccoliqn
pip install -e .
pytest -q
```


## Requirements
anndata >= 0.9
pandas
numpy
scipy
