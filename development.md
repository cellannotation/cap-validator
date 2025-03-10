### Setting Up the Environment

```bash
git clone https://github.com/cellannotation/cap-validator.git
cd cap-validator
python -m pip install uv
uv venv --python 3.13  # or any python version 3.9+
source .venv/bin/activate
uv pip install -r pyproject.toml
```

### Running Tests
```bash
pytest test/
```

### Linting & Formatting
```bash
black .
flake8
```
