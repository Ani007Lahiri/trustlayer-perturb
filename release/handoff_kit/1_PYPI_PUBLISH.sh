# ============================================================
# PyPI PUBLISH — run these on YOUR machine (needs a PyPI account + token)
# ============================================================
# 1. Create a PyPI account at https://pypi.org/account/register/ and an API token
#    (Account settings -> API tokens -> "Add API token", scope: entire account first time).
# 2. From the package root (release/trustlayer/), build and upload:

cd release/trustlayer
python -m pip install --upgrade build twine
python -m build                      # produces dist/*.whl and dist/*.tar.gz
python -m twine upload dist/*         # paste your token as the password (username = __token__)

# 3. Verify from a clean environment:
python -m venv /tmp/tlcheck && /tmp/tlcheck/bin/pip install trustlayer-perturb
/tmp/tlcheck/bin/python -c "from trustlayer import TrustLayer; print('OK', TrustLayer(alpha=0.1))"

# 4. Tag a GitHub release v0.1.0 so the DOI/citation line resolves:
git tag -a v0.1.0 -m "TrustLayer v0.1.0 — calibrated trust for perturbation prediction"
git push origin v0.1.0
# ============================================================
