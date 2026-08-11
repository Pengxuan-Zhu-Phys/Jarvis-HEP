# Releasing Jarvis-HEP

## Package identity

- PyPI distribution: `Jarvis-HEP`
- Python import package: `jarvishep2`
- Command-line entry point: `Jarvis`
- V1 releases stop at the `1.x` series; V2 starts at `2.0.0` on the same
  `Jarvis-HEP` PyPI project.

PyPI normalizes distribution names, but it does not infer them from the import
directory. Keep `[project].name = "Jarvis-HEP"` in `pyproject.toml` even though
the source directory is named `jarvishep2`.

## One-time Trusted Publishing setup

The release workflow uses PyPI Trusted Publishing (GitHub OIDC), so no PyPI API
token is stored in GitHub.

1. In GitHub, open **Settings → Environments**, create an environment named
   `pypi`, and require a maintainer's approval for deployment if the repository
   plan supports it.
2. In the PyPI `Jarvis-HEP` project, open **Manage → Publishing**, add a GitHub
   Actions trusted publisher, and enter:

   | Field | Value |
   | --- | --- |
   | Owner | `Pengxuan-Zhu-Phys` |
   | Repository | `Jarvis-HEP` |
   | Workflow filename | `publish-pypi.yml` |
   | Environment | `pypi` |

The workflow file must be merged into the repository's default branch before
the first release is published.

## Release procedure

1. Update the version in both `pyproject.toml` and
   `jarvishep2/__init__.py`. PyPI versions are immutable, so choose a version
   that does not already exist under `Jarvis-HEP`.
2. Run the local release checks:

   ```bash
   python3 -m pip install --upgrade build twine
   python3 -m pip install '.[distributed,dev]'
   python3 -m pytest -q
   python3 -m build
   python3 -m twine check dist/*
   ```

   The default test command intentionally skips the long-running
   `tests/test_adaptive_bridson.py`, `tests/test_ensemble_samplers.py`,
   `tests/test_distributed_acceptance.py`, `tests/test_distributed_resume.py`,
   `tests/test_mcmc_sampler.py`, `tests/test_worker_pool.py`,
   `tests/test_variable_distributions.py`, and `tests/test_worker_failure.py`;
   run the relevant files separately before sampler- or distributed-runtime-
   specific releases.

3. Merge the release commit into `master`.
4. Create a GitHub Release whose tag is exactly `v<version>`, for example:

   ```bash
   gh release create v2.0.0 --target master --generate-notes --title "Jarvis-HEP 2.0.0"
   ```

Publishing the GitHub Release triggers `.github/workflows/publish-pypi.yml`.
The workflow checks that the tag, both source version declarations, V2 major
version, distribution name, and wheel contents agree. It then runs the default
test suite, builds one wheel and one source distribution, and publishes those
exact artifacts to PyPI using a short-lived OIDC credential.

Do not upload with a local token as part of the normal release process. Do not
enable `skip-existing`: a duplicate version should fail visibly instead of
hiding a release mistake.

## Migration from the temporary `jarvishep2` distribution

The temporary `jarvishep2` distribution and `Jarvis-HEP` V2 both contain the
same `jarvishep2` Python package, so they should not coexist in one environment:

```bash
python3 -m pip uninstall jarvishep2
python3 -m pip install --upgrade 'Jarvis-HEP[distributed]'
```

Existing V1 users can upgrade directly from `Jarvis-HEP` 1.x to 2.x, but V2 is
a major-version migration: the old `jarvishep` import package is replaced by
`jarvishep2`.
