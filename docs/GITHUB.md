# Uploading the project to GitHub

The upload should include the three source trees and Makefiles, runnable job
inputs, Python utilities, documentation, GitHub Actions workflow, and the scripts
and recorded results under `validation/`. The root [.gitignore](../.gitignore) excludes
compiled programs and objects, regenerated `sim1` output, Python caches,
local IDE settings and credential/configuration directories.

[.gitattributes](../.gitattributes) keeps source/documentation line endings
consistent while preserving CSV and example-config bytes. This avoids
changing recorded input hashes during a checkout on another platform.
CSV CRLF endings and the fixtures' final config blank lines are recognized
by the per-path whitespace rules. See the
[Git attributes manual](https://git-scm.com/docs/gitattributes/2.43.0).

## Validate the folder

From the project root:

```bash
python3 validation/run_tests.py --build-mode releaseO3
```

Supply `--eigen-include /path/to/eigen` if needed. The
[validation guide](../validation/README.md) gives the other build modes, sanitizer
checks and optional scientific analyses.

The project is released under the [MIT License](../LICENSE), selected by the
author. Include `LICENSE` in the first commit. Its text follows the
[Open Source Initiative MIT template](https://opensource.org/license/mit).

## Create the repository and first commit

If this folder is not yet a Git repository, initialize it locally:

```bash
git init -b main
```

For older Git versions without `init -b`, use `git init` followed by
`git symbolic-ref HEAD refs/heads/main` before the first commit. If Git is
already initialized, keep its history and use the intended branch instead.

Review and commit the selected files:

```bash
git status --short --untracked-files=all
git add .
git diff --cached --stat
git diff --cached --check
git commit -m "Add Tuckerman disparate-mass MTS implementation and validation"
```

The staged file list should include `initial_conf.csv` and `config.txt`, but
exclude `bin/`, `build/`, generated `OutputTimes/`, local settings and
credential files. `.gitignore` affects untracked files; it does not remove
anything already tracked in an existing repository.

## Connect the GitHub repository

Create an empty repository on your GitHub account. Leave automatic README,
license and `.gitignore` creation unchecked, since these are managed locally.
Then substitute your own account and repository name below:

```bash
git remote add origin git@github.com:YOUR_ACCOUNT/YOUR_REPOSITORY.git
git push -u origin main
```

If using HTTPS authentication, use
`https://github.com/YOUR_ACCOUNT/YOUR_REPOSITORY.git` as the remote URL.
For an existing `origin`, inspect `git remote -v` before changing it.

After the first push, inspect the **Simulation checks** workflow under
GitHub Actions. It runs the regular suite in three build modes
and a sanitizer job. Local validation and hosted CI
results are reported separately.

## Suggested repository description

> C++17 reproduction of Tuckerman et al.'s disparate-mass multiple-time-step
> molecular dynamics algorithm, with serial/OpenMP implementations and a
> uniform velocity-Verlet baseline.

The [algorithm guide](ALGORITHM.md) supplies the reference, equation mapping,
odd-ratio extension and scope of reproduction. Cite that paper for the
numerical method; it is distinct from citing this implementation.
