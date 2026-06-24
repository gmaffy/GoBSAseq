# Code review

## Overall assessment
The project is functional and the current test suite passes. The recent QTL detection work in [stats/detect.go](stats/detect.go) is a clear improvement, but the codebase still has a few structural issues that will make future changes more error-prone.

Verified with:
- `go test ./...`
- `go vet ./...`

## Findings

### 1. CLI wiring is overly procedural and hard to maintain
The command entrypoint in [cmd/root.go](cmd/root.go) parses and validates a large number of flags inline, then passes them through many intermediate variables. This makes the flow difficult to read and easy to break when new options are added.

Suggested direction:
- Introduce a dedicated config/option struct for CLI parsing.
- Move validation into small helper functions.
- Keep the command layer focused on orchestration rather than argument massaging.

### 2. Pipeline orchestration and user interaction are mixed together
The workflow in [run/run.go](run/run.go) combines file handling, sample prompting, VCF/BAM selection, filtering, stats, and execution logic in one place. That makes the code hard to test and harder to reuse in non-interactive contexts.

Suggested direction:
- Split the workflow into smaller stages such as input discovery, sample selection, analysis execution, and output writing.
- Replace direct terminal prompts with a more explicit interface where possible.

### 3. Several functions encode the same domain rules in multiple places
The BSA analysis type logic is effectively duplicated across [run/run.go](run/run.go), [stats/smoothing.go](stats/smoothing.go), and related stats code. This creates drift risk: one helper may accept a mode that another helper does not understand.

Suggested direction:
- Centralize the BSA type / bulk combination rules in one shared helper.
- Reuse that helper everywhere instead of re-implementing the same branching logic.

### 4. Some correctness checks use panics instead of returned errors
The invariant check in [stats/smoothing.go](stats/smoothing.go) uses `panic` when the CompositeZ statistic count does not match expectations. That is reasonable as a guardrail, but it makes the library less composable and can crash the whole process instead of returning a controlled error.

Suggested direction:
- Prefer returning an error from public entrypoints when invariants are violated.
- Reserve panics for truly unrecoverable internal failures.

### 5. Legacy and backup files increase maintenance risk
The repository contains several duplicate or historical files such as [stats/detect.go.1](stats/detect.go.1), [stats/detect.go.2](stats/detect.go.2), [stats/detect.go.working](stats/detect.go.working), and [run/run.go.1](run/run.go.1). These files can confuse contributors and make it unclear which implementation is the current one.

Suggested direction:
- Remove obsolete copies or move them into an archive folder.
- Keep only the active implementation in the main source tree.

### 6. Test coverage is still thin for the end-to-end workflow
The repository has tests in [stats/detect_test.go](stats/detect_test.go), which is positive, but most of the pipeline layers are not covered. That means regressions in CLI parsing, sample selection, and workflow glue are more likely to slip through.

Suggested direction:
- Add tests for the command/config parsing path.
- Add tests for the high-level workflow steps in [run/run.go](run/run.go) where practical.
- Keep regression tests next to the logic that changed.

## Priority recommendations
1. Reduce duplication in the pipeline and CLI layers.
2. Consolidate BSA type handling in one place.
3. Remove or archive legacy source files.
4. Add focused tests for the non-stats workflow components.
