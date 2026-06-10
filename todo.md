Code Review Summary

I've completed a comprehensive review of the GoBSAseq codebase. Here are the findings:

─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

🐛 Bugs Found

1. Missing return statement - cmd/root.go:153

Severity: High

lowParentDepth, err = strconv.Atoi(parentsDepthLst[1])
if err != nil {
    color.Red("parentsDepth is supposed to be in the form a,b (where a and b are integers)")
}
// Missing return here - execution continues despite error

Impact: When parentsDepth parsing fails, the error is printed but the program continues execution with invalid/uninitialized values, leading to undefined behavior.

Fix: Add return after the error message (similar to lines 128, 134, 166, 171).

<<ref_file file="/Users/godwinmafireyi/Developer/GitHub/GoBSAseq/cmd/root.go" lines="150-155" />

─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

⚡ Optimization Opportunities

1. Code Duplication in Sample Selection - run/run.go

Severity: Medium

The sample selection logic (lines 232-474) contains massive code duplication. The same pattern is repeated 4 times for high parent, low parent, high bulk, and low bulk selection
(~240 lines of duplicated logic).

Recommendation: Extract into a reusable function:

func selectSample(sampleName *string, sampleIdx *int, sampleMap map[int]string,
                 sampleNames []string, prompt string) error {
    // Unified sample selection logic
}

<<ref_file file="/Users/godwinmafireyi/Developer/GitHub/GoBSAseq/run/run.go" lines="232-474" />

─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

2. Duplicate Function - utils/utils.go and filter/filter.go

Severity: Low

The getFloat() function is defined identically in both files:

  • <ref_file file="/Users/godwinmafireyi/Developer/GitHub/GoBSAseq/utils/utils.go" lines="12-27" />
  • <ref_file file="/Users/godwinmafireyi/Developer/GitHub/GoBSAseq/filter/filter.go" lines="70-85" />

Recommendation: Keep in utils package only and import from filter.

─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

3. Dead Code - stats/stats.go

Severity: Low

Lines 54-269 contain a large block of commented-out code for the old RawStats implementation. If truly obsolete, it should be removed to improve code maintainability.

<<ref_file file="/Users/godwinmafireyi/Developer/GitHub/GoBSAseq/stats/stats.go" lines="54-269" />

─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

4. Memory Allocation in Smoothing - stats/smoothing.go

Severity: Medium

In gaussianSmooth(), multiple slices are allocated per chromosome per statistic. Could be optimized by:

  • Pre-allocating with known capacity
  • Reusing buffers where possible
  • Using sync.Pool for frequently allocated slices

<<ref_file file="/Users/godwinmafireyi/Developer/GitHub/GoBSAseq/stats/smoothing.go" lines="115-275" />

─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

5. Large Functions Need Refactoring

Severity: Low-Medium

Several functions exceed reasonable complexity:

  • cmd/root.go Run function: 230+ lines
  • filter/filter.go HardFilterVcf: 230+ lines
  • run/run.go Run function: 320+ lines
  • plots/plots.go GeneratePlots: 270+ lines

Recommendation: Break down into smaller, testable functions with single responsibilities.

─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

6. CLI Argument Parsing - cmd/root.go

Severity: Low

Lines 26-174 contain repetitive argument parsing and validation. Could be refactored into helper functions to improve readability and reduce error-prone duplication.

<<ref_file file="/Users/godwinmafireyi/Developer/GitHub/GoBSAseq/cmd/root.go" lines="26-174" />

─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

7. Plot Generation Memory - plots/plots.go

Severity: Low

In GeneratePlots, many temporary slices are allocated per chromosome. Could optimize by:

  • Reusing slice capacity across chromosomes
  • Using more efficient data structures for threshold arrays

<<ref_file file="/Users/godwinmafireyi/Developer/GitHub/GoBSAseq/plots/plots.go" lines="115-234" />

─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

8. Error Handling Inconsistency - stats/genespace.go

Severity: Low

Line 114: Gene space errors are logged but not propagated. This may be intentional (non-fatal), but inconsistent with the rest of the codebase which typically returns errors.

<<ref_file file="/Users/godwinmafireyi/Developer/GitHub/GoBSAseq/stats/genespace.go" lines="112-116" />

─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

📊 Code Quality Observations

Strengths:

  • Good use of parallel processing with goroutines (filter, threshold simulation)
  • Well-structured package organization (stats, filter, plots, run, utils)
  • Comprehensive statistical analysis pipeline
  • Good use of interfaces for gene space intervals
  • Proper error handling in most places

Areas for Improvement:

  • Add unit tests for core functions
  • Consider adding integration tests
  • Add more documentation for complex statistical functions
  • Standardize error message formats
  • Consider using a linter like golangci-lint

─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

🎯 Priority Recommendations

  1. Fix the missing return bug in cmd/root.go:153 immediately (high severity)
  2. Refactor sample selection in run.go to reduce duplication
  3. Remove duplicate getFloat function
  4. Remove dead/commented code in stats.go
  5. Break down large functions for better testability and maintainability
