package run

import (
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"slices"
	"strings"

	"github.com/AlecAivazis/survey/v2"
	"github.com/brentp/vcfgo"
	"github.com/fatih/color"
)

type AnalysisConfig struct {
	VCF string
	// Parameters for the analysis
	Population      string
	WindowSize      int
	StepSize        int
	Rep             int
	Alpha           float64
	MinQTLWidth     int64
	MergeDistance   int64
	OutputFile      string
	HighParentIdx   int
	HighParentName  string
	HighParentDepth int

	LowParentIdx   int
	LowParentName  string
	LowParentDepth int

	OneParentIdx   int
	OneParentName  string
	OneParentDepth int

	HighBulkIdx   int
	HighBulkName  string
	HighBulkDepth int
	HighBulkSize  int

	LowBulkIdx   int
	LowBulkName  string
	LowBulkDepth int
	LowBulkSize  int
	OneBulkIdx   int
	OneBulkName  string
	OneBulkDepth int
	OneBulkSize  int
}

func Run(cfg AnalysisConfig) error {
	type row struct {
		label string
		value string
	}

	bold := color.New(color.Bold).SprintFunc()
	cyan := color.New(color.FgCyan).SprintFunc()

	printTable := func(title string, width int, rows []row, showDashForEmpty bool) {
		fmt.Printf("\n%s\n\n", bold(cyan(title)))
		for _, r := range rows {
			value := r.value
			if showDashForEmpty && value == "" {
				value = "-"
			}
			fmt.Printf("  %-*s %s\n", width, r.label+":", value)
		}
		fmt.Println()
	}

	printTable("  =================== Parameters ===================", 22, []row{
		{"VCF", cfg.VCF},
		{"Population", cfg.Population},
		{"Window size", fmt.Sprintf("%d", cfg.WindowSize)},
		{"Step size", fmt.Sprintf("%d", cfg.StepSize)},
		{"Simulations", fmt.Sprintf("%d", cfg.Rep)},
		{"Alpha", fmt.Sprintf("%v", cfg.Alpha)},
		{"Min QTL length", fmt.Sprintf("%d", cfg.MinQTLWidth)},
		{"Merge distance", fmt.Sprintf("%d", cfg.MergeDistance)},
		{"Output dir/prefix", cfg.OutputFile},
		{"High parent depth", fmt.Sprintf("%d", cfg.HighParentDepth)},
		{"Low parent depth", fmt.Sprintf("%d", cfg.LowParentDepth)},
		{"One parent depth", fmt.Sprintf("%d", cfg.OneParentDepth)},
		{"High bulk depth", fmt.Sprintf("%d", cfg.HighBulkDepth)},
		{"Low bulk depth", fmt.Sprintf("%d", cfg.LowBulkDepth)},
		{"One bulk depth", fmt.Sprintf("%d", cfg.OneBulkDepth)},
		{"High bulk size", fmt.Sprintf("%d", cfg.HighBulkSize)},
		{"Low bulk size", fmt.Sprintf("%d", cfg.LowBulkSize)},
		{"One bulk size", fmt.Sprintf("%d", cfg.OneBulkSize)},
	}, false)

	file, err := os.Open(cfg.VCF)
	if err != nil {
		return err
	}

	var (
		reader  io.Reader = file
		gz      *gzip.Reader
		cleanup = func() {
			if gz != nil {
				_ = gz.Close()
			}
			_ = file.Close()
		}
	)
	defer cleanup()

	if strings.HasSuffix(cfg.VCF, ".gz") {
		gz, err = gzip.NewReader(file)
		if err != nil {
			return err
		}
		reader = gz
	}

	vcfReader, err := vcfgo.NewReader(reader, false)
	if err != nil {
		return err
	}

	sampleNames := vcfReader.Header.SampleNames
	availableSamples := map[int]string{0: "None"}
	for i, name := range sampleNames {
		availableSamples[i+1] = name
	}

	removeSample := func(name string) {
		if name == "" || name == "None" {
			return
		}
		for idx, sample := range availableSamples {
			if sample == name {
				delete(availableSamples, idx)
				return
			}
		}
	}

	chooseSample := func(label string, exclude ...string) (string, error) {
		keys := make([]int, 0, len(availableSamples))
		for idx := range availableSamples {
			keys = append(keys, idx)
		}
		slices.Sort(keys)

		choices := make([]string, 0, len(keys))
		for _, idx := range keys {
			name := availableSamples[idx]
			if name == "None" || !slices.Contains(exclude, name) {
				choices = append(choices, name)
			}
		}
		if len(choices) == 0 {
			return "", fmt.Errorf("no samples available to assign to %s", label)
		}

		var answer string
		prompt := &survey.Select{
			Message: "Select " + label + ":",
			Options: choices,
		}
		if err := survey.AskOne(prompt, &answer, survey.WithValidator(survey.Required)); err != nil {
			return "", err
		}

		removeSample(answer)
		color.Green("  + %s -> %s\n", label, answer)
		return answer, nil
	}

	useConfiguredSample := func(label, name string, exclude ...string) (string, error) {
		if slices.Contains(sampleNames, name) {
			removeSample(name)
			color.Green("  + %s -> %s\n", label, name)
			return name, nil
		}

		color.Yellow("  ! %s %q not found in VCF -- please re-select.\n", label, name)
		return chooseSample(label, exclude...)
	}

	color.Cyan("\n ==================================== Parent selection =========================================== \n\n")

	switch {
	case cfg.HighParentName == "" && cfg.LowParentName == "":
		cfg.HighParentName, err = chooseSample("HIGH PARENT")
		if err != nil {
			return err
		}
		cfg.LowParentName, err = chooseSample("LOW PARENT", cfg.HighParentName)
		if err != nil {
			return err
		}
	case cfg.HighParentName != "" && cfg.LowParentName != "":
		cfg.HighParentName, err = useConfiguredSample("HIGH PARENT", cfg.HighParentName)
		if err != nil {
			return err
		}
		cfg.LowParentName, err = useConfiguredSample("LOW PARENT", cfg.LowParentName, cfg.HighParentName)
		if err != nil {
			return err
		}
	default:
		return fmt.Errorf("incomplete parent configuration: both parents must be set or both left empty")
	}

	color.Cyan("\n =========================================== Bulk selection ============================================= \n\n")

	if cfg.HighBulkName != "" || cfg.LowBulkName != "" || cfg.OneBulkName != "" {
		return fmt.Errorf("bulk samples must be selected interactively")
	}

	cfg.HighBulkName, err = chooseSample("HIGH BULK", cfg.HighParentName, cfg.LowParentName)
	if err != nil {
		return err
	}
	cfg.LowBulkName, err = chooseSample("LOW BULK", cfg.HighParentName, cfg.LowParentName, cfg.HighBulkName)
	if err != nil {
		return err
	}

	printTable("  ================ Selection summary ================", 16, []row{
		{"High parent", cfg.HighParentName},
		{"Low parent", cfg.LowParentName},
		{"One parent", cfg.OneParentName},
		{"High bulk", cfg.HighBulkName},
		{"Low bulk", cfg.LowBulkName},
		{"One bulk", cfg.OneBulkName},
	}, true)

	if cfg.OutputFile != "" {
		dir := filepath.Dir(cfg.OutputFile)
		if dir != "." {
			if err := os.MkdirAll(dir, 0o755); err != nil {
				return fmt.Errorf("failed to create output directory: %w", err)
			}
		}
	}

	highBulkSelected := cfg.HighBulkName != "" && cfg.HighBulkName != "None"
	lowBulkSelected := cfg.LowBulkName != "" && cfg.LowBulkName != "None"
	highParentSelected := cfg.HighParentName != "" && cfg.HighParentName != "None"
	lowParentSelected := cfg.LowParentName != "" && cfg.LowParentName != "None"

	fmt.Println()
	switch {
	case highBulkSelected && lowBulkSelected && !highParentSelected && !lowParentSelected:
		color.Cyan("  -> Bulks only (no parents)\n")
	case highBulkSelected && !lowBulkSelected && highParentSelected && lowParentSelected:
		color.Cyan("  -> One-bulk BSA-seq (HIGH bulk)\n")
	case !highBulkSelected && lowBulkSelected && highParentSelected && lowParentSelected:
		color.Cyan("  -> One-bulk BSA-seq (LOW bulk)\n")
	default:
		color.Cyan("  -> Two-bulk analysis\n")
	}
	fmt.Println()

	return nil
}
