package tools;

import ru.ifmo.genetics.Runner;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.BoolParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;

public class DistanceMatrixBuilderMain extends Tool {
    public static final String NAME = "matrix-builder";
    public static final String DESCRIPTION = "Builds the distance matrix for input sequences";


    static {
        forceParameter.set(true);
        launchOptions.remove(forceParameter);
    }

    // input parameters
    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .important()
            .withShortOpt("k")
            .withDefaultValue(31)
            .withDescription("k-mer size (in nucleotides, maximum 31 due to realization details)")
            .create());

    public final Parameter<File[]> inputFiles = addParameter(new FileMVParameterBuilder("reads")
            .withShortOpt("i")
            .mandatory()
            .withDescription("list of reads files from single environment. FASTQ, BINQ, FASTA")
            .create());

    public final Parameter<Integer> maximalBadFrequency = addParameter(new IntParameterBuilder("maximal-bad-frequency")
            .important()
            .withShortOpt("b")
            .withDescription("maximal frequency for a k-mer to be assumed erroneous")
            .withDefaultValue(1)
            .create());

    public final Parameter<Integer> minLen = addParameter(new IntParameterBuilder("min-seq-len")
            .important()
            .withShortOpt("l")
            .withDescription("minimum sequence length to be added to a component (in nucleotides)")
            .withDefaultValue(100)
            .create());

    public final Parameter<Boolean> useReadsForCalculatingFeatures = addParameter(new BoolParameterBuilder("use-reads-for-calculating-features")
            .withDescription("use reads instead of kmers for calculating features (characteristic vectors)")
            .create());

    public final Parameter<File> matrixFile = addParameter(new FileParameterBuilder("matrix-file")
            .withDefaultValue(workDir.append("matrices").append("dist_matrix_$DT.txt"))
            .withDefaultComment("<workDir>/matrices/dist_matrix_<date>_<time>.txt")
            .withDescription("resulting distance matrix file")
            .create());

    public final Parameter<File> heatmapFile = addParameter(new FileParameterBuilder("heatmap-file")
            .withDefaultValue(workDir.append("matrices").append("dist_matrix_$DT_heatmap.png"))
            .withDefaultComment("<workDir>/matrices/dist_matrix_<date>_<time>_heatmap.png")
            .withDescription("resulting heatmap file")
            .create());

    public final File[] outputDescFiles = new File[]{
                    new File("output_description.txt"),
                    workDir.append("output_description.txt").get()};



    // adding sub tools
    public KmersCounterForManyFilesMain kmersCounter = new KmersCounterForManyFilesMain();
    {
        setFix(kmersCounter.k, k);
        setFix(kmersCounter.inputFiles, inputFiles);
        setFix(kmersCounter.maximalBadFrequency, maximalBadFrequency);
        setFixDefault(kmersCounter.outputDir);
        kmersCounter.outputDescFiles = outputDescFiles;
        addSubTool(kmersCounter);
    }

    public SeqBuilderForManyFilesMain seqBuilder = new SeqBuilderForManyFilesMain();
    {
        setFix(seqBuilder.k, k);
        setFix(seqBuilder.inputFiles, kmersCounter.resultingKmerFiles);
        setFix(seqBuilder.maximalBadFrequency, maximalBadFrequency);
        setFix(seqBuilder.sequenceLen, minLen);
        setFixDefault(seqBuilder.outputDir);
        seqBuilder.outputDescFiles = outputDescFiles;
        addSubTool(seqBuilder);
    }

    public ComponentCutterMain compCutter = new ComponentCutterMain();
    {
        setFix(compCutter.k, k);
        setFix(compCutter.sequencesFiles, seqBuilder.outputFilesOut);
        setFix(compCutter.minLen, minLen);
        setFixDefault(compCutter.componentsFile);
        compCutter.outputDescFiles = outputDescFiles;
        addSubTool(compCutter);
    }

    public FeaturesCalculatorMain featuresCalculator = new FeaturesCalculatorMain();
    {
        setFix(featuresCalculator.k, k);
        setFix(featuresCalculator.componentsFile, compCutter.componentsFileOut);
        setFixDefault(featuresCalculator.readsFiles);
        setFix(featuresCalculator.kmersFiles, kmersCounter.resultingKmerFiles);
        setFixDefault(featuresCalculator.threshold);
        featuresCalculator.outputDescFiles = outputDescFiles;
        addSubTool(featuresCalculator);
    }

    public DistanceMatrixCalculatorMain distMatrixCalculator = new DistanceMatrixCalculatorMain();
    {
        setFix(distMatrixCalculator.featuresFiles, featuresCalculator.featuresFilesOut);
        setFix(distMatrixCalculator.matrixFile,
                workDir.append("matrices").append("dist_matrix_$DT_original_order.txt"));
        distMatrixCalculator.outputDescFiles = outputDescFiles;
        addSubTool(distMatrixCalculator);
    }

    public HeatMapMakerMain heatMapMaker = new HeatMapMakerMain();
    {
        setFix(heatMapMaker.matrixFile, distMatrixCalculator.matrixFile);
        setFix(heatMapMaker.newMatrixFile, matrixFile);
        setFix(heatMapMaker.heatmapFile, heatmapFile);
        heatMapMaker.outputDescFiles = outputDescFiles;
        addSubTool(heatMapMaker);
    }


    @Override
    protected void runImpl() throws ExecutionFailedException {
        // preparing
        info("Found " + inputFiles.get().length + " libraries to process");
        if (useReadsForCalculatingFeatures.get()) {
            setFix(featuresCalculator.kmersFiles, new File[]{});
            setFix(featuresCalculator.readsFiles, inputFiles);
        }
        createOutputDescFiles();

        // running steps
        addStep(kmersCounter);
        addStep(seqBuilder);
        addStep(compCutter);
        addStep(featuresCalculator);
        addStep(distMatrixCalculator);
        addStep(heatMapMaker);
    }

    private void createOutputDescFiles() {
        for (File f : outputDescFiles) {
            try {
                PrintWriter out = new PrintWriter(f);
                out.println("# Output files' description " +
                        "for run started at " +
                        new SimpleDateFormat("dd-MMM-yyyy (EEE) HH:mm:ss").format(startDate));
                out.println();
                for (File ff : getLogFiles()) {
                    out.println(ff);
                }
                out.println("   Identical files with run log");

                out.println();
                for (File ff : outputDescFiles) {
                    out.println(ff);
                }
                out.println("   Identical files with output files' description");
                out.close();
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
        }
    }


    @Override
    protected void cleanImpl() {
    }


    @Override
    public void mainImpl(String[] args) {
        if (Runner.containsOption(args, Runner.getOptKeys(continueParameter)) ||
                Runner.containsOption(args, Runner.getOptKeys(Tool.startParameter))) {
            forceParameter.set(false);
        }
        super.mainImpl(args);
    }

    public static void main(String[] args) {
        new DistanceMatrixBuilderMain().mainImpl(args);
    }

    public DistanceMatrixBuilderMain() {
        super(NAME, DESCRIPTION);
    }
}
