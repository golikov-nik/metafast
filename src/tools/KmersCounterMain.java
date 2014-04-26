package tools;

import io.IOUtils;
import ru.ifmo.genetics.dna.kmers.ShortKmerIteratorFactory;
import ru.ifmo.genetics.structures.map.ArrayLong2IntHashMap;
import ru.ifmo.genetics.utils.tool.*;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.BoolParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;

import java.io.File;
import java.io.IOException;

/**
 * Created by ulyantsev on 17.04.14.
 *
 */
public class KmersCounterMain extends Tool {

    public static final String NAME = "kmer-counter";
    public static final String DESCRIPTION = "Read files with FASTQ reads and count k-mers with De Bruijn graph";

    static final int LOAD_TASK_SIZE = 1 << 15;

    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size")
            .create());

    public final Parameter<File[]> inputFiles = addParameter(new FileMVParameterBuilder("reads")
            .withShortOpt("i")
            .mandatory()
            .withDescription("list of input files with FASTQ reads from single environment")
            .create());

    public final Parameter<Boolean> plainText = addParameter(new BoolParameterBuilder("plainText")
            .withShortOpt("pt")
            .optional()
            .withDescription("print k-mers in plain text, not in binary (64 + 32) format")
            .withDefaultValue(false)
            .create());

    public final Parameter<Boolean> notPrintNames = addParameter(new BoolParameterBuilder("not-print-names")
            .withShortOpt("npn")
            .optional()
            .withDescription("print only k-mers counts for all possible k-mers (k < 15)")
            .withDefaultValue(false)
            .create());

    @Override
    protected void runImpl() throws ExecutionFailedException {
        int LEN = k.get();

        ArrayLong2IntHashMap hm;
        try {
            hm = IOUtils.loadReads(inputFiles.get(), LEN, LOAD_TASK_SIZE,
                    new ShortKmerIteratorFactory(), availableProcessors.get(), this.logger);
        } catch (IOException e) {
            throw new ExecutionFailedException("Couldn't load kmers", e);
        }

        try {
            String fp = workDir + File.separator + "kmers";
            IOUtils.printKmers(hm, fp, 2);
            info("k-mers printed to " + fp);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Override
    protected void cleanImpl() {
    }

    public static void main(String[] args) {
        new SeqBuilderMain().mainImpl(args);
    }

    public KmersCounterMain() {
        super(NAME, DESCRIPTION);
    }

}