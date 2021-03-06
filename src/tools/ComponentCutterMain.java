package tools;

import algo.ComponentsBuilder;
import algo.KmerOperations;
import algo.SequencesFinders;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.tool.values.InMemoryValue;
import ru.ifmo.genetics.utils.tool.values.InValue;
import structures.ConnectedComponent;
import structures.Sequence;
import io.IOUtils;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.*;
import ru.ifmo.genetics.utils.NumUtils;

import java.io.*;
import java.util.Deque;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;

import org.apache.commons.lang.mutable.MutableLong;

public class ComponentCutterMain extends Tool {
	public static final String NAME = "component-cutter";
	public static final String DESCRIPTION = "Build graph components from tangled graph";

	public final Parameter<Integer> k = addParameter(
			new IntParameterBuilder("k").mandatory().withShortOpt("k").withDescription("k-mer size").create());

	public final Parameter<Integer> maximalBadFrequency = addParameter(new IntParameterBuilder("maximal-bad-frequency")
			.optional().withShortOpt("b").withDescription("maximal frequency for a k-mer to be assumed erroneous")
			.withDefaultValue(1).create());

	public final Parameter<Integer> alg = addParameter(new IntParameterBuilder("alg").important().withShortOpt("alg")
			.withDescription("number of an algorithm").withDefaultValue(0).create());

	public final Parameter<Integer> minLen = addParameter(new IntParameterBuilder("min-seq-len").important()
			.withShortOpt("l").withDescription("minimum sequence length to be added").withDefaultValue(100).create());

	public final Parameter<Integer> minComponentSize = addParameter(new IntParameterBuilder("min-component-size")
			.important().withShortOpt("b1").withDescription("minimum component size in component-cutter (in k-mers)")
			.withDefaultValue(1000).withDescriptionShort("Minimal component size")
			.withDescriptionRu("Минимальный размер компоненты для использования")
			.withDescriptionRuShort("Минимальный размер компоненты").create());

	public final Parameter<Integer> maxComponentSize = addParameter(new IntParameterBuilder("max-component-size")
			.important().withShortOpt("b2").withDescription("maximum component size in component-cutter (in k-mers)")
			.withDefaultValue(10000).withDescriptionShort("Maximal component size")
			.withDescriptionRu("Максимальный размер компоненты для использования")
			.withDescriptionRuShort("Максимальный размер компоненты").create());

	public final Parameter<File[]> sequencesFiles = addParameter(new FileMVParameterBuilder("sequences")
			.withShortOpt("i").mandatory().withDescription("list of input files").create());

	public final Parameter<File> componentsFile = addParameter(
			new FileParameterBuilder("components-file").withDescription("file to write found components to")
					.withDefaultValue(workDir.append("components.bin")).create());

	public File[] outputDescFiles = null;

	// output values
	public final InValue<File> componentsFileOut = addOutput("components-file", componentsFile, File.class);
	private final InMemoryValue<File> componentsStatPr = new InMemoryValue<File>();
	public final InValue<File> componentsStatOut = addOutput("components-stat", componentsStatPr, File.class);

	@Override
	protected void runImpl() throws ExecutionFailedException, IOException {
		Timer t = new Timer();
		debug("Loading sequences from files...");
		File[] seqFiles = sequencesFiles.get();
		BigLong2ShortHashMap hm = IOUtils.loadReads(seqFiles, k.get(), minLen.get(), availableProcessors.get(), logger);
		Deque<Sequence> sequences = null;
		if (alg.get() == 5)
			try {
				sequences = SequencesFinders.thresholdStrategy(hm, availableProcessors.get(), maximalBadFrequency.get(),
						maxComponentSize.get() + k.get() - 1, k.get());
			} catch (InterruptedException e1) {
				e1.printStackTrace();
			}
		debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);
		if (hm.size() == 0) {
			throw new ExecutionFailedException(
					"No sequences were found in input files! The following steps will be useless");
		}
		int solo = 0;
		int root = 0;
		int link = 0;
		int multi = 0;
		for (MutableLong kmer: hm) {
			if (hm.get(kmer.toLong()) <= 0) continue;
			int neighbors = ComponentsBuilder.neighboursInGraph(hm, kmer.toLong(), k.get()).length;
			if (neighbors == 0) solo++;
			else if (neighbors == 1) root++;
			else if (neighbors == 2) link++;
			else multi++;
		}
		info(String.format("Found %d solo kmers, %d roots, %d links, %d multi kmers\n", solo, root, link, multi));
		info("Searching for components...");
		List<ConnectedComponent> components;
		try {
			String statFP = workDir + File.separator + "components-stat-" + minComponentSize.get() + "-"
					+ maxComponentSize.get() + ".txt";
			components = ComponentsBuilder.splitStrategy(alg.get(), hm, sequences, k.get(), minComponentSize.get(),
					maxComponentSize.get(), statFP, logger, availableProcessors.get());

			componentsStatPr.set(new File(statFP));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			return;
		}
		info("Total " + NumUtils.groupDigits(components.size()) + " components were found");
		if (components.size() == 0) {
			warn("No components were extracted! Perhaps you should decrease --min-component-size value");
		}
		
		try {
			ConnectedComponent.saveComponents(components, componentsFile.get().getAbsolutePath());
			// ConnectedComponent.printComponents(logger, components, k.get(),
			// workDir.get().gettAbsolutePath());
			long size = components.stream().map(e -> e.size).reduce(Long::sum).get();
			info(String.format("Average component size %d\n", Math.round((double)size/components.size())));
			info("Components saved to " + componentsFile.get());
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		debug("Components-cutter has finished! Time = " + t);
	}

	@Override
	protected void cleanImpl() {
	}

	@Override
	protected void postprocessing() {
		IOUtils.tryToAppendDescription(outputDescFiles, componentsStatOut.get(),
				"File with components' statistics (in text format)");
		IOUtils.tryToAppendDescription(outputDescFiles, componentsFileOut.get(),
				"File with extracted components (in binary format)");
	}

	public static void main(String[] args) {
		new ComponentCutterMain().mainImpl(args);
	}

	public ComponentCutterMain() {
		super(NAME, DESCRIPTION);
	}
}
