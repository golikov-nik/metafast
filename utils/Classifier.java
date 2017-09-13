import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;

import org.apache.log4j.Logger;

import io.IOUtils;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import structures.ConnectedComponent;

public class Classifier {
	final static int k = 31;

	public static BigLong2ShortHashMap[] loadHashMaps(File[] input, Logger logg)
			throws ExecutionFailedException, IOException {
		BigLong2ShortHashMap[] ans = new BigLong2ShortHashMap[input.length];
		for (int i = 0; i < ans.length; i++) {
			ans[i] = IOUtils.loadReads(new File[] { input[i] }, k, 0, Runtime.getRuntime().availableProcessors(), logg);
		}
		return ans;
	}

	public static void main(String[] args) throws ExecutionFailedException, IOException {
		if (args.length != 2) {
			System.out.println("Usage: Classifier <genome directory> <components.bin file>");
			return;
		}
		File genomeDir = new File(args[0]);
		assert genomeDir.isDirectory();
		File[] inputGenomeFiles = genomeDir.listFiles();
		File compFile = new File(args[1]);
		Logger logg = Logger.getLogger("");
		BigLong2ShortHashMap[] genHms = loadHashMaps(inputGenomeFiles, logg);
		List<ConnectedComponent> components = ConnectedComponent.loadComponents(compFile);
		double[][] coeffs = new double[components.size()][inputGenomeFiles.length];
		PrintWriter out = new PrintWriter("coeffs.txt");
		for (int i = 0; i < components.size(); i++) {
			System.out.printf("Processing component %d\n", (i + 1));
			ConnectedComponent comp = components.get(i);
			for (int j = 0; j < genHms.length; j++) {
				long commonKmers = 0;
				for (long kmer : comp.kmers) {
					if (genHms[j].contains(kmer))
						commonKmers++;
				}
				coeffs[i][j] = ((double) commonKmers) / comp.kmers.size();
			}
		}
		for (int i = 0; i < coeffs.length; i++) {
			out.println(Arrays.toString(coeffs[i]));
		}
		out.close();
		out = new PrintWriter(new File("classifier.txt"));
		int perfect = 0, ninety = 0, eighty = 0, seventy = 0, sixty = 0, bad = 0;
		for (int i = 0; i < coeffs.length; i++) {
			double maximum = Arrays.stream(coeffs[i]).max().getAsDouble();
			if (maximum == 1.0) {
				perfect++;
			} else if (maximum >= 0.9) {
				ninety++;
			} else if (maximum >= 0.8) {
				eighty++;
			} else if (maximum >= 0.7) {
				seventy++;
			} else if (maximum >= 0.6) {
				sixty++;
			} else {
				bad++;
			}
		}
		out.printf("1.0: %d\n", perfect);
		out.printf("0.9 - 1.0: %d\n", ninety);
		out.printf("0.8 - 0.9: %d\n", eighty);
		out.printf("0.7 - 0.8: %d\n", seventy);
		out.printf("0.6 - 0.7: %d\n", sixty);
		out.printf("0.6-: %d\n", bad);
		out.close();
	}

}
