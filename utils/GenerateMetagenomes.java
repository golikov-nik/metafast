import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Scanner;

public class GenerateMetagenomes {
	public static boolean isValidInputFilename(String filename) {
		return filename.endsWith(".fa") || filename.endsWith(".fasta") || filename.endsWith(".fastq") || filename.endsWith(".fna");
	}

	public static void main(String[] args) throws IOException {
		if (args.length < 5 || args.length > 6) {
			System.out.println(
					"Usage: java GenerateMetagenomes <genomes.fa|folder> <out_prefix> <genome number> <metagenome number> <coeffs> (coefficient)");
			return;
		}

		System.out.println("Loading genomes...");

		String[] genome = new String[Integer.parseInt(args[2])];
		int n = -1;

		BufferedReader in;
		
		long averageGenomeLength = 0;
		String input = args[0];
		if (isValidInputFilename(input)) {
			in = new BufferedReader(new FileReader(input));
			String s = in.readLine();
			while (s != null) {
				if (s.startsWith(">")) {
					n++;
					genome[n] = "";
				} else {
					genome[n] += s.replaceAll("[^AGCT]", "");
					averageGenomeLength += genome[n].length();
				}
				s = in.readLine();
			}
			in.close();
		} else { //one directory for all input files
			File f = new File(input);
			assert f.isDirectory();
			for (File file : f.listFiles()) {
				if (isValidInputFilename(file.getName())) {
					n++;
					in = new BufferedReader(new FileReader(file));
					String s = in.readLine();
					StringBuilder ans = new StringBuilder();
					while (s != null) {
						if (s.startsWith(">")) {
							s = in.readLine();
							continue;
						} else {
							String toAppend = s.replaceAll("[^AGCT]", "");
							ans.append(toAppend);
							averageGenomeLength += toAppend.length();
						}
						s = in.readLine();
					}
					in.close();
					genome[n] = ans.toString();
				}
			}
		}
		n++;
		averageGenomeLength = averageGenomeLength / n;
		System.out.println("Loaded " + n + " genomes");
		for (int i = 0; i < n; i++) {
			System.out.println("Genome " + (i + 1) + ": length = " + genome[i].length());
			if (genome[i].length() < 100) {
				System.out.println("ERROR: Too small genome " + (i + 1));
				return;
			}
		}
		System.out.println();
		
 		int m = Integer.parseInt(args[3]); // metagenomes' number
		double[][] coef = new double[m][n];
		PrintWriter out;
		File coeffsFile = new File(args[4]);
		Scanner sc = new Scanner(coeffsFile);
		for (int i = 0; i < m; i++) { 
			for (int j = 0; j < n; j++) {
				coef[i][j] = Double.parseDouble(sc.next());
			}
		}
		sc.close();
		
		int cov = 50;
		if (args.length > 5)
			cov = Integer.parseInt(args[5]);
		int readLength = 90;
		
		System.out.println("Generating test...");
		List<String>[] reads = new List[m];
		for (int i = 0; i < m; i++) {
			System.out.println("Generating metagenome " + (i + 1) + ":");
			reads[i] = new ArrayList<String>();
			for (int j = 0; j < n; j++) {
				System.out.print("Genome " + (j + 1) + "...  ");
				int readsNumber = (int) Math.round(coef[i][j] * cov * averageGenomeLength / readLength);
				for (int r = 0; r < readsNumber; r++) {
					int pos = (int) ((genome[j].length() - readLength + 1) * Math.random());
					String read = genome[j].substring(pos, pos + readLength);
					reads[i].add(read);
				}
			}

			Collections.shuffle(reads[i]);
			out = new PrintWriter(args[1] + "_" + (i + 1) + ".fa");
			for (int j = 0; j < reads[i].size(); j++) {
				out.println(">read_" + (j + 1) + " length=" + readLength);
				out.println(reads[i].get(j));
			}
			out.close();
		}
		System.out.println("Done!");
	}

}
