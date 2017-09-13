import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Scanner;

import org.apache.log4j.Logger;

import io.IOUtils;
import ru.ifmo.genetics.dna.Dna;
import ru.ifmo.genetics.dna.kmers.ShortKmer;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;

public class Test {
	final static int k = 31;
	
	public static String determineColor() {
		int value1 = (int) (256 * Math.random());
		int value2 = (int) (256 * Math.random());
		int value3 = (int) (256 * Math.random());
        return String.format("#%02X%02X%02X", value1, value2, value3);
	}
	
	public static void main(String[] args) throws ExecutionFailedException, IOException {
		Scanner in = new Scanner(new File("C:\\Users\\golik\\workspace\\metafast\\src\\comp_1.gfa"));
		PrintWriter out = new PrintWriter(new File("C:\\Users\\golik\\workspace\\metafast\\src\\comp_1_modified.gfa"));
		out.println(in.nextLine());
		out.print(in.next() + "\t");
		out.print(in.next() + "\t");
		String dna = in.next();
		out.print(dna + "\n");
		in.close();

		File genomesDir = new File("C:\\Users\\golik\\workspace\\metafast\\src\\genomes");
		File[] inputFiles = genomesDir.listFiles();
		Logger logg = Logger.getLogger("1");
		BigLong2ShortHashMap[] hms = new BigLong2ShortHashMap[inputFiles.length];
		for (int i = 0; i < hms.length; i++) {
			hms[i] = IOUtils.loadReads(new File[] { inputFiles[i] }, k, 0, Runtime.getRuntime().availableProcessors(),
					logg);
		}
		boolean[] genomes = new boolean[inputFiles.length];
		for (int i = 0; i < genomes.length; i++) {
			boolean mark = true;
			for (ShortKmer kmer : ShortKmer.kmersOf(new Dna(dna), k)) {
				if (!hms[i].contains(kmer.toLong())) {
					mark = false;
					break;
				}
			}
			genomes[i] = mark;
		}
		
		String color = determineColor();
        out.println("\tCL:z:" + color + "\tC2:z:" + color);
        out.close();
	}

}
