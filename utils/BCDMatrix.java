import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Locale;
import java.util.Scanner;

public class BCDMatrix {
	
	static ArrayList<Double>[] readDoublesWithExceptions(File f, int m) throws FileNotFoundException {
		Scanner in = new Scanner(f);
		ArrayList<Double>[] ans = new ArrayList[m];
		for (int i = 0; i < ans.length; i++) {
			ans[i] = new ArrayList<>();
		}
		Scanner sc;
		for (int i = 0; i < m; i++) {
			String s = in.nextLine();
			sc = new Scanner(s);
			while (sc.hasNext()) {
				try {
					ans[i].add(Double.parseDouble(sc.next()));
				} catch (NumberFormatException e) {
					continue;
				}
			}
		}
		in.close();
		return ans;
	}
	
	public static void main(String[] args) throws FileNotFoundException {
		if (args.length != 3) {
			System.out.println("Usage: BCDMatrix <matrix> <vector_number> <out_prefix>");
			return;
		}
		int m = Integer.parseInt(args[1]);
		ArrayList<Double>[] metagenomes = readDoublesWithExceptions(new File(args[0]), m);
		PrintWriter out = new PrintWriter(new File(args[2] + ".txt"));
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				out.printf(Locale.US, "%.4f", Statistics.brayCurtisDistance(metagenomes[i], metagenomes[j]));
				out.print(" ");
			}
			out.println();
		}
		out.close();
	}

}
