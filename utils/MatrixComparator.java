import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class MatrixComparator {
	static List<Double> readDoublesWithExceptions(File f) throws FileNotFoundException {
		Scanner in = new Scanner(f);
		List<Double> ans = new ArrayList<>();
		while (in.hasNext()) {
			try {
				ans.add(Double.parseDouble(in.next()));
			} catch (NumberFormatException e) {
				continue;
			}
		}
		in.close();
		return ans;
	}
	
	public static void main(String[] args) throws FileNotFoundException {
		if (args.length != 2) {
			System.out.println("Usage: java MatrixComparator <matrix1> <matrix2>");
			return;
		}
		File matrix1File = new File(args[0]);
		File matrix2File = new File(args[1]);
		List<Double> matrix1 = readDoublesWithExceptions(matrix1File);
		List<Double> matrix2 = readDoublesWithExceptions(matrix2File);
		if (matrix1.size() != matrix2.size()) {
			System.out.println("Error: Matrices must be the same size");
			return;
		}
		System.out.println("Calculating statistics...\n");
		double bcd = Statistics.brayCurtisDistance(matrix1, matrix2);
		double spear = Statistics.spearmansCorr(matrix1, matrix2);
		double pears = Statistics.pearsonsCorr(matrix1, matrix2);
		System.out.println("Results:");
		System.out.println("Bray-Curtis distance: " + bcd);
		System.out.println("Spearman's correlation: " + spear);
		System.out.println("Pearson's correlation: " + pears);
	}

}
