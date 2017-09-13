import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;

public class Statistics {
	private static double[] doubleListAsArr(List<Double> list) {
		return list.stream().mapToDouble(Double::doubleValue).toArray();
	}
	
	public static double spearmansCorr(List<Double> vector1, List<Double> vector2) {
		SpearmansCorrelation spear = new SpearmansCorrelation();
		double[] arr1 = doubleListAsArr(vector1);
		double[] arr2 = doubleListAsArr(vector2);
		return spear.correlation(arr1, arr2);
	}
	
	public static double pearsonsCorr(List<Double> vector1, List<Double> vector2) {
		PearsonsCorrelation spear = new PearsonsCorrelation();
		double[] arr1 = doubleListAsArr(vector1);
		double[] arr2 = doubleListAsArr(vector2);
		return spear.correlation(arr1, arr2);
	}
	
	public static double brayCurtisDistance(List<Double> vector1, List<Double> vector2) {
		assert vector1.size() == vector2.size();

		double sumdiff = 0.0, sum = 0.0;

		for (int pos = 0; pos < vector1.size(); pos++) {
			sumdiff += Math.abs(vector1.get(pos) - vector2.get(pos));
			sum += Math.abs(vector1.get(pos)) + Math.abs(vector2.get(pos));
		}

		assert sum > 0;
		return sumdiff / (double) sum;
	}
	
}
