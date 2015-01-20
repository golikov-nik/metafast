package tools;

import ru.ifmo.genetics.utils.FileUtils;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.BoolParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;

import java.io.*;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

public class DistanceMatrixCalculatorMain extends Tool {
    public static final String NAME = "dist-matrix-calculator";
    public static final String DESCRIPTION = "Calculates distance matrix using features values";

    public static final String SEPARATOR = "\t";


    public final Parameter<File[]> featuresFiles = addParameter(new FileMVParameterBuilder("features")
            .mandatory()
            .withDescription("features values files (for same components)")
            .create());

    public final Parameter<Boolean> withoutNames = addParameter(new BoolParameterBuilder("without-names")
            .withShortOpt("wn")
            .optional()
            .withDescription("do not print matrix row and column names (as given file names)")
            .create());

    public final Parameter<File> matrixFile = addParameter(new FileParameterBuilder("matrix-file")
            .optional()
            .withDefaultValue(workDir.append("dist_matrix_$DT.txt"))
            .withDefaultComment("<workDir>/dist_matrix_<date>_<time>.txt")
            .withDescription("resulting distance matrix file")
            .create());


    @Override
    protected void runImpl() throws ExecutionFailedException {
        List<List<Long>> features = new ArrayList<List<Long>>();

        for (File featuresFile : featuresFiles.get()) {
            try {
                features.add(readVector(featuresFile));
            } catch (IOException e) {
                throw new ExecutionFailedException("Failed to read features from " + featuresFile);
            }
        }

        int cnt = features.size();
        double[][] distMatrix = new double[cnt][cnt];

        for (int i = 0; i < cnt; i++) {
            for (int j = i + 1; j < cnt; j++) {
                distMatrix[i][j] = brayCurtisDistance(features.get(i), features.get(j));
                distMatrix[j][i] = distMatrix[i][j];
            }
        }

        String timestamp = new SimpleDateFormat("yyyy.MM.dd_HH.mm.ss").format(new Date());
        String matrixPath = matrixFile.get().getPath().replace("$DT", timestamp);
        String[] names = new String[featuresFiles.get().length];
        for (int i = 0; i < names.length; i++) {
            names[i] = featuresFiles.get()[i].getName();
        }

        try {
            printMatrix(distMatrix, matrixPath, names, null);
            info("Distance matrix printed to " + matrixPath);
        } catch (FileNotFoundException e) {
            throw new ExecutionFailedException("Failed to print matrix to " + matrixPath);
        }
        matrixFile.set(new File(matrixPath));
    }

    public static void printMatrix(double[][] matrix, String fp, String[] names, int[] perm) throws FileNotFoundException {
        File f = new File(fp);
        FileUtils.makeSubDirsOnly(f);
        PrintWriter out = new PrintWriter(f);

        if (names != null) {
            out.print("#");
            for (int i = 0; i < names.length; i++) {
                out.print(SEPARATOR);
                out.print(names[(perm == null) ? i : perm[i]]);
            }
            out.println();
        }

        for (int i = 0; i < matrix.length; i++) {
            if (names != null) {
                out.print(names[(perm == null) ? i : perm[i]] + SEPARATOR);
            }
            for (int j = 0; j < matrix[i].length; j++) {
                if (j > 0) {
                    out.print(SEPARATOR);
                }
                if (perm == null) {
                    out.print(matrix[i][j]);
                } else {
                    out.print(matrix[perm[i]][perm[j]]);
                }
            }
            out.println();
        }

        out.close();
    }

    private List<Long> readVector(File file) throws IOException {
        List<Long> ans = new ArrayList<Long>();

        BufferedReader reader = new BufferedReader(new FileReader(file));
        String line;
        while ((line = reader.readLine()) != null) {
            if (!line.isEmpty()) {
                ans.add(Long.parseLong(line));
            }
        }
        reader.close();

        return ans;
    }

    private double brayCurtisDistance(List<Long> vector1, List<Long> vector2) {
        assert vector1.size() == vector2.size();

        double sumdiff = 0, sum = 0;

        for (int pos = 0; pos < vector1.size(); pos++) {
            sumdiff += Math.abs(vector1.get(pos) - vector2.get(pos));
            sum += Math.abs(vector1.get(pos)) + Math.abs(vector2.get(pos));
        }

        assert sum > 0;
        return sumdiff / sum;
    }


    @Override
    protected void cleanImpl() {

    }

    public static void main(String[] args) {
        new SeqBuilderMain().mainImpl(args);
    }

    public DistanceMatrixCalculatorMain() {
        super(NAME, DESCRIPTION);
    }
}
