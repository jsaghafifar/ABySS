package abyss;

import lphy.base.evolution.eigensystem.ComplexColtEigenSystem;
import lphy.base.evolution.eigensystem.EigenDecompositionExt;
import lphy.base.evolution.eigensystem.EigenSystem;
import lphy.core.logger.LoggerUtils;
import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;

/**
 * @author Jasmine Saghafifar
 */

public class ComputeEquilibrium extends DeterministicFunction<Double[]> {
    protected static final String qParamName = "Q";
    private static final double DEFAULT_BRANCH_LENGTH = 10000;

    public ComputeEquilibrium(@ParameterInfo(name = qParamName,
            description = "the Q matrix to compute equilibrium frequencies from.") Value<Double[][]> Q) {
        setParam(qParamName, Q);
    }

    @GeneratorInfo(name = "computeEquilibrium", description = "The sum of the elements of the given array")
    public Value<Double[]> apply() {
        Double[][] Q = (Double[][])getParams().get(qParamName).value();
        return new Value<>(null, computeEquilibrium(Q), this);
    }

    private Double[] computeEquilibrium(Double[][] Q) {
        int numStates = Q.length;

        // Q assumed to be nonreversible, so complex eigen decomposition
        EigenSystem complexEigenSystem = new ComplexColtEigenSystem(numStates);
        EigenDecompositionExt complexDecomposition = complexEigenSystem.decomposeMatrix(Q);
        double[] Eval = complexDecomposition.getEigenValues();
        double[] evec = complexDecomposition.getEigenVectors();
        double[] ievc = complexDecomposition.getInverseEigenVectors();

        double[] EvalImag = new double[numStates];
        Double[] freqs = new Double[numStates];
        double[][] iexp = new double[numStates][numStates];
        double[][] transProbs = new double[numStates][numStates];
        System.arraycopy(Eval, numStates, EvalImag, 0, numStates);

        double temp;
        int k;
        double branchLength = DEFAULT_BRANCH_LENGTH;
        boolean equilibrium = false;

        // get transition probabilities with bigger branch lengths until equilibrium is reached
        while (!equilibrium) {
            for (int i = 0; i < numStates; i++) {
                if (EvalImag[i] == 0) {
                    // 1x1 block
                    temp = Math.exp(branchLength * Eval[i]);
                    for (int j = 0; j < numStates; j++) {
                        iexp[i][j] = ievc[i * numStates + j] * temp;
                    }
                } else {
                    // 2x2 conjugate block
                    // If A is 2x2 with complex conjugate pair eigenvalues a +/- bi, then
                    // exp(At) = exp(at)*( cos(bt)I + \frac{sin(bt)}{b}(A - aI)).
                    int i2 = i + 1;
                    double b = EvalImag[i];
                    double expat = Math.exp(branchLength * Eval[i]);
                    double expatcosbt = expat * Math.cos(branchLength * b);
                    double expatsinbt = expat * Math.sin(branchLength * b);

                    for (int j = 0; j < numStates; j++) {
                        iexp[i][j] = expatcosbt * ievc[i * numStates + j] +
                                expatsinbt * ievc[i2 * numStates + j];
                        iexp[i2][j] = expatcosbt * ievc[i2 * numStates + j] -
                                expatsinbt * ievc[i * numStates + j];
                    }
                    i++; // processed two conjugate rows
                }
            }
            for (int i = 0; i < numStates; i++) {
                for (int j = 0; j < numStates; j++) {
                    temp = 0.0;
                    for (k = 0; k < numStates; k++) {
                        temp += evec[i * numStates + k] * iexp[k][j];
                    }
                    transProbs[i][j] = Math.abs(temp);
                }
            }
            boolean reached = true;
            for (int i = 0; i < numStates; i++) {
                freqs[i] = transProbs[0][i];
                for (int j = 1; j < numStates; j++) {
                    if (Math.abs(transProbs[0][i] - transProbs[j][i]) > 1e-5) {
                        reached = false;
                        break;
                    }
                }
                if (!reached) break;
            }
            if (reached) equilibrium = true;
            branchLength *= 10;
            if (branchLength > 1e10) LoggerUtils.log.severe("Could not compute equilibrium.");
        }
        return freqs;
    }

}
