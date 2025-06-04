package abyss;

import lphy.base.evolution.eigensystem.*;
import lphy.base.evolution.substitutionmodel.RateMatrix;
import lphy.base.evolution.substitutionmodel.SubstModelParamNames;
import lphy.core.logger.LoggerUtils;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import lphy.core.model.datatype.DoubleArray2DValue;

import java.util.Arrays;

public class NonReversible extends RateMatrix {

    private Value<Boolean> symmetric;
    private static final double BRANCH_LENGTH = 1e6;

    public static final String indicatorsParamName = "indicators";
    protected static final String ratesParamName = SubstModelParamNames.RatesParamName;
    protected static final String freqParamName = SubstModelParamNames.FreqParamName;
    public static final String symmetricParamName = "symmetric";

    public NonReversible(@ParameterInfo(name = ratesParamName, narrativeName = "relative rates", description = "the relative rates of the substitution process.") Value<Double[]> rates,
                         @ParameterInfo(name = freqParamName, narrativeName = "base frequencies", description = "the base frequencies.", optional = true) Value<Double[]> freq,
                         @ParameterInfo(name = indicatorsParamName, narrativeName = "rate indicators", description = "a boolean for each rate to indicate the presence or absence of transition matrix entries", optional=true) Value<Boolean[]> indicators,
                         @ParameterInfo(name = RateMatrix.meanRateParamName, description = "the mean rate of the process. default 1.0", optional=true) Value<Number> meanRate,
                         @ParameterInfo(name = symmetricParamName, narrativeName = "", description = "", optional = true) Value<Boolean> symmetric) {
        super(meanRate);

        if (symmetric == null) symmetric = new Value(null,false);
        setParam(symmetricParamName, symmetric);

        int numStates = getNumStates(rates, symmetric);
        setParam(ratesParamName, rates);

        if (freq == null) {
            double f = (double) 1 /numStates;
            Double[] freqs = new Double[numStates];
            Arrays.fill(freqs, f);
            freq = new Value(null, freqs);
        }
        setParam(freqParamName, freq);
        if (indicators != null) {
            if (indicators.value().length != rates.value().length)
                throw new IllegalArgumentException("Indicators must have same number of dimensions as rates.");
        } else {
            Boolean[] b = new Boolean[rates.value().length];
            Arrays.fill(b, true);
            indicators = new Value(null, b);
        }
        setParam(indicatorsParamName, indicators);


    }

    private static int getNumStates(Value<Double[]> rates, Value<Boolean> symmetric) {
        int numStates;
        if (symmetric.value()) {
            double root = (-1 - Math.sqrt(1+8* rates.value().length))/2;
            if (root - (int) root < 1e-6) numStates = (int) Math.abs(root);
            else throw new IllegalArgumentException("Rates must have have (n²-n)/2 number of dimensions as frequencies.");
        } else {
            double root = (-1 - Math.sqrt(1+4* rates.value().length))/2;
            if (root - (int) root < 1e-6) numStates = (int) Math.abs(root);
            else throw new IllegalArgumentException("Rates must have have (n²-n) number of dimensions as frequencies.");
        }
        return numStates;
    }

    @GeneratorInfo(name = "nonReversible", verbClause = "is", narrativeName = "Estimated nonreversible substitution model",
            category = GeneratorCategory.RATE_MATRIX,
            description = "A customisable rate matrix using stochastic variable selection for nonreversible substitution.")
    public Value<Double[][]> apply() {
        Double[][] Q = getQ();
        return new DoubleArray2DValue(Q, this);
    }

    @Override
    public boolean canReturnComplexDiagonalization() { return !getSymmetric().value(); }

    protected Double[][] getQ() {
        Double[] r = getRates().value();
        Double[] f = getFreq().value();
        Boolean[] b = getIndicators().value();
        Boolean sym = getSymmetric().value();

        int numStates = f.length;
        Double[][] Q = new Double[numStates][numStates];
        int x = 0;
        for (int i = 0; i < numStates; i++) {
            Q[i][i] = 0.0;
            if (!sym) {
                for (int j = 0; j < i; j++) {
                    Q[i][j] = b[x] ? r[x] : Double.valueOf(0.0);
                    Q[i][i] -= Q[i][j];
                    x++;
                }
                for (int j = i + 1; j < numStates; j++) {
                    Q[i][j] = b[x] ? r[x] : Double.valueOf(0.0);
                    Q[i][i] -= Q[i][j];
                    x++;
                }
            } else {
                for (int j = i + 1; j < numStates; j++) {
                    if (b[x]) {
                        Q[i][j] = r[x] * f[j];
                        Q[j][i] = r[x] * f[i];
                    } else {
                        Q[i][j] = 0.0;
                        Q[j][i] = 0.0;
                    }
                    x++;
                }
                for (int j = 0; j < numStates; j++) {
                    if (i==j) continue;
                    Q[i][i] -= Q[i][j];
                }
            }
            if (Q[i][i]==0) LoggerUtils.log.severe("Empty row." +
                    "Rate indicators must have at least one true per row. " +
                    "Try increasing p in Bernoulli distribution, or using ConnectedSVS to ensure a connected graph.");
        }
        for (int i = 0; i < numStates; i++) {
            double test = 0.0;
            for (int j = 0; j < numStates; j++) {
                if (i==j) continue;
                test += Q[j][i];
            }
            if (test == 0.0) LoggerUtils.log.severe("Empty column." +
                    "Rate indicators must have at least one true per column." +
                    "Try increasing p in Bernoulli distribution, or using ConnectedSVS to ensure a connected graph.");
        }
        if (!sym) f = getEquilibriumFrequencies(numStates, Q);
        double subst = 0.0;
        for (int i = 0; i < Q.length; i++) {
            subst += -Q[i][i] * f[i];
        }

        for (int i = 0; i < Q.length; i++) {
            for (int j = 0; j < Q.length; j++) {
                Q[i][j] = Q[i][j] / subst;
            }
        }
        return Q;
    }

    private Double[] getEquilibriumFrequencies(int numStates, Double[][] Qm) {
        double temp;

        EigenSystem complexEigenSystem = new ComplexColtEigenSystem(numStates);
        EigenDecompositionExt complexDecomposition = complexEigenSystem.decomposeMatrix(Qm);
        double[] Eval = complexDecomposition.getEigenValues();
        double[] evec = complexDecomposition.getEigenVectors();
        double[] ievc = complexDecomposition.getInverseEigenVectors();

        double[][] iexp = new double[numStates][numStates];
        double[] EvalImag = new double[numStates];
        double[][] transProbs = new double[numStates][numStates];

        System.arraycopy(Eval, numStates, EvalImag, 0, numStates);
        for (int i = 0; i < numStates; i++) {
            if (EvalImag[i] == 0) {
                // 1x1 block
                temp = Math.exp(BRANCH_LENGTH * Eval[i]);
                for (int j = 0; j < numStates; j++) {
                    iexp[i][j] = ievc[i * numStates + j] * temp;
                }
            } else {
                // 2x2 conjugate block
                // If A is 2x2 with complex conjugate pair eigenvalues a +/- bi, then
                // exp(At) = exp(at)*( cos(bt)I + \frac{sin(bt)}{b}(A - aI)).
                int i2 = i + 1;
                double b = EvalImag[i];
                double expat = Math.exp(BRANCH_LENGTH * Eval[i]);
                double expatcosbt = expat * Math.cos(BRANCH_LENGTH * b);
                double expatsinbt = expat * Math.sin(BRANCH_LENGTH * b);

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
                for (int k = 0; k < numStates; k++) {
                    temp += evec[i * numStates + k] * iexp[k][j];
                }
                transProbs[i][j] = Math.abs(temp);
            }
        }
        Double[] freqs = new Double[transProbs.length];
        for (int i = 0; i < freqs.length; i++) {
            freqs[i] = transProbs[0][i];
            for (int j = 1; j < freqs.length; j++) {
                if (Math.abs(transProbs[0][i] - transProbs[j][i]) > 1e-12) {
                    System.out.println("WARNING: branch length used to get equilibrium distribution was not long enough!");
                }

            }
        }
        return freqs;
    }

    public Value<Double[]> getRates() {
        return getParams().get(ratesParamName);
    }

    public Value<Double[]> getFreq() {
        return getParams().get(freqParamName);
    }

    public Value<Boolean[]> getIndicators() {
        return getParams().get(indicatorsParamName);
    }

    public Value<Boolean> getSymmetric() {
        return getParams().get(symmetricParamName);
    }

}