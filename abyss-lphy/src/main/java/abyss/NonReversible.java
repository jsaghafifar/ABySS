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

import static abyss.ComputeEquilibrium.computeEquilibrium;

/**
 * Q matrix with optional nonreversibility and SVS rate indicators
 * @author Jasmine Saghafifar
 */
public class NonReversible extends RateMatrix {
    public static final String indicatorsParamName = "indicators";
    protected static final String ratesParamName = SubstModelParamNames.RatesParamName;
    protected static final String freqParamName = SubstModelParamNames.FreqParamName;
    public static final String symmetricParamName = "symmetric";

    public NonReversible(@ParameterInfo(name = ratesParamName, narrativeName = "relative rates", description = "the relative rates of the substitution process.") Value<Double[]> rates,
                         @ParameterInfo(name = freqParamName, narrativeName = "base frequencies", description = "the base frequencies.", optional = true) Value<Double[]> freq,
                         @ParameterInfo(name = indicatorsParamName, narrativeName = "rate indicators", description = "a boolean for each rate to indicate the presence or absence of transition matrix entries", optional=true) Value<Boolean[]> indicators,
                         @ParameterInfo(name = RateMatrix.meanRateParamName, description = "the mean rate of the process. default 1.0", optional=true) Value<Number> meanRate,
                         @ParameterInfo(name = symmetricParamName, narrativeName = "Q matrix symmetry", description = "whether Q is time-reversible. default false", optional = true) Value<Boolean> symmetric) {
        super(meanRate);


        setParam(ratesParamName, rates);

        if (symmetric == null) symmetric = new Value(null,false);
        setParam(symmetricParamName, symmetric);

        Boolean sym = symmetric.value();
        if (sym) {
            if (freq == null) throw new IllegalArgumentException("Frequencies must be specified when reversible Q matrix");
            else setParam(freqParamName, freq);
        } else if (freq != null) throw new IllegalArgumentException("Frequencies cannot be specified when nonreversible Q matrix");

        if (indicators != null) {
            int indicatorsLength = indicators.value().length;
            int ratesLength = rates.value().length;
            if (indicatorsLength != ratesLength)
                throw new IllegalArgumentException("Indicators must have same number of dimensions as rates.");
            if (!SVS.Utils.isStronglyConnected(indicators.value(),getNumStates(indicatorsLength,sym),sym))
                LoggerUtils.log.severe("Q matrix is not strongly connected. " +
                        "Try adjusting Bernoulli or using ConnectedSVS to ensure connectivity.");
            setParam(indicatorsParamName, indicators);
        }


    }

    private static int getNumStates(int len, Boolean symmetric) {
        int numStates;
        if (symmetric) {
            double root = (-1 - Math.sqrt(1+8* len))/2;
            if (root - (int) root < 1e-6) numStates = (int) Math.abs(root);
            else throw new IllegalArgumentException("Rates must have have (n²-n)/2 number of dimensions as frequencies.");
        } else {
            double root = (-1 - Math.sqrt(1+4* len))/2;
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

    private static void setupRelativeRates(Double[][] Q, Double[] r, Boolean[] b, int numStates) {
        int x = 0;
        for (int i = 0; i < numStates; i++) {
            Q[i][i] = 0.0;
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
        }
    }

    protected Double[][] getQ() {
        Double[] r = getRates().value();
        Boolean sym = getSymmetric().value();
        int numStates = getNumStates(r.length, sym);

        Boolean[] b = new Boolean[r.length];
        if (getIndicators() != null) b = getIndicators().value();
        else Arrays.fill(b, true);

        Double[] f;
        Double[][] Q = new Double[numStates][numStates];
        if (!sym) {
            setupRelativeRates(Q, r, b, numStates);
            f = computeEquilibrium(Q);
        } else {
            f = getFreq().value();
            int x = 0;
            for (int i = 0; i < numStates; i++) {
                Q[i][i] = 0.0;
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
                    if (i == j) continue;
                    Q[i][i] -= Q[i][j];
                }
                if (Q[i][i] == 0) LoggerUtils.log.severe("Empty row." +
                        "Rate indicators must have at least one true per row. " +
                        "Try increasing p in Bernoulli distribution, or using ConnectedSVS to ensure a connected graph.");
            }
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