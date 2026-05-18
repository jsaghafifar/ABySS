package abyss;

import lphy.base.evolution.substitutionmodel.RateMatrix;
import lphy.base.evolution.substitutionmodel.SubstModelParamNames;
import lphy.core.logger.LoggerUtils;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import lphy.core.model.datatype.DoubleArray2DValue;

/**
 * Q matrix with SVS rate indicators and optional nonreversibility
 * Based on bugged version of BEASTClassic model
 * @author Jasmine Saghafifar
 */
public class ClassicNonReversible extends RateMatrix {
    public static final String indicatorsParamName = "indicators";
    protected static final String ratesParamName = SubstModelParamNames.RatesParamName;
    protected static final String freqParamName = SubstModelParamNames.FreqParamName;
    public static final String symmetricParamName = "symmetric";

    public ClassicNonReversible(@ParameterInfo(name = ratesParamName, narrativeName = "relative rates", description = "the relative rates of the substitution process.") Value<Double[]> rates,
                                @ParameterInfo(name = freqParamName, narrativeName = "base frequencies", description = "the base frequencies.") Value<Double[]> freq,
                                @ParameterInfo(name = indicatorsParamName, narrativeName = "rate indicators", description = "a boolean for each rate to indicate the presence or absence of transition matrix entries") Value<Boolean[]> indicators,
                                @ParameterInfo(name = RateMatrix.meanRateParamName, description = "the mean rate of the process. default 1.0", optional=true) Value<Number> meanRate,
                                @ParameterInfo(name = symmetricParamName, narrativeName = "Q matrix symmetry", description = "whether Q is time-reversible. default false", optional = true) Value<Boolean> symmetric) {
        super(meanRate);


        setParam(ratesParamName, rates);
        setParam(freqParamName, freq);
        if (symmetric == null) symmetric = new Value(null,false);
        if (!symmetric.value())
            LoggerUtils.log.warning("classicNonReversible is intentionally left bugged! " +
                    "Use nonReversible for valid time-nonreversible substitution simulation.");

        setParam(symmetricParamName, symmetric);

        if (indicators.value().length != rates.value().length)
            throw new IllegalArgumentException("Indicators must have same number of dimensions as rates.");
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

    @GeneratorInfo(name = "classicNonReversible", verbClause = "is", narrativeName = "Bugged nonreversible substitution model",
            category = GeneratorCategory.RATE_MATRIX,
            description = "Bugged, maintained for presenting implementation process. " +
                    "A cautionary tale in why frequencies should not inform time-nonreversible substitution.")
    public Value<Double[][]> apply() {
        Double[][] Q = getQ();
        return new DoubleArray2DValue(Q, this);
    }

    @Override
    public boolean canReturnComplexDiagonalization() { return true; }

    private static void setupRelativeRates(Double[][] Q, Double[] r, Boolean[] b, int numStates, Boolean sym) {
        int x = 0;
        for (int i = 0; i < numStates; i++) {
            Q[i][i] = 0.0;
            if (sym) {
                for (int j = 0; j < i; j++) {
                    Q[i][j] = b[x] ? r[x] : Double.valueOf(0.0);
                    Q[j][i] = b[x] ? r[x] : Double.valueOf(0.0);
                    x++;
                }
            } else {
                for (int j = 0; j < i; j++) {
                    Q[i][j] = b[x] ? r[x] : Double.valueOf(0.0);
                    x++;
                }
                for (int j = i + 1; j < numStates; j++) {
                    Q[i][j] = b[x] ? r[x] : Double.valueOf(0.0);
                    x++;
                }
            }
        }
    }

    protected Double[][] getQ() {
        int numStates = getNumStates(getRates(), getSymmetric());
        Double[] r = getRates().value();
        Boolean[] b = getIndicators().value();

        Double[] f = getFreq().value();
        Double[][] Q = new Double[numStates][numStates];
        setupRelativeRates(Q, r, b, numStates, getSymmetric().value());

        for (int i = 0; i < numStates; i++) {
            for (int j = 0; j < numStates; j++) {
                if (i == j) continue;
                Q[i][j] *= f[j];
                Q[i][i] -= Q[i][j];
            }
        }
        normalize(f, Q, totalRateDefault1());
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