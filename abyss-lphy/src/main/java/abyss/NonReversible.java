package abyss;

import lphy.base.evolution.substitutionmodel.RateMatrix;
import lphy.base.evolution.substitutionmodel.SubstModelParamNames;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import lphy.core.model.datatype.DoubleArray2DValue;

import java.util.Arrays;

public class NonReversible extends RateMatrix {

    protected static final String indicatorsParamName = "indicators";
    protected static final String ratesParamName = SubstModelParamNames.RatesParamName;
    protected static final String freqParamName = SubstModelParamNames.FreqParamName;

    public NonReversible(@ParameterInfo(name = ratesParamName, narrativeName = "relative rates", description = "the relative rates of the substitution process.") Value<Double[]> rates,
                         @ParameterInfo(name = freqParamName, narrativeName = "base frequencies", description = "the base frequencies.") Value<Double[]> freq,
                         @ParameterInfo(name = indicatorsParamName, narrativeName = "rate indicators", description = "a boolean for each rate to indicate the presence or absence of transition matrix entries", optional=true) Value<Boolean[]> indicators,
                         @ParameterInfo(name = RateMatrix.meanRateParamName, description = "the mean rate of the process. default 1.0", optional=true) Value<Number> meanRate) {
        super(meanRate);

        setParam(freqParamName, freq);
        int numStates = freq.value().length;

        if (rates.value().length != (numStates * numStates - numStates)) {
            throw new IllegalArgumentException("Rates must have have (n²-n) number of dimensions as frequencies.");
        }
        setParam(ratesParamName, rates);

        if (indicators != null) {
            if (indicators.value().length != rates.value().length)
                throw new IllegalArgumentException("Indicators must have same number of dimensions as rates.");
//            int x = 0;
//            for (int i = 0; i < numStates; i++) {
//                boolean test = false;
//                for (int j = 0; j < (numStates - 1); j++) {
//                    if (indicators.value()[x]) {
//                        test = true;
//                    } x++;
//                } if (!test) {
//                    throw new IllegalArgumentException("Rate indicators must have at least one true per row. " +
//                            "Try increasing p or minSuccesses in Bernoulli distribution, or resampling.");
//                }
//            }
        } else {
            Boolean[] b = new Boolean[rates.value().length];
            Arrays.fill(b, true);
            indicators = new Value("indicators", b);
        }
        setParam(indicatorsParamName, indicators);
    }

    @GeneratorInfo(name = "nonReversible", verbClause = "is", narrativeName = "Estimated nonreversible substitution model",
            category = GeneratorCategory.RATE_MATRIX,
            description = "A variable asymmetric rate matrix using stochastic variable selection for nonreversible substitution.")
    public Value<Double[][]> apply() {
        Double[][] Q = getQ();
        return new DoubleArray2DValue(Q, this);
    }

    @Override
    public boolean canReturnComplexDiagonalization() { return true; }

    protected Double[][] getQ() {
        Double[] r = getRates().value();
        Double[] f = getFreq().value();
        Boolean[] b = getIndicators().value();

        int numStates = f.length;
        Double[][] Q = new Double[numStates][numStates];
        int x = 0;
        for (int i = 0; i < numStates; i++) {
            Q[i][i] = 0.0;
            for (int j = 0; j < i; j++) {
                if (b[x]) {
                    Q[i][j] = r[x] * f[j];
                    Q[i][i] -= Q[i][j];
                } else Q[i][j] = 0.0;
                x++;
            }
            for (int j = i + 1; j < numStates; j++) {
                if (b[x]) {
                    Q[i][j] = r[x] * f[j];
                    Q[i][i] -= Q[i][j];
                } else Q[i][j] = 0.0;
                x++;
            }
            if (Q[i][i]==0) throw new IllegalArgumentException("Empty row." +
                    "Rate indicators must have at least one true per row. " +
                    "Try increasing p or minSuccesses in Bernoulli distribution, or resampling.");
        }
        for (int i = 0; i < numStates; i++) {
            double test = 0.0;
            for (int j = 0; j < numStates; j++) {
                if (i==j) continue;
                test += Q[j][i];
            }
            if (test == 0.0) throw new IllegalArgumentException("Empty column." +
                    "Rate indicators must have at least one true per column." +
                    "Try increasing p or minSuccesses in Bernoulli distribution, or resampling.");
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

}