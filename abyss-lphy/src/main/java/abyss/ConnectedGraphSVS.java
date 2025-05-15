package abyss;

import lphy.base.distribution.ParametricDistribution;
import lphy.base.distribution.RandomBooleanArray;
import lphy.base.evolution.substitutionmodel.SubstModelParamNames;
import lphy.core.logger.LoggerUtils;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import lphy.core.model.datatype.IntegerValue;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.util.FastMath;

import java.util.Map;
import java.util.TreeMap;

/**
 *
 */
public class ConnectedGraphSVS extends ParametricDistribution<Boolean[]> {
    private Value<Double[]> rates;
    private Value<Double> k;
    private Value<Integer> numStates;

    private static final int MAX_TRIES = 10000;

    public static final String ratesParamName = SubstModelParamNames.RatesParamName;
    public static final String scaleParamName = "scale";
    public static final String numStatesParamName = "numStates";

    private final int NUM_RATES;

    public ConnectedGraphSVS(@ParameterInfo(name = ratesParamName, narrativeName = "rates", description = "") Value<Double[]> rates,
                             @ParameterInfo(name = scaleParamName, narrativeName = "scale", description = "") Value<Double> k,
                             @ParameterInfo(name = numStatesParamName, narrativeName = "numStates", description = "") Value<Integer> numStates) {
        super();
        this.rates = rates;
        this.k = k;
        this.numStates = numStates;
        this.NUM_RATES = numStates.value() * numStates.value() - numStates.value();

        constructDistribution(random);
    }

    @Override
    protected void constructDistribution(RandomGenerator random) {
        //TODO
    }

    @GeneratorInfo(name = "ConnectedSVS", verbClause = "has", narrativeName = "distribution prior",
            category = GeneratorCategory.PRIOR,
            description = "")
    public RandomVariable<Boolean[]> sample() {
        double k = getScale().value();
        Double[] r = getRates().value();
        int n = getNumStates().value();

        Boolean[] b = new Boolean[NUM_RATES];
        boolean connectedGraph = false;
        int iter = 0;

        while (!connectedGraph) {
            double[] p = new double[NUM_RATES];
            double probSum = 0.0;
            for (int i = 0; i < p.length; i++) {
                double probability = 1/(1+Math.exp(-k * Math.log(r[i]))); //TODO fix
                if (probability == Double.NEGATIVE_INFINITY) {
                    p[i] = 0.0F;
                } else p[i] = FastMath.exp(probability);
                probSum += p[i];
            }
            double U = random.nextDouble() * probSum;
            probSum = 0.0;
            int index = 0;
            for (int i = 0; i < p.length; i++) {
                probSum += p[i];
                if (probSum > U) {
                    index = i;
                    break;
                }
            }
            RandomBooleanArray randomBooleanArray = new RandomBooleanArray(
                    new Value<>(null,NUM_RATES), new IntegerValue(index, null));
            b = randomBooleanArray.sample().value();

            boolean rowTest = false;
            boolean colTest = false;
            for (int i = 0; i < n; i++) {
                rowTest = false;

                for (int j = 0; j < (n-1); j++) {
                    if (b[i * (n-1) + j]) {
                        rowTest = true;
                        break;
                    }
                }
                if (!rowTest) break;
            }
            for (int j = 0; j < (n-1); j++) {
                colTest = false;
                for (int i = 0; i < n; i++) {
                    if (b[i * (n-1) + j]) {
                        colTest = true;
                        break;
                    }
                }
                if (!colTest) break;
            }
            if (rowTest && colTest) connectedGraph = true;
            iter++;
            if (!connectedGraph && iter % (MAX_TRIES/10)==0 && iter < MAX_TRIES) {
                LoggerUtils.log.warning("Graph not connected after " + iter + " iterations. Consider increasing "+
                        ".");
            }
            if (!connectedGraph && iter >= MAX_TRIES) {
                LoggerUtils.log.severe("Max iterations exceeded! Try increasing " +
                        " to keep graph connected.");
                throw new RuntimeException("Max iterations exceeded.");
            }
        }
        return new RandomVariable<>(null, b, this);
    }

    public double logDensity(Boolean[] successes) {
        //TODO
        throw new UnsupportedOperationException("TODO");
    }

    @Override
    public Map<String, Value> getParams() {
        return new TreeMap<>() {{
            put(ratesParamName, rates);
            put(scaleParamName, k);
            put(numStatesParamName, numStates);
        }};
    }

    @Override
    public void setParam(String paramName, Value value) {
        switch (paramName) {
            case ratesParamName:
                rates = value;
                break;
            case scaleParamName:
                k = value;
                break;
            case numStatesParamName:
                numStates = value;
                break;
            default:
                throw new RuntimeException("Unrecognised parameter name: " + paramName);
        }

        super.setParam(paramName, value); // constructDistribution
    }

    public Value<Double[]> getRates() {
        return getParams().get(ratesParamName);
    }

    public Value<Double> getScale() {
        return getParams().get(scaleParamName);
    }

    public Value<Integer> getNumStates() {
        return getParams().get(numStatesParamName);
    }


}
