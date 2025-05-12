package abyss;

import lphy.base.distribution.ParametricDistribution;
import lphy.base.distribution.RandomBooleanArray;
import lphy.core.logger.LoggerUtils;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import lphy.core.model.datatype.IntegerValue;
import lphy.core.vectorization.IID;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.random.RandomGenerator;

import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

import static lphy.base.distribution.DistributionConstants.pParamName;

/**
 * A Binomial process ensuring a connected
 * graph when applied to rate matrices.
 *
 */
public class ConnectedGraphBinomial extends ParametricDistribution<Boolean[]> {
    private Value<Double> p;
    private Value<Integer> numStates;

    public static final String numParamName = "numStates";//

    BinomialDistribution binomialDistribution;

    public ConnectedGraphBinomial(@ParameterInfo(name = pParamName, description = "the probability of success.") Value<Double> p,
                                  @ParameterInfo(name = numParamName, description = "the number of states.") Value<Integer> numStates){
        super();
        this.numStates = numStates;
        if (p.value() >= ((double) (numStates.value() - 2) / (numStates.value() - 1))) {
            LoggerUtils.log.severe(pParamName + " is too high. Must be < (numStates-2)/(numStates-1)");
        } else {
            this.p = p;
        }

        constructDistribution(random);
    }

    @Override
    protected void constructDistribution(RandomGenerator random) {
        int numRates = numStates.value() * numStates.value() - numStates.value();
        binomialDistribution = new BinomialDistribution(random, numRates, p.value());
    }

    @GeneratorInfo(name = "ConnectedBinomial", verbClause = "has", narrativeName = "connected graph prior on boolean array",
            category = GeneratorCategory.PRIOR,
            description = "The binomial distribution with success probability of p, producing a connected graph of booleans, size n².")
    public RandomVariable<Boolean[]> sample() {
        int numRates = numStates.value() * numStates.value() - numStates.value();
        int minSuccesses = numStates.value();
        Boolean[] b = new Boolean[numRates];
        boolean connectedGraph = false;
        int iter = 0;

        while (!connectedGraph) {
            double[] p = new double[numRates - minSuccesses];
            double probSum = 0.0;

            for (int i = 0; i < p.length; i++) {
                p[i] = binomialDistribution.probability(i + minSuccesses);
                probSum += p[i];
            }
            if (probSum > 0.0) {
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
                RandomBooleanArray randomBooleanArray = new RandomBooleanArray(new Value<>(null, numRates),
                        new IntegerValue(index + minSuccesses, null));
                b = randomBooleanArray.sample().value();
            } else {
                LoggerUtils.log.severe(pParamName + " is too high.");
            }

            boolean rowTest = false;
            boolean colTest = false;
            for (int i = 0; i < numStates.value(); i++) {
                rowTest = false;

                for (int j = 0; j < (numStates.value()-1); j++) {
                    if (b[i * (numStates.value()-1) + j]) {
                        rowTest = true;
                        break;
                    }
                }
                if (!rowTest) break;
            }
            for (int j = 0; j < (numStates.value()-1); j++) {
                colTest = false;
                for (int i = 0; i < numStates.value(); i++) {
                    if (b[i * (numStates.value()-1) + j]) {
                        colTest = true;
                        break;
                    }
                }
                if (!colTest) break;
            }
            if (rowTest && colTest) connectedGraph = true;
            iter++;
            if (!connectedGraph && iter % 10==0 && iter < 100) {
                LoggerUtils.log.warning("Graph not connected after " + iter + " iterations. Consider increasing "+
                         pParamName + ".");
            }
            if (!connectedGraph && iter >= 100) {
                LoggerUtils.log.severe("Max iterations exceeded! Try increasing " +
                        pParamName + " to keep graph connected.");
                throw new RuntimeException("Max iterations exceeded.");
            }

        }

        return new RandomVariable<>(null, b, this);
    }

//    private int hammingWeight(Boolean[] b) {
//        int sum = 0;
//        for (Boolean i : b) {
//            if (i) sum += 1;
//        }
//        return sum;
//    }
//
//    private Boolean[] bernoulli(double p, int n) {
//        Boolean[] successes = new Boolean[n * n];
//        for (int i = 0; i < successes.length; i++) {
//            successes[i] = (random.nextDouble() < p);
//        }
//        return successes;
//    }

    public double logDensity(Boolean[] successes) {
        int sum = 0;
        int numRates = numStates.value() * numStates.value() - numStates.value();
        for (int i = 0; i < numRates; i ++) {
            if (successes[i]) {
                sum ++;
            }
        }

        double logFactorialn=0, logFactorial1=0, logFactorial2=0;
        for (int i = 2; i <= numRates;  i ++) logFactorialn += Math.log(i);
        for (int i = 2; i <= sum; i ++) logFactorial1 += Math.log(i);
        for (int i = 2; i <= (numRates-sum); i ++) logFactorial2 += Math.log(i);

        double logP = logFactorialn - (logFactorial1 + logFactorial2) + sum*Math.log(p.value()) +
                (numRates-sum)*Math.log(1-p.value());

        return logP;
    }

    @Override
    public Map<String, Value> getParams() {
        return new TreeMap<>() {{
            put(pParamName, p);
            put(numParamName, numStates);
        }};
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(pParamName)) p = value;
        else if (paramName.equals(numParamName)) numStates = value;
        else throw new RuntimeException("Unrecognised parameter name: " + paramName);

        super.setParam(paramName, value); // constructDistribution
    }

//    public void setSuccessProbability(double p) {
//        this.p.setValue(p);
//    }

    public Value<Double> getP() {
        return getParams().get(pParamName);
    }

}
