package abyss;

import lphy.base.distribution.ParametricDistribution;
import lphy.base.evolution.substitutionmodel.SubstModelParamNames;
import lphy.core.logger.LoggerUtils;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.util.FastMath;

import java.util.Map;
import java.util.TreeMap;

/**
 * A rate-informed Bernoulli distribution
 * for obtaining connected SVS rate indicators
 * @author jsaghafifar
 */
public class ConnectedSVS extends ParametricDistribution<Boolean[]> implements SVS {
    private Value<Double[]> rates;
    private Value<Number> a;
    private Value<Number> k;

    private static final int MAX_TRIES = 10000;
    private double[] p;

    public static final String ratesParamName = SubstModelParamNames.RatesParamName;
    public static final String scaleParamName = "scale";
    public static final String shapeParamName = "shape";

    public ConnectedSVS(@ParameterInfo(name = ratesParamName, narrativeName = "rates", description = "normalised nonreversible empirical rate matrix") Value<Double[]> rates,
                        @ParameterInfo(name = scaleParamName, narrativeName = "scale", description = "scale of active rates, default 1", optional=true) Value<Number> a,
                        @ParameterInfo(name = shapeParamName, narrativeName = "shape", description = "indicators' sensitivity to rates, default 1", optional=true) Value<Number> k) {
        super();
        this.rates = rates;
        if (a != null) {
            if (a.value().doubleValue() > 0) {
                this.a = a;
            } else throw new IllegalArgumentException(scaleParamName+" must be > 0.");
        } else this.a = new Value(1.0,null);

        if (k != null) {
            if (k.value().doubleValue() > 0) {
                this.k = k;
            } else throw new IllegalArgumentException(shapeParamName+" must be > 0.");
        } else this.k = new Value(1.0,null);

        constructDistribution(random);
    }

    @Override
    protected void constructDistribution(RandomGenerator random) { }

    @GeneratorInfo(name = "ConnectedSVS", verbClause = "has", narrativeName = "Bernoulli distribution prior",
            category = GeneratorCategory.PRIOR,
            description = "rate-informed Bernoulli distribution with checks to ensure connected graphs for SVS rate matrices")
    public RandomVariable<Boolean[]> sample() {
        double a = getScale().value().doubleValue();
        double k = getShape().value().doubleValue();
        Double[] r = getRates().value();
        double sum = 0;
        for (int i = 0; i < r.length; i++) {
            sum += r[i];
        }
        if (Math.abs(sum - 1.0) > 1e-6) throw new IllegalArgumentException("Rates must be normalised to sum to 1.");

        int n;
        double root = (-1 - Math.sqrt(1+4*r.length))/2;
        if (root % 1 == 0) {
            n = (int) Math.abs(root);
        } else throw new IllegalArgumentException(
                "Wrong number of rates. Must be such that n(n-1) (n = number of states). " +
                        "e.g. 4 nucleotide states = 12 rates.");

        Boolean[] b = new Boolean[r.length];
        boolean connectedGraph = false;
        int iter = 0;

        p = new double[r.length];
        for (int i = 0; i < p.length; i++) {
            double x = Math.log(p.length * r[i]);
            x *= -k; // scaled
            p[i] = a/(a + FastMath.exp(x));
        }
        while (!connectedGraph) {
            int successes = 0;
            for (int i = 0; i < p.length; i++) {
                double U = random.nextDouble();
                if (p[i] > U) {
                    b[i] = Boolean.TRUE;
                    successes++;
                } else b[i] = Boolean.FALSE;
            }
            if (successes >= n &&
                    SVS.Utils.isStronglyConnected(b, n, false) &&
                    SVS.Utils.connectedAndWellConditioned(p)) {
                connectedGraph = true;
            }
            iter++;
            if (!connectedGraph && iter % (MAX_TRIES/10)==0 && iter < MAX_TRIES) {
                LoggerUtils.log.warning("Graph not connected after " + iter + " iterations. Consider increasing "+scaleParamName+".");
            }
            if (!connectedGraph && iter >= MAX_TRIES) {
                LoggerUtils.log.severe("Max iterations exceeded! Try increasing "+scaleParamName+" to keep the graph connected.");
                throw new RuntimeException("Max iterations exceeded.");
            }
        }
        return new RandomVariable<>(null, b, this);
    }

    public double logDensity(Boolean[] successes) { //TODO
        double logP = 0.0;
        for (int i = 0; i < successes.length; i++) {
            logP += successes[i] ? Math.log(p[i]) : Math.log(1-p[i]);
        }
        return logP;
    }

    @Override
    public Map<String, Value> getParams() {
        return new TreeMap<>() {{
            put(ratesParamName, rates);
            put(scaleParamName, a);
            put(shapeParamName, k);
        }};
    }

    @Override
    public void setParam(String paramName, Value value) {
        switch (paramName) {
            case ratesParamName:
                rates = value;
                break;
            case scaleParamName:
                a = value;
                break;
            case shapeParamName:
                k = value;
                break;
            default:
                throw new RuntimeException("Unrecognised parameter name: " + paramName);
        }

        super.setParam(paramName, value); // constructDistribution
    }

    public Value<Double[]> getRates() {
        return getParams().get(ratesParamName);
    }

    public Value<Number> getScale() {
        return getParams().get(scaleParamName);
    }

    public Value<Number> getShape() {
        return getParams().get(shapeParamName);
    }


}
