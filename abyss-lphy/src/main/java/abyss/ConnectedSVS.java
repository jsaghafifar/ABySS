package abyss;

import lphy.base.distribution.ParametricDistribution;
import lphy.core.logger.LoggerUtils;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.util.FastMath;

import java.util.Arrays;
import java.util.Map;
import java.util.Objects;
import java.util.TreeMap;

/**
 * A rate-informed Bernoulli distribution
 * for obtaining connected SVS rate indicators
 * @author jsaghafifar
 */
public class ConnectedSVS extends ParametricDistribution<Boolean[]> implements SVS {
    private Value<Integer> numStates;
    private Value<Double[][]> Q;
    private Value<Number> a;
    private Value<Number> k;
    private Value<Boolean> symmetric;

    private static final int MAX_TRIES = 10000;
    private double[] p;

    public static final String numStatesParamName = "states";
    public static final String empiricalQParamName = "empiricalQ";
    public static final String sensitivityParamName = "sens";
    public static final String scaleParamName = "scale";
    public static final String symmetricParamName = "symmetric";

    public ConnectedSVS(@ParameterInfo(name = numStatesParamName, narrativeName = "states", description = "number of states in Q matrix") Value<Integer> numStates,
                        @ParameterInfo(name = empiricalQParamName, narrativeName = "empiricalQ", description = "normalised empirical Q matrix", optional = true) Value<Double[][]> Q,
                        @ParameterInfo(name = sensitivityParamName, narrativeName = "sens", description = "indicators' sensitivity to rates, default 1", optional=true) Value<Number> k,
                        @ParameterInfo(name = scaleParamName, narrativeName = "scale", description = "scale of active rates, default 1", optional=true) Value<Number> a,
                        @ParameterInfo(name = symmetricParamName, narrativeName = "symmetric", description = "whether rate matrix is symmetric", optional = true) Value<Boolean> symmetric) {
        super();
        this.numStates = numStates;

        if (Q != null) this.Q = Q;

        if (k != null) {
            if (Q != null) {
                if (k.value().doubleValue() > 0) this.k = k;
                else throw new IllegalArgumentException(sensitivityParamName +" must be over 0.");
            } else throw new IllegalArgumentException(sensitivityParamName +" should not be specified when "+
                    empiricalQParamName+" is not specified.");
        } else this.k = new Value(1.0,null);

        if (a != null) {
            if (a.value().doubleValue() > 0) this.a = a;
            else throw new IllegalArgumentException(scaleParamName +" must be over 0.");
        } else this.a = new Value(1.0,null);

        setParam(symmetricParamName, Objects.requireNonNullElseGet(symmetric, () -> new Value(null, false)));

        constructDistribution(random);
    }

    @Override
    protected void constructDistribution(RandomGenerator random) { }

    @GeneratorInfo(name = "ConnectedSVS", verbClause = "has", narrativeName = "Bernoulli distribution prior",
            category = GeneratorCategory.PRIOR,
            description = "rate-informed Bernoulli distribution with checks to ensure connected graphs for SVS rate matrices")
    public RandomVariable<Boolean[]> sample() {
        double k = getSensitivity().value().doubleValue();
        double a = getScale().value().doubleValue();
        int n = getNumStates().value();
        Double[][] q = new Double[n][n];
        if (getEmpiricalQ() != null) {
            q = getEmpiricalQ().value();
        } else {
            // uninformed Bernoulli: each indicator will have equal probability of being active
            for (int i = 0; i < n; i++) {
                Arrays.fill(q[i], 1.0);
            }
        }
        Double[] r = getQValues(q);
        Boolean[] b = new Boolean[r.length];
        Boolean sym = getSymmetric().value();
        
        boolean connectedGraph = false;
        int iter = 0;

        p = new double[r.length];
        for (int i = 0; i < p.length; i++) {
            double x = Math.log(p.length * r[i]);
            x *= -k; // sensitivity to empirical Q input is applied
            p[i] = a/(a + FastMath.exp(x)); // probability of indicator scaled
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
                    SVS.Utils.isStronglyConnected(b, n, sym)) {
                connectedGraph = true;
            }
            iter++;
            if (!connectedGraph && iter % (MAX_TRIES/10)==0 && iter < MAX_TRIES) {
                LoggerUtils.log.warning("Graph not connected after " + iter + " iterations. " +
                        "Consider increasing "+ scaleParamName +".");
            }
            if (!connectedGraph && iter >= MAX_TRIES) {
                LoggerUtils.log.severe("Max iterations exceeded! " +
                        "Try increasing "+ scaleParamName +" to keep the graph connected.");
                throw new RuntimeException("Max iterations exceeded.");
            }
        }
        return new RandomVariable<>(null, b, this);
    }

    public double logDensity(Boolean[] successes) {
        double logP = 0.0;
        for (int i = 0; i < successes.length; i++) {
            logP += successes[i] ? Math.log(p[i]) : Math.log(1-p[i]);
        }
        return logP;
    }

    @Override
    public Map<String, Value> getParams() {
        return new TreeMap<>() {{
            put(numStatesParamName, numStates);
            if (Q != null) put(empiricalQParamName, Q);
            put(sensitivityParamName, k);
            put(scaleParamName, a);
            put(symmetricParamName, symmetric);
        }};
    }

    @Override
    public void setParam(String paramName, Value value) {
        switch (paramName) {
            case numStatesParamName:
                numStates = value;
                break;
            case empiricalQParamName:
                Q = value;
                break;
            case sensitivityParamName:
                k = value;
                break;
            case scaleParamName:
                a = value;
                break;
            case symmetricParamName:
                symmetric = value;
                break;
            default:
                throw new RuntimeException("Unrecognised parameter name: " + paramName);
        }

        super.setParam(paramName, value); // constructDistribution
    }

    public Value<Double[][]> getEmpiricalQ() {
        if (Q != null) return getParams().get(empiricalQParamName);
        else return null;
    }

    public Value<Integer> getNumStates() {
        return getParams().get(numStatesParamName);
    }

    public Value<Number> getSensitivity() {
        return getParams().get(sensitivityParamName);
    }

    public Value<Number> getScale() {
        return getParams().get(scaleParamName);
    }

    public Value<Boolean> getSymmetric() {
        return getParams().get(symmetricParamName);
    }

    private Double[] getQValues(Double[][] q) {
        Boolean sym = getSymmetric().value();
        int symQ = sym ? 2 : 1;
        int n = getNumStates().value();

        // construct rate matrix from empirical Q input
        Double[] r = new Double[(n * n - n)/symQ];
        int x = 0;
        double sum = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i==j) continue;
                if (sym & j>i) continue;
                r[x] = q[i][j];
                sum += r[x];
                x++;
            }
        }

        // normalise rates
        x = 0;
        for (int i = 0; i < r.length; i++) {
            r[x] /= sum;
            x++;
        }

        return r;
    }


}
