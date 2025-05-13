package abyss;

import lphy.base.distribution.ParametricDistribution;
import lphy.base.math.MathUtils;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.Map;
import java.util.TreeMap;

/**
 * Dirichlet distribution prior informed by empirical rates/frequencies.
 */
public class InformedDirichlet extends ParametricDistribution<Double[]> {

    private Value<Number[]> inf;
    public static final String infParamName = "information";
    private Value<Number> omega;
    public static final String omegaParamName = "omega";
    private Value<Number> mean;
    private static final String meanParamName = "mean";


    public InformedDirichlet(
            @ParameterInfo(name=infParamName, narrativeName = "information",
                    description="informing concentration.") Value<Number[]> inf,
            @ParameterInfo(name=omegaParamName, narrativeName = "omega", description = "scale factor") Value<Number> omega,
            @ParameterInfo(name = meanParamName,
                    description = "The expected mean per element. By default, the sampled values sum to 1.",
                    optional = true) Value<Number> mean)
             {
        super();
        this.inf = inf;
        this.omega = omega;
        if (mean != null) {
            this.mean = mean;
        } else this.mean = new Value<>(null,1.0);

    }

    @Override
    protected void constructDistribution(RandomGenerator random) {  }

    @GeneratorInfo(name="InformedDirichlet", verbClause = "have", narrativeName = "Rate- or frequency-informed Dirichlet distribution prior",
            category = GeneratorCategory.PRIOR,
            description="The informed dirichlet probability distribution.")
    public RandomVariable<Double[]> sample() {

        double mean = getMean().value().doubleValue(); // means and omega fits where?
        double omega = getOmega().value().doubleValue();
        Double[] dirichlet = new Double[inf.value().length];
        double sum = 0.0;
        for (int i = 0; i < dirichlet.length; i++) {
            double val = MathUtils.randomGamma(inf.value()[i].doubleValue(), 1.0, random); // should (*omega) be here? instead
            dirichlet[i] = val * omega;
            sum += val;
        }

        for (int i = 0; i < dirichlet.length; i++) {
                dirichlet[i] /= sum; // (* mean) here?
            }


        return new RandomVariable<>("x", dirichlet, this);
    }

    public double density(Double[] d) {
        Number[] infVal = getInf().value();
        Double[] alpha = (Double[]) infVal;
        if (alpha.length != d.length) {
            throw new IllegalArgumentException("Dimensions don't match");
        }

        double sumAlpha = 0.0;
        for (Double a: alpha){
            a *= omega.value().doubleValue(); // scale
            sumAlpha += a;
        }

        double sumD = 0;
        for (double a : d){
            sumD += a;
        }

        // calc gamma(sumAlpha)
        double gammaSumAlpha = org.apache.commons.math3.special.Gamma.gamma(sumAlpha);
        // calc ∏ Gamma(alpha_i)
        double gammaAlphaProd = 1.0;
        for (Number a : alpha) {
            gammaAlphaProd *= org.apache.commons.math3.special.Gamma.gamma(a.doubleValue());
        }

        //∏ d_i^(alpha_i - 1)
        double dProd = 1.0;
        for (int i = 0; i < d.length; i++) {
            double x = d[i];
            Number a = alpha[i];

            if (x <= 1e-15) {
                return 0.0;
            }

            dProd *= Math.pow(x / sumD, a.doubleValue() - 1);
        }

        // check scaling factor
        double sFactor = Math.pow(sumD, -(sumAlpha - d.length));

        double density = (gammaSumAlpha / gammaAlphaProd) * dProd * sFactor;
        return density;
    }

    @Override
    public Map<String,Value> getParams() {
        return new TreeMap<>() {{
            if (inf!= null) put(infParamName, inf);
            if (omega != null) put(omegaParamName, omega);
            if (mean != null) put(meanParamName, mean);
        }};
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(infParamName)) inf = value;
        else if (paramName.equals(omegaParamName)) omega = value;
        else if (paramName.equals(meanParamName)) mean = value;
        else throw new RuntimeException("Unrecognised parameter name: " + paramName);

        super.setParam(paramName, value); // constructDistribution
    }

    public Value<Number[]> getInf() {
        return inf;
    }

    public Value<Number> getOmega() {
        return omega;
    }

    public Value<Number> getMean() {
        return mean;
    }
}