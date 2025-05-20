package abyss;

import lphy.base.distribution.ParametricDistribution;
import lphy.base.math.MathUtils;
import lphy.core.logger.LoggerUtils;
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

    private Value<Number[]> conc;
    public static final String concParamName = "conc";
    private Value<Number> omega;
    public static final String scaleParamName = "scale";

    public InformedDirichlet(
            @ParameterInfo(name= concParamName, narrativeName = "conc",
                    description="input rates/freqs informing concentration.") Value<Number[]> conc,
            @ParameterInfo(name= scaleParamName, narrativeName = "scale",
                    description = "scale factor") Value<Number> omega) {
        super();
        this.conc = conc;
        this.omega = omega;

    }

    @Override
    protected void constructDistribution(RandomGenerator random) {  }

    @GeneratorInfo(name="InformedDirichlet", verbClause = "have", narrativeName = "Rate- or frequency-informed Dirichlet distribution prior",
            category = GeneratorCategory.PRIOR,
            description="The informed dirichlet probability distribution.")
    public RandomVariable<Double[]> sample() {
        Double[] normConc = normaliseConcentration(conc);
        Double[] dirichlet = new Double[normConc.length];
        double sum = 0.0;

        for (int i = 0; i < normConc.length; i++) {
            dirichlet[i] = MathUtils.randomGamma(normConc[i], 1.0, random);
            sum += dirichlet[i];
        }

        // dirichlet normalised to sum to 1
        for (int i = 0; i < dirichlet.length; i++) {
            dirichlet[i] /= sum;
        }

        int x = 0;
        for (int i = 0; i < dirichlet.length; i++) {
            if (dirichlet[i] < 1e-6) x++;
            if (x > dirichlet.length/10) throw new IllegalArgumentException("At least one dirichlet close to zero, increase " +
                        scaleParamName +" ("+ omega +").");

        }

        return new RandomVariable<>("x", dirichlet, this);
    }

    public double density(Double[] d) {
        Double[] alpha = normaliseConcentration(conc);
        if (alpha.length != d.length) {
            throw new IllegalArgumentException("Dimensions don't match");
        }

        double sumAlpha = 0.0;
        for (Double a: alpha){
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
        for (Double a : alpha) {
            gammaAlphaProd *= org.apache.commons.math3.special.Gamma.gamma(a);
        }

        //∏ d_i^(alpha_i - 1)
        double dProd = 1.0;
        for (int i = 0; i < d.length; i++) {
            double x = d[i];
            Double a = alpha[i];

            if (x <= 1e-15) {
                return 0.0;
            }

            dProd *= Math.pow(x / sumD, a - 1);
        }

        // check scaling factor
        double sFactor = Math.pow(sumD, -(sumAlpha - d.length));

        double density = (gammaSumAlpha / gammaAlphaProd) * dProd * sFactor;
        return density;
    }

    @Override
    public Map<String,Value> getParams() {
        return new TreeMap<>() {{
            if (conc != null) put(concParamName, conc);
            if (omega != null) put(scaleParamName, omega);
        }};
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(concParamName)) conc = value;
        else if (paramName.equals(scaleParamName)) omega = value;
        else throw new RuntimeException("Unrecognised parameter name: " + paramName);

        super.setParam(paramName, value); // constructDistribution
    }

    public Value<Number[]> getConc() {
        return conc;
    }

    public Value<Number> getScale() {
        return omega;
    }

    protected Double[] normaliseConcentration(Value<Number[]> conc) {
        Double[] normConc = (Double[]) this.conc.value();
        double sum = 0.0;
        for (int i = 0; i < normConc.length; i++) {
            sum += normConc[i];
        }

        double omega = getScale().value().doubleValue();
        double finalSum = 0.0;

        // conc normalised and scaled to sum to omega
        for (int i = 0; i < normConc.length; i++) {
            normConc[i] /= sum;
            normConc[i] *= omega * normConc.length;
            finalSum += normConc[i];
        }
        if (Math.abs(finalSum - omega*normConc.length) > 1e-6) throw new IllegalArgumentException("Conc must be normalised to sum to "+ scaleParamName + "*numRates = "+omega*normConc.length+".");
        int x = 0;
        for (int i = 0; i < normConc.length; i++) {
            if (normConc[i] < 1.0) x++;
            if (x > normConc.length/10) LoggerUtils.log.warning("At least one concentration is less than 1 after scaling. Consider increasing "+ scaleParamName + " by "+(double)Math.round(1/normConc[i] * 100d) / 100d+ " times to avoid zeroed rates.");
        }
        return normConc;
    }
}