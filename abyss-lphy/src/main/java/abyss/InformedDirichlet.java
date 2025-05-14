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

    private Value<Number[]> in;
    public static final String rfInputsParamName = "rfInputs";
    private Value<Number> omega;
    public static final String omegaParamName = "omega";

    public InformedDirichlet(
            @ParameterInfo(name= rfInputsParamName, narrativeName = "rfInputs",
                    description="input rates/freqs informing concentration.") Value<Number[]> in,
            @ParameterInfo(name=omegaParamName, narrativeName = "omega",
                    description = "scale factor") Value<Number> omega) {
        super();
        this.in = in;
        this.omega = omega;

    }

    @Override
    protected void constructDistribution(RandomGenerator random) {  }

    @GeneratorInfo(name="InformedDirichlet", verbClause = "have", narrativeName = "Rate- or frequency-informed Dirichlet distribution prior",
            category = GeneratorCategory.PRIOR,
            description="The informed dirichlet probability distribution.")
    public RandomVariable<Double[]> sample() {
        Double[] conc = normaliseConcentration();
        Double[] dirichlet = new Double[conc.length];
        double sum = 0.0;

        for (int i = 0; i < conc.length; i++) {
            dirichlet[i] = MathUtils.randomGamma(conc[i], 1.0, random);
            sum += dirichlet[i];
        }

        // dirichlet normalised to sum to 1
        for (int i = 0; i < dirichlet.length; i++) {
            dirichlet[i] /= sum;
        }

        return new RandomVariable<>("x", dirichlet, this);
    }

    public double density(Double[] d) {
        Double[] alpha = normaliseConcentration();
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
            if (in != null) put(rfInputsParamName, in);
            if (omega != null) put(omegaParamName, omega);
        }};
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(rfInputsParamName)) in = value;
        else if (paramName.equals(omegaParamName)) omega = value;
        else throw new RuntimeException("Unrecognised parameter name: " + paramName);

        super.setParam(paramName, value); // constructDistribution
    }

    public Value<Number[]> getRFInputs() {
        return in;
    }

    public Value<Number> getOmega() {
        return omega;
    }

    protected Double[] normaliseConcentration() {
        Double[] conc = (Double[]) in.value();
        double sum = 0.0;
        for (int i = 0; i < conc.length; i++) {
            sum += conc[i];
        }

        double omega = getOmega().value().doubleValue();

        // conc normalised and scaled to sum to omega
        for (int i = 0; i < conc.length; i++) {
            conc[i] /= sum;
            conc[i] *= omega;
        }
        return conc;
    }
}