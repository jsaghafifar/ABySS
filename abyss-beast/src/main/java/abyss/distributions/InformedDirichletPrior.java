package abyss.distributions;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.distribution.TensorDistribution;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.type.RealVector;
import beast.base.spec.type.Simplex;
import org.apache.commons.numbers.gamma.LogGamma;
import org.apache.commons.statistics.distribution.GammaDistribution;

import java.util.Arrays;
import java.util.List;

/**
 * @author Jasmine Saghafifar
 */

@Description("Informed Dirichlet distribution. Implementation not yet verified.")
public class InformedDirichletPrior extends TensorDistribution<Simplex, Double> {
    final public Input<RealVector<PositiveReal>> alphaInput = new Input<>("alpha", "coefficients of the Dirichlet distribution", Validate.REQUIRED);
    final public Input<RealScalarParam<PositiveReal>> scaleInput = new Input<>("scale", "scale of coefficients", Validate.REQUIRED);
    private GammaDistribution[] gammas;

    public InformedDirichletPrior() {}

    public InformedDirichletPrior(Simplex param, RealVector<PositiveReal> alpha, PositiveReal scale) {
        try {
            initByName("param", param, "alpha", alpha, "scale", scale);
        } catch (Exception e) {
            throw new RuntimeException("Failed to initialize " + getClass().getSimpleName() +
                    " via initByName in constructor.", e );
        }
    }

    @Override
    public void initAndValidate() {
        refresh();
        super.initAndValidate(); // load param here
        if (alphaInput.get().size() != dimension())
            throw new IllegalArgumentException("Dimensions of alpha and param should be the same, " +
                    "but dim(alpha)=" + alphaInput.get().size() + " and dim(x)=" + dimension());

    }

    @Override
    public void refresh() {
        List<Double> alpha = alphaInput.get().getElements();
        double scale = scaleInput.get().get();
        gammas = new GammaDistribution[alpha.size()];

        for (int i = 0; i < alpha.size(); i++) {
            if (gammas[i] == null)
                gammas[i] = GammaDistribution.of(alpha.get(i), scale);

            // Floating point comparison
            if (Math.abs(gammas[i].getShape() - alpha.get(i)) > EPS)
                gammas[i] = GammaDistribution.of(alpha.get(i), scale);
        }
    }

    @Override
    protected double calcLogP(Double... value) {
        return this.calcLogP(Arrays.asList(value));
    }

    @Override
    public double calculateLogP() {
        // Avoid unnecessary conversions, use List<> directly for better performance
        logP = this.calcLogP(param.getElements());
        return logP;
    }

    private double calcLogP(List<Double> value) {
        refresh(); // this make sure distribution parameters are updated if they are sampled during MCMC

        List<Double> alpha = alphaInput.get().getElements();
        if (alpha.size() != value.size())
            throw new IllegalArgumentException("Dimensions of alpha and param should be the same, " +
                    "but dim(alpha)=" + alpha.size() + " and dim(x)=" + value.size());
        double scale = scaleInput.get().get();

        double logP = 0;
        for (int i = 0; i < value.size(); i++) {
            logP += (alpha.get(i)*scale - 1) * Math.log(value.get(i)); // a-scale or a*scale-1? check
            logP -= LogGamma.value(alpha.get(i)*scale); // check
        }
        double alphaSum = alpha.stream().mapToDouble(Double::doubleValue).sum();
        alphaSum += scale*alpha.size(); // check
        logP += LogGamma.value(alphaSum);

        // area = sumX^(dim-1)
        double sumX = value.stream()              // Stream<Double>
                .mapToDouble(Double::doubleValue) // unbox to double
                .sum();
        final double expectedSum = 1.0;
        if (Math.abs(sumX - expectedSum) > 1e-6) {
            Log.trace("sum of values (" + sumX +") differs significantly from the expected sum of values (" + expectedSum +")");
            return Double.NEGATIVE_INFINITY;
        }
        logP -= (value.size() - 1) * Math.log(sumX);

        return logP;
    }

    @Override
    public List sample() {
        return List.of();
    }

    @Override
    public Double getLowerBoundOfParameter() {
        // all gammas have the same bounds
        return gammas[0].getSupportLowerBound();
    }

    @Override
    public Double getUpperBoundOfParameter() {
        return gammas[0].getSupportUpperBound();
    }

}
