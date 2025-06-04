package abyss.distributions;

import abyss.inference.AbyssSVS;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math3.util.FastMath;

import java.util.List;
import java.util.Random;


@Description("Bernoulli distribution, used as prior for SVS rate indicator booleans in rate matrices.")
public class SVSPrior extends Distribution {

    final public Input<RealParameter> scaleInput = new Input<>("scale", "probability scale parameter.", Input.Validate.REQUIRED);
    final public Input<RealParameter> shapeInput = new Input<>("shape", "probability shape parameter.", Input.Validate.REQUIRED);
    final public Input<BooleanParameter> indicatorsInput = new Input<>("indicators", "the results of a series of bernoulli trials.");
    final public Input<RealParameter> empiricalRatesInput = new Input<>("rates", "rates informing bernoulli trials", Input.Validate.REQUIRED);
    final public Input<Boolean> isSymmetricInput = new Input<>("symmetric", "whether rates come from symmetric matrix", Input.Validate.REQUIRED);

    public double calculateLogP() {
        this.logP = 0.0;
        double a = this.scaleInput.get().getValue();
        double k = this.shapeInput.get().getValue();
        Boolean[] indicators = this.indicatorsInput.get().getValues();
        Double[] rates = this.empiricalRatesInput.get().getValues();
        Boolean sym = this.isSymmetricInput.get();

        double[] p = new double[rates.length];
        for (int i = 0; i < p.length; i++) {
            double x = Math.log(p.length * rates[i]);
            x *= -k; // scaled
            p[i] = a/(a + FastMath.exp(x));
            if (p[i] < 0 || p[i] > 1)
                return Double.NEGATIVE_INFINITY;
        }
        for (int i = 0; i < p.length; i++) {
            logP += indicators[i] ? Math.log(p[i]) : Math.log(1-p[i]);
        }

        int nrOfStates;
        if (sym) nrOfStates = (int) Math.abs((-1 - Math.sqrt(1+8*rates.length))/2);
        else nrOfStates = (int) Math.abs((-1 - Math.sqrt(1+4*rates.length))/2);

        int sum = 0;
        for (int i = 0; i < indicators.length; i ++) {
            if (indicators[i]) {
                sum ++;
            }
        }
        if (sum < nrOfStates)
            return Double.NEGATIVE_INFINITY;

        Double[] indicatorValues = new Double[indicators.length];
        for (int i = 0; i < indicators.length; i++) {
            indicatorValues[i] = indicators[i] ? 1.0 : 0.0;
        }
        if (!AbyssSVS.Utils.connectedAndWellConditioned(p) ||
                !AbyssSVS.Utils.isStronglyConnected(indicatorValues, nrOfStates, sym))
            return Double.NEGATIVE_INFINITY;

        return logP;
    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void initAndValidate() {
        if (indicatorsInput.get().getDimension() != empiricalRatesInput.get().getDimension()) {
            throw new RuntimeException("Indicators must be same size as empirical rates parameter but it was dimension " + empiricalRatesInput.get().getDimension());
        }
    }
}
