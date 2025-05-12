package abyss.distributions;

import abyss.inference.AbyssSVS;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;

import java.util.List;
import java.util.Random;


@Description("Bernoulli distribution, used as prior for SVS rate indicator booleans in rate matrices.")
public class SVSPrior extends Distribution {

    final public Input<RealParameter> pInput = new Input<>("p", "probability p parameter.", Input.Validate.REQUIRED);
    final public Input<BooleanParameter> trialsInput = new Input<>("parameter", "the results of a series of bernoulli trials.");

//    final public Input<IntegerParameter> minSuccessesInput = new Input<>("minSuccesses",
//            "Optional condition: the minimum number of ones in the boolean array.");
    protected SubstitutionModel substitutionModel;
    private int numStates = substitutionModel.getStateCount();

    public double calculateLogP() {
        logP = 0.0;

        BooleanParameter trials = trialsInput.get();
        double p = pInput.get().getValue();
//        IntegerParameter minSuccesses = minSuccessesInput.get();


        if (p < 0 | p > 1)
            return Double.NEGATIVE_INFINITY;

//        for (int i = 0; i < trials.getDimension(); i++) {
//            logP += Math.log(trials.getValue(i) ? prob : 1.0 - prob);
//        }

        // uniform false not needed once impl

        int n = trials.getDimension();
        int sum = 0;
        for (int i = 0; i < n; i ++) {
            if (trials.getValue(i)) {
                sum ++;
            }
        }

        // reject if < minSuccesses
        if (sum < numStates)
            return Double.NEGATIVE_INFINITY;

        Boolean[] trialValues = trials.getValues();
        Double[] indicatorValues = new Double[n];
        for (int i = 0; i < n; i++) {
            indicatorValues[i] = trialValues[i] ? 1.0 : 0.0;
        }
        if (!AbyssSVS.Utils.connectedAndWellConditioned(null, substitutionModel) ||
                !AbyssSVS.Utils.isStronglyConnected(indicatorValues, numStates, false))
            return Double.NEGATIVE_INFINITY;
        // Boolean[] trials != Double[] indicatorVals

        // Binomial distribution
        double logFactorialn=0, logFactorial1=0, logFactorial2=0;
        for (int i = 2; i <= n;  i ++) logFactorialn += Math.log(i);
        for (int i = 2; i <= sum; i ++) logFactorial1 += Math.log(i);
        for (int i = 2; i <= (n-sum); i ++) logFactorial2 += Math.log(i);

        logP = logFactorialn - (logFactorial1 + logFactorial2) + sum*Math.log(p) + (n-sum)*Math.log(1-p);

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

//    @Override
//    public void initAndValidate() {
//        if (pInput.get().getDimension() != 1 && pInput.get().getDimension() != trialsInput.get().getDimension()) {
//            throw new RuntimeException("p parameter must be size 1 or the same size as trials parameter but it was dimension " + pInput.get().getDimension());
//        }
//    }
}
