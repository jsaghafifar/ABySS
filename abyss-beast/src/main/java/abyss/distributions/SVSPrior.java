package abyss.distributions;

import abyss.inference.AbyssSVS;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.parameter.*;
import beast.base.spec.type.Simplex;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * @author Jasmine Saghafifar
 */

@Description("Bernoulli distribution, used as prior for SVS rate indicator booleans in rate matrices.")
public class SVSPrior extends Distribution {

    final public Input<RealScalarParam<PositiveReal>> scaleInput = new Input<>("scale", "probability scale parameter.", Input.Validate.REQUIRED);
    final public Input<RealScalarParam<PositiveReal>> sensitivityInput = new Input<>("sens", "probability sensitivity parameter.", Input.Validate.OPTIONAL);
    final public Input<BoolVectorParam> indicatorsInput = new Input<>("indicators", "the results of a series of bernoulli trials.", Input.Validate.REQUIRED);
    final public Input<Integer> nrOfStatesInput = new Input<>("nrOfStates", "number of states being used in Q matrix.", Input.Validate.REQUIRED);
    final public Input<SimplexParam> empiricalQInput = new Input<>("empiricalQ", "empirical Q matrix informing bernoulli trials.", Input.Validate.OPTIONAL);
    final public Input<Boolean> isSymmetricInput = new Input<>("symmetric", "whether Q is reversible/symmetric rates.", Input.Validate.REQUIRED);

    public double calculateLogP() {
        this.logP = 0.0;
        double a = this.scaleInput.get().get();
        double k;
        if (this.sensitivityInput.get() != null)
            k = this.sensitivityInput.get().get();
        else k = 1.0;

        Boolean[] indicators = this.indicatorsInput.get().getElements().toArray(new Boolean[0]);
        int nrOfStates = this.nrOfStatesInput.get();
        int nrOfRates = isSymmetricInput.get() ? (nrOfStates*nrOfStates-nrOfStates)/2 : nrOfStates*nrOfStates-nrOfStates;
        Double[] rates = new Double[nrOfRates];
        if (this.empiricalQInput.get() != null) {
            if (this.empiricalQInput.get().size() != nrOfRates) throw new IllegalArgumentException("Number of empirical rates must be "+nrOfRates+".");
            rates = this.empiricalQInput.get().getElements().toArray(new Double[0]);
        }
        else {
            Arrays.fill(rates, (double) 1 /nrOfRates); // uninformed Bernoulli: each indicator will have equal probability of being active
        }

        double[] p = new double[rates.length];
        for (int i = 0; i < p.length; i++) {
            double x = Math.log(p.length * rates[i]);
            x *= -k; // sensitivity to empirical Q input is applied
            p[i] = a/(a + FastMath.exp(x)); // probability of indicator scaled
            if (p[i] < 0 || p[i] > 1)
                return Double.NEGATIVE_INFINITY;
        }
        for (int i = 0; i < p.length; i++) {
            logP += indicators[i] ? Math.log(p[i]) : Math.log(1-p[i]);
        }

        int sum = 0;
        for (Boolean indicator : indicators) {
            if (indicator) {
                sum++;
            }
        }
        if (sum < nrOfStates)
            return Double.NEGATIVE_INFINITY;

        Double[] indicatorValues = new Double[indicators.length];
        for (int i = 0; i < indicators.length; i++) {
            indicatorValues[i] = indicators[i] ? 1.0 : 0.0;
        }
        Boolean sym = this.isSymmetricInput.get();
        if (!AbyssSVS.Utils.connectedAndWellConditioned(p) ||
                !AbyssSVS.Utils.isStronglyConnected(indicatorValues, nrOfStates, sym))
            return Double.NEGATIVE_INFINITY;

        return logP;
    }

    @Override
    public List<String> getArguments() {
        List<String> args = new ArrayList<>();
        args.add(indicatorsInput.get().getID());
        return args;
    }

    @Override
    public List<String> getConditions() {
        List<String> conds = new ArrayList<>();
        conds.add(scaleInput.get().getID());
        conds.add(sensitivityInput.get().getID());
        return conds;
    }

    @Override
    public void sample(State state, Random random) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void initAndValidate() {
        int nrOfStates = nrOfStatesInput.get();
        int symQ = isSymmetricInput.get() ? 2 : 1;
        if (indicatorsInput.get().size() != (nrOfStates * nrOfStates - nrOfStates)/symQ) {
            throw new RuntimeException("Indicators must have " +
                    (nrOfStates * nrOfStates - nrOfStates)/symQ + " dimension, based on specified number of states.");
        }

        if (this.sensitivityInput.get() != null)
            if (this.empiricalQInput.get() == null) throw new RuntimeException("SVS probability sensitivity should" +
                    " not be specified if it is not informed by empirical Q.");
    }
}
