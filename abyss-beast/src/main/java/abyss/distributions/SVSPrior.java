package abyss.distributions;

import abyss.inference.AbyssSVS;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
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

    final public Input<RealParameter> scaleInput = new Input<>("scale", "probability scale parameter.", Input.Validate.REQUIRED);
    final public Input<RealParameter> sensitivityInput = new Input<>("sens", "probability sensitivity parameter.", Input.Validate.OPTIONAL);
    final public Input<BooleanParameter> indicatorsInput = new Input<>("indicators", "the results of a series of bernoulli trials.", Input.Validate.REQUIRED);
    final public Input<Integer> nrOfStatesInput = new Input<>("nrOfStates", "number of states being used in Q matrix.", Input.Validate.REQUIRED);
    final public Input<RealParameter> empiricalQInput = new Input<>("empiricalQ", "empirical Q matrix informing bernoulli trials.", Input.Validate.OPTIONAL);
    final public Input<Boolean> isSymmetricInput = new Input<>("symmetric", "whether Q is reversible/symmetric rates.", Input.Validate.REQUIRED);

    public double calculateLogP() {
        this.logP = 0.0;
        double a = this.scaleInput.get().getValue();
        double k;
        if (this.sensitivityInput.get() != null)
            k = this.sensitivityInput.get().getValue();
        else k = 1.0;

        Boolean[] indicators = this.indicatorsInput.get().getValues();
        int nrOfStates = this.nrOfStatesInput.get();

        Double[][] q = new Double[nrOfStates][nrOfStates];
        for (int i = 0; i < nrOfStates; i++) {
            if (this.empiricalQInput.get() != null) q[i] = this.empiricalQInput.get().getRowValues(i);
            else Arrays.fill(q[i], 1.0); // uninformed Bernoulli: each indicator will have equal probability of being active
        }
        Double[] rates = getQValues(q);

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

    private Double[] getQValues(Double[][] q) {
        Boolean sym = this.isSymmetricInput.get();
        int symQ = sym ? 2 : 1;
        int n = this.nrOfStatesInput.get();

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

    @Override
    public void initAndValidate() {
        int nrOfStates = nrOfStatesInput.get();
        int symQ = isSymmetricInput.get() ? 2 : 1;
        if (indicatorsInput.get().getDimension() != (nrOfStates * nrOfStates - nrOfStates)/symQ) {
            throw new RuntimeException("Indicators must have " +
                    (nrOfStates * nrOfStates - nrOfStates)/symQ + " dimension, based on specified number of states.");
        }

        if (this.sensitivityInput.get() != null)
            if (this.empiricalQInput.get() == null) throw new RuntimeException("SVS probability sensitivity should" +
                    " not be specified if it is not informed by empirical Q.");
    }
}
