package abyss.distributions;


import beast.base.core.*;
import beast.base.evolution.alignment.Alignment;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;


@Description("Takes a collection of tree likelihoods with different substitution models.")
public class MixedTreeLikelihood extends Distribution {
    final public Input<List<Distribution>> pLikelihoods =
            new Input<>("likelihood", "individual tree likelihoods", new ArrayList<>());
    final public Input<RealParameter> weightVectorInput =
            new Input<>("weights", "likelihood weights. Equal if not specified", Input.Validate.OPTIONAL);

    protected RealParameter weightVector;

    @Override
    public void initAndValidate() {
        // System.setProperty("java.only", "true");
        super.initAndValidate();

        if (pLikelihoods.get().isEmpty()) logP = 0;

        if (pLikelihoods.get().size() > 1 &&
                pLikelihoods.get().get(0).getInput("data").get() instanceof Alignment) {
            Alignment alignment = (Alignment) pLikelihoods.get().get(0).getInput("data").get();
            for (int i = 1; i < pLikelihoods.get().size(); i++) {
                if (alignment != pLikelihoods.get().get(i).getInput("data").get())
                    throw new IllegalArgumentException("Data should be same between likelihoods for MixedTreeLikelihood");

            }
        }

        if (weightVectorInput.get() != null) {
            weightVector = weightVectorInput.get();
            if (weightVector.getDimension() != pLikelihoods.get().size())
                throw new IllegalArgumentException("Must be same number of likelihood weights as subst models.");
            double sum = 0;
            for (int i = 0; i < weightVector.getDimension(); i++) {
                sum += weightVector.getArrayValue(i);
            }
            if (Math.abs(sum - 1) > 1e-6) throw new IllegalArgumentException("Likelihood weights must sum to 1.");
        } else {
            Double[] weights = new Double[pLikelihoods.get().size()];
            Arrays.fill(weights, 1.0 / weights.length);
            weightVector = new RealParameter(weights);
        }
    }


    @Override
    public double calculateLogP() {
        logP = 0;
        double[] p = new double[pLikelihoods.get().size()];

        for (int i = 0; i < p.length; i++) {
            Distribution likelihood = pLikelihoods.get().get(i);
            p[i] = likelihood.isDirtyCalculation() ? likelihood.calculateLogP() : likelihood.getCurrentLogP();
            p[i] += Math.log(weightVector.getArrayValue(i));
            if (Double.isInfinite(p[i]) || Double.isNaN(p[i])) {
                logP += p[i];
                return logP;
            }
        }

        double max = Arrays.stream(p).max().getAsDouble();

        for (int i = 0; i < p.length; i++) {
            logP += FastMath.exp(p[i]-max);
        }
        logP = max + Math.log(logP);

        return logP;
    }

    @Override
    public void sample(State state, Random random) {
        if (sampledFlag)
            return;

        sampledFlag = true;

        for (Distribution likelihood : pLikelihoods.get()) {
            likelihood.sample(state, random);
        }
    }

    @Override
    public List<String> getArguments() {
        List<String> arguments = new ArrayList<>();
        for (Distribution likelihood : pLikelihoods.get()) {
            arguments.addAll(likelihood.getArguments());
        }
        return arguments;
    }

    @Override
    public List<String> getConditions() {
        List<String> conditions = new ArrayList<>();
        for (Distribution likelihood : pLikelihoods.get()) {
            conditions.addAll(likelihood.getConditions());
        }
        conditions.removeAll(getArguments());

        return conditions;
    }

    @Override
    public boolean isStochastic() {
        for (Distribution likelihood : pLikelihoods.get()) {
            if (likelihood.isStochastic())
                return true;
        }
        
        return false;
    }
    
    @Override
    public double getNonStochasticLogP() {
        double logP = 0;
        double[] p = new double[pLikelihoods.get().size()];

        for (int i = 0; i < p.length; i++) {
            Distribution likelihood = pLikelihoods.get().get(i);
            p[i] = likelihood.getNonStochasticLogP() + Math.log(weightVector.getArrayValue(i));
            if (Double.isInfinite(p[i]) || Double.isNaN(p[i])) {
                logP += p[i];
                return logP;
            }
        }

        double max = Arrays.stream(p).max().getAsDouble();

        for (int i = 0; i < p.length; i++) {
            logP += FastMath.exp(p[i]-max);
        }
        logP = max + Math.log(logP);

        return logP;
    }
    
} // class MixedTreeLikelihood
