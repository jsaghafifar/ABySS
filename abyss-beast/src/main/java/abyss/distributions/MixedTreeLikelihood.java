package abyss.distributions;


import beast.base.core.*;
import beast.base.evolution.alignment.Alignment;
import beast.base.spec.evolution.likelihood.TreeLikelihood;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.inference.parameter.RealVectorParam;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * @author Jasmine Saghafifar
 */

@Description("Takes a collection of tree likelihoods with different substitution models.")
public class MixedTreeLikelihood extends Distribution {
    final public Input<List<Distribution>> pLikelihoods =
            new Input<>("likelihood", "individual tree likelihoods", new ArrayList<>());
    final public Input<String> modeInput = new Input<>("mode",
            "how likelihoods should be handled. " +
                    "Site mixture (mix), model averaging (avg), or (both).", Input.Validate.REQUIRED);
    final public Input<RealVectorParam<NonNegativeReal>> metaWeightsInput =
            new Input<>("metaWeights",
                    "tree likelihood weights. Should be an indicator " +
                            "if averaging, or uniform if mixing. Default uniform.", Input.Validate.OPTIONAL);
    final public Input<RealVectorParam<NonNegativeReal>> siteModelWeightsInput = new Input<>("siteModelWeights",
            "estimated site likelihood weights. Required for site mixture (one for each model)", Input.Validate.OPTIONAL);

    protected String mode;
    protected RealVectorParam<NonNegativeReal> metaWeights;
    protected RealVectorParam<NonNegativeReal> siteModelWeights;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        if (pLikelihoods.get().isEmpty()) logP = 0;
        if (pLikelihoods.get().size() > 1 &&
                pLikelihoods.get().getFirst().getInput("data").get() instanceof Alignment alignment) {
            for (int i = 1; i < pLikelihoods.get().size(); i++) {
                if (alignment != pLikelihoods.get().get(i).getInput("data").get())
                    throw new IllegalArgumentException("Data should be same between likelihoods for MixedTreeLikelihood");

            }
        }

        mode = modeInput.get();
        if (metaWeightsInput.get() != null) {
            if (mode.equalsIgnoreCase("avg"))
                throw new IllegalArgumentException("Model averaging does not support custom meta weights");
            metaWeights = metaWeightsInput.get();
            if (mode.equalsIgnoreCase("mix") && metaWeights.size() != pLikelihoods.get().size())
                throw new IllegalArgumentException(
                        "Site mixture model must have same number of meta weights as given tree likelihoods.");
            if (mode.equalsIgnoreCase("both") && metaWeights.size() != (pLikelihoods.get().size() + 1))
                throw new IllegalArgumentException("Site mixture with model averaging must have "+
                        (pLikelihoods.get().size() + 1)+" dimensions (site mixture model weight last).");
            double sum = 0;
            for (int i = 0; i < metaWeights.size(); i++) {
                sum += metaWeights.get(i);
            }
            if (Math.abs(sum - 1) > 1e-6) throw new IllegalArgumentException("Meta weights must sum to 1.");
        } else {
            double[] weights = mode.equalsIgnoreCase("both") ?
                    new double[pLikelihoods.get().size() + 1] : new double[pLikelihoods.get().size()];
            Arrays.fill(weights, 1.0 / weights.length);
            metaWeights = new RealVectorParam<>(weights, NonNegativeReal.INSTANCE);
        }

        if (siteModelWeightsInput.get() != null) {
            siteModelWeights = siteModelWeightsInput.get();
            if (siteModelWeights.size() != pLikelihoods.get().size())
                throw new IllegalArgumentException("Site weights must have same dimensions as given tree likelihoods.");
            double sum = 0;
            for (int i = 0; i < siteModelWeights.size(); i++) {
                sum += siteModelWeights.get(i);
            }
            if (Math.abs(sum - 1) > 1e-6) throw new IllegalArgumentException("Site weights must sum to 1.");
        } else if (!mode.equalsIgnoreCase("avg")) {
            throw new IllegalArgumentException("Site weights must be specified for site mixture model.");
        }
    }


    @Override
    public double calculateLogP() {
        logP = 0;
        double[] p = getPartialLogLikelihoods();

        for (double v : p) {
            if (Double.isInfinite(v) || Double.isNaN(v)) {
                logP += v;
                return logP;
            }
        }

        if (mode.equalsIgnoreCase("mix")) logP = p[0];
        else logP = logSumExp(p);

        return logP;
    }

    @Override
    public void sample(State state, Random random) {    }

    @Override
    public List<String> getArguments() {
        List<String> arguments = new ArrayList<>();
        if (siteModelWeightsInput.get() != null) arguments.add(siteModelWeightsInput.get().getID());
        for (Distribution likelihood : pLikelihoods.get()) {
            arguments.addAll(likelihood.getArguments());
        }
        return arguments;
    }

    @Override
    public List<String> getConditions() {
        List<String> conditions = new ArrayList<>();
        if (metaWeightsInput.get() != null) conditions.add(metaWeightsInput.get().getID());
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

    public double[] getPartialLogLikelihoods() {
        double[] p;
        int nrOfLikelihoods = pLikelihoods.get().size();

        if (mode.equalsIgnoreCase("mix")) {
            p = new double[1];
            for (int i = 0; i < nrOfLikelihoods; i++) {
                Distribution likelihood = pLikelihoods.get().get(i);
                if (likelihood.somethingIsDirty()) likelihood.calculateLogP();
            }
            p[0] = calculateSiteMixtureLogP();

        } else if (mode.equalsIgnoreCase("both")) {
            p = new double[nrOfLikelihoods + 1];
            for (int i = 0; i < nrOfLikelihoods; i++) {
                Distribution likelihood = pLikelihoods.get().get(i);
                p[i] = likelihood.somethingIsDirty() ? likelihood.calculateLogP() : likelihood.getCurrentLogP();
                p[i] += Math.log(metaWeights.get(i));
            }

            p[nrOfLikelihoods] = calculateSiteMixtureLogP();
            p[nrOfLikelihoods] += Math.log(metaWeights.get(nrOfLikelihoods));

        } else if (mode.equalsIgnoreCase("avg")) {
            p = new double[nrOfLikelihoods];
            for (int i = 0; i < nrOfLikelihoods; i++) {
                Distribution likelihood = pLikelihoods.get().get(i);
                p[i] = likelihood.somethingIsDirty() ? likelihood.calculateLogP() : likelihood.getCurrentLogP();
                p[i] += Math.log(metaWeights.get(i));
            }

        } else throw new UnsupportedOperationException();

        return p;
    }

    public double[] getMetaWeights() {
        return metaWeights.getValues();
    }

    private double calculateSiteMixtureLogP() {
        int nrOfLikelihoods = pLikelihoods.get().size();
        int nrOfPatterns = ((TreeLikelihood) pLikelihoods.get().getFirst()).getPatternLogLikelihoods().length;
        double[] pSite = new double[nrOfLikelihoods];
        double logPMixture = 0;
        Alignment data = ((TreeLikelihood) pLikelihoods.get().getFirst()).dataInput.get();

        for (int j = 0; j < nrOfPatterns; j++) {
            for (int i = 0; i < nrOfLikelihoods; i++) {
                Distribution likelihood = pLikelihoods.get().get(i);
                double[] patternLogLikelihoods = ((TreeLikelihood) likelihood).getPatternLogLikelihoods();
                pSite[i] = patternLogLikelihoods[j];
                pSite[i] += Math.log(siteModelWeights.get(i));
            }
            logPMixture += logSumExp(pSite) * data.getPatternWeight(j);
        }

        return logPMixture;
    }

    private double logSumExp(double[] p) {
        double logP = 0;
        double max = Arrays.stream(p).max().getAsDouble();
        for (double v : p) {
            logP += FastMath.exp(v - max);
        }
        logP = max + Math.log(logP);

        return logP;
    }
    
} // class MixedTreeLikelihood
