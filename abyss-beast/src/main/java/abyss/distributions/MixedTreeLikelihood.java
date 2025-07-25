package abyss.distributions;


import beast.base.core.*;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
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
    final public Input<RealParameter> metaWeightsInput =
            new Input<>("metaWeights",
                    "tree likelihood weights. Should be an indicator " +
                            "if averaging, or uniform if mixing. Default uniform.", Input.Validate.OPTIONAL);
    final public Input<RealParameter> siteModelWeightsInput = new Input<>("siteModelWeights",
            "estimated site likelihood weights. Required for site mixture (one for each model)", Input.Validate.OPTIONAL);

    protected String mode;
    protected RealParameter metaWeights;
    protected RealParameter siteModelWeights;

    @Override
    public void initAndValidate() {
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

        mode = modeInput.get();
        if (metaWeightsInput.get() != null) {
            if (mode.equalsIgnoreCase("avg"))
                throw new IllegalArgumentException("Model averaging does not support custom meta weights");
            metaWeights = metaWeightsInput.get();
            if (mode.equalsIgnoreCase("mix") && metaWeights.getDimension() != pLikelihoods.get().size())
                throw new IllegalArgumentException(
                        "Site mixture model must have same number of meta weights as given tree likelihoods.");
            if (mode.equalsIgnoreCase("both") && metaWeights.getDimension() != (pLikelihoods.get().size() + 1))
                throw new IllegalArgumentException("Site mixture with model averaging must have "+
                        (pLikelihoods.get().size() + 1)+" dimensions (site mixture model weight last).");
            double sum = 0;
            for (int i = 0; i < metaWeights.getDimension(); i++) {
                sum += metaWeights.getArrayValue(i);
            }
            if (Math.abs(sum - 1) > 1e-6) throw new IllegalArgumentException("Meta weights must sum to 1.");
        } else {
            Double[] weights = mode.equalsIgnoreCase("both") ?
                    new Double[pLikelihoods.get().size() + 1] : new Double[pLikelihoods.get().size()];
            Arrays.fill(weights, 1.0 / weights.length);
            metaWeights = new RealParameter(weights);
        }

        if (siteModelWeightsInput.get() != null) {
            siteModelWeights = siteModelWeightsInput.get();
            if (siteModelWeights.getDimension() != pLikelihoods.get().size())
                throw new IllegalArgumentException("Site weights must have same dimensions as given tree likelihoods.");
            double sum = 0;
            for (int i = 0; i < metaWeights.getDimension(); i++) {
                sum += metaWeights.getArrayValue(i);
            }
            if (Math.abs(sum - 1) > 1e-6) throw new IllegalArgumentException("Meta weights must sum to 1.");
        } else if (!mode.equalsIgnoreCase("avg")) {
            throw new IllegalArgumentException("Site weights must be specified for site mixture model.");
        }
    }


    @Override
    public double calculateLogP() {
        logP = 0;
        double[] p = getPartialLogLikelihoods();

        for (int i = 0; i < p.length; i++) {
            if (Double.isInfinite(p[i]) || Double.isNaN(p[i])) {
                logP += p[i];
                return logP;
            }
        }

        if (mode.equalsIgnoreCase("mix")) logP = p[0];
        else logP = logSumExp(p);

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
        if (siteModelWeightsInput.get() != null) arguments.add(siteModelWeightsInput.get().getID());
        for (Distribution likelihood : pLikelihoods.get()) {
            arguments.addAll(likelihood.getArguments());
        }
        return arguments;
    }

    @Override
    public List<String> getConditions() {
        List<String> conditions = new ArrayList<>();
        conditions.add(modeInput.get());
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
            p[0] = calculateSiteMixtureLogP();

        } else if (mode.equalsIgnoreCase("both")) {
            p = new double[nrOfLikelihoods + 1];
            for (int i = 0; i < nrOfLikelihoods; i++) {
                Distribution likelihood = pLikelihoods.get().get(i);
                p[i] = likelihood.isDirtyCalculation() ? likelihood.calculateLogP() : likelihood.getCurrentLogP();
                p[i] += Math.log(metaWeights.getArrayValue(i));
            }

            p[nrOfLikelihoods] = calculateSiteMixtureLogP();
            p[nrOfLikelihoods] += Math.log(metaWeights.getArrayValue(nrOfLikelihoods));

        } else if (mode.equalsIgnoreCase("avg")) {
            p = new double[nrOfLikelihoods];
            for (int i = 0; i < nrOfLikelihoods; i++) {
                Distribution likelihood = pLikelihoods.get().get(i);
                p[i] = likelihood.isDirtyCalculation() ? likelihood.calculateLogP() : likelihood.getCurrentLogP();
                p[i] += Math.log(metaWeights.getArrayValue(i));
            }

        } else throw new UnsupportedOperationException();

        return p;
    }

    private double calculateSiteMixtureLogP() {
        int nrOfLikelihoods = pLikelihoods.get().size();
        int nrOfPatterns = ((TreeLikelihood) pLikelihoods.get().get(0)).getPatternLogLikelihoods().length;
        double[] pSite = new double[nrOfLikelihoods];
        double logPMixture = 0;
        Alignment data = ((TreeLikelihood) pLikelihoods.get().get(0)).dataInput.get();

        for (int j = 0; j < nrOfPatterns; j++) {
            for (int i = 0; i < nrOfLikelihoods; i++) {
                Distribution likelihood = pLikelihoods.get().get(i);
                double[] patternLogLikelihoods = ((TreeLikelihood) likelihood).getPatternLogLikelihoods();
                pSite[i] = patternLogLikelihoods[j];
                pSite[i] += Math.log(siteModelWeights.getArrayValue(i));
            }
            logPMixture += logSumExp(pSite) * data.getPatternWeight(j);
        }

        return logPMixture;
    }

    private double logSumExp(double[] p) {
        double logP = 0;
        double max = Arrays.stream(p).max().getAsDouble();
        for (int i = 0; i < p.length; i++) {
            logP += FastMath.exp(p[i] - max);
        }
        logP = max + Math.log(logP);

        return logP;
    }
    
} // class MixedTreeLikelihood
