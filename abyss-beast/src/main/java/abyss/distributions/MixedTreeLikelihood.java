package abyss.distributions;


import abyss.evolution.substitution.ABySSModelAveraging;
import beast.base.core.*;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;


@Description("Takes a collection of tree likelihoods with different substitution models.")
public class MixedTreeLikelihood extends TreeLikelihood {
    final public Input<List<TreeLikelihood>> pLikelihoods =
            new Input<>("likelihood", "individual tree likelihoods", new ArrayList<>());
    final public Input<RealParameter> weightVectorInput =
            new Input<>("weights", "likelihood weights. Equal if not specified", Input.Validate.OPTIONAL);

    protected List<SubstitutionModel> substModels;
    protected RealParameter weightVector;

    @Override
    public void initAndValidate() {
        // System.setProperty("java.only", "true");
        super.initAndValidate();

        if (pLikelihoods.get().isEmpty()) logP = 0;

        if (m_siteModel.substModelInput.get() instanceof ABySSModelAveraging) { // TODO change name to reflect new purpose
            List<GeneralSubstitutionModel> modelList = ((ABySSModelAveraging) m_siteModel.substModelInput.get()).substModelInput.get();
            substModels.addAll(modelList);
        } else {
            throw new IllegalArgumentException("Expected subst model to be of type ABySSModelAveraging");
        }

        if (weightVectorInput.get() != null) {
            weightVector = weightVectorInput.get();
            if (weightVector.getDimension() != substModels.size())
                throw new IllegalArgumentException("Must be same number of likelihood weights as subst models.");
            double sum = 0;
            for (int i = 0; i < weightVector.getDimension(); i++) {
                sum += weightVector.getArrayValue(i);
            }
            if (Math.abs(sum - 1) > 1e-6) throw new IllegalArgumentException("Likelihood weights must sum to 1.");
        } else {
            Double[] weights = new Double[substModels.size()];
            Arrays.fill(weights, 1.0 / weights.length);
            weightVector = new RealParameter(weights);
        }
//        for(Distribution dists : pDistributions.get()) {
//        	logP += dists.calculateLogP();
//        }
    }


    @Override
    public double calculateLogP() {
        logP = 0;

        for (int i = 0; i < pLikelihoods.get().size(); i++) {
            TreeLikelihood likelihood = pLikelihoods.get().get(i);
            double p;
            if (likelihood.isDirtyCalculation()) {
                p = likelihood.calculateLogP();
            } else {
                p = likelihood.getCurrentLogP();
            }
            if (Double.isInfinite(p) || Double.isNaN(p)) {
                logP += p;
                return logP;
            } else {
                logP += p * weightVectorInput.get().getArrayValue(i);
            }
        }

        return logP;
    }

    @Override
    public void sample(State state, Random random) {
        if (sampledFlag)
            return;

        sampledFlag = true;

        for (TreeLikelihood likelihood : pLikelihoods.get()) {
            likelihood.sample(state, random);
        }
    }

    @Override
    public List<String> getArguments() {
        List<String> arguments = new ArrayList<>();
        for (TreeLikelihood likelihood : pLikelihoods.get()) {
            arguments.addAll(likelihood.getArguments());
        }
        return arguments;
    }

    @Override
    public List<String> getConditions() {
        List<String> conditions = new ArrayList<>();
        for (TreeLikelihood likelihood : pLikelihoods.get()) {
            conditions.addAll(likelihood.getConditions());
        }
        conditions.removeAll(getArguments());

        return conditions;
    }

    @Override
    public boolean isStochastic() {
        for (TreeLikelihood likelihood : pLikelihoods.get()) {
            if (likelihood.isStochastic())
                return true;
        }
        
        return false;
    }
    
    @Override
    public double getNonStochasticLogP() {
        double logP = 0;
        for (TreeLikelihood likelihood : pLikelihoods.get()) {
            logP += likelihood.getNonStochasticLogP();
            if (Double.isInfinite(logP) || Double.isNaN(logP)) {
                return logP;
            }
        }
        return logP;
    }
    
} // class MixedTreeLikelihood
