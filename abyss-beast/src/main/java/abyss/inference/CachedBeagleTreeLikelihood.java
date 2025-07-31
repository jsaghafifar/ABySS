package abyss.inference;

import beast.base.evolution.likelihood.BeagleTreeLikelihood;

public class CachedBeagleTreeLikelihood extends BeagleTreeLikelihood {

    protected double[] storedPatternLogLikelihoods;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        if (patternLogLikelihoods == null) patternLogLikelihoods = new double[patternCount];
        beagle.getSiteLogLikelihoods(patternLogLikelihoods);
        storedPatternLogLikelihoods = new double[alignment.getPatternCount()];
    }

    @Override
    public void store() {
        super.store();
        System.arraycopy(patternLogLikelihoods, 0, storedPatternLogLikelihoods, 0, patternLogLikelihoods.length);
    }

    @Override
    public void restore() {
        super.restore();
        double[] tmp = patternLogLikelihoods;
        patternLogLikelihoods = storedPatternLogLikelihoods;
        storedPatternLogLikelihoods = tmp;
    }

}
